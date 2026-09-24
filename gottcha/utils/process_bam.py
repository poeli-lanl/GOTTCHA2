#!/usr/bin/env python3
"""
process_bam.py

Compute per-region coverage and consensus-mismatch metrics from a BAM file
(with .bai index present), **without** a reference FASTA.

Assumptions / notes:
- Alignments are from minimap2 with `--eqx`, so mismatches are encoded as CIGAR op `X`
  and matches as `=`.
- Depth is computed from aligned query bases (CIGAR ops M/= /X). Deletions (D) and
  refskips (N) do not contribute to depth in this implementation.
- "mismatches" counts total mismatched aligned bases across all reads (sum of `X`).
- "pileup_mismatch" counts reference positions where >50% of aligned reads at that
  position are mismatches (i.e., #positions where X_depth / depth > 0.5).

Parallelization:
- References are split into fixed-size chunks along their length. Each chunk is processed
  independently in a worker process.
- Each worker opens the BAM once (via Pool initializer) for performance.

Output columns (TSV):
- rname
- startpos
- endpos
- numreads
- covbases
- coverage
- mismatches
- pileup_mismatch
- meandepth
"""

from __future__ import annotations

import argparse
import multiprocessing as mp
import os
import sys
import logging
from collections import Counter, defaultdict
from typing import Iterable, List, Optional, Tuple

from pathlib import Path
import numpy as np
import pysam

from gottcha.utils import extract_reads

# Global BAM handle and config for worker processes
_BAM: Optional[pysam.AlignmentFile] = None
_CFG = {}

def _init_worker(
    bam_path: str,
    htslib_threads: int,
    min_mapq: int,
    min_frac: float,
    min_idt: float,
    min_alen: int,
    include_secondary: bool,
    include_supplementary: bool,
    include_duplicates: bool,
    include_qcfail: bool,
    split_read_flag: Optional[bool] = False,
) -> None:
    """Initializer for each worker process: open BAM once and stash filters."""
    global _BAM, _CFG
    _BAM = pysam.AlignmentFile(bam_path, "rb", threads=htslib_threads)
    _CFG = {
        "min_mapq": min_mapq,
        "min_frac": min_frac,
        "min_idt": min_idt,
        "min_alen": min_alen,
        "include_secondary": include_secondary,
        "include_supplementary": include_supplementary,
        "include_duplicates": include_duplicates,
        "include_qcfail": include_qcfail,
        "split_read_flag": split_read_flag,
    }


def _is_aln_valid(aln, rname) -> Tuple[bool, Optional[str]]:
    global _CFG, _BAM
    min_mapq = _CFG["min_mapq"]
    min_frac = _CFG["min_frac"]
    min_idt = _CFG["min_idt"]
    inc_sec = _CFG["include_secondary"]
    inc_sup = _CFG["include_supplementary"]
    inc_dup = _CFG["include_duplicates"]
    inc_qcf = _CFG["include_qcfail"]
    min_alen = _CFG["min_alen"]
    bam = _BAM

    if aln.is_unmapped:
        return False, "aln_type"
    if (not inc_sec) and aln.is_secondary:
        return False, "aln_type"
    if (not inc_sup) and aln.is_supplementary:
        return False, "aln_type"
    if (not inc_dup) and aln.is_duplicate:
        return False, "aln_type"
    if (not inc_qcf) and aln.is_qcfail:
        return False, "aln_type"
    if aln.mapping_quality < min_mapq:
        return False, "aln_quality"

    # Note: aln.reference_start is 0-based leftmost coordinate of the alignment on the reference.
    # Only count reads that have their aligned portion starting within the chunk towards numreads, to avoid double-counting reads that span multiple chunks.
    if min_idt > 0.0 and aln.has_tag('de'):
        mm_idt = 1-aln.get_tag('de')
        if min_idt > mm_idt:
            return False, "aln_quality"

    if min_frac > 0.0:
        if aln.query_length <= 0:
            #for hard clips, query_length can be 0, recover it from CIGAR
            query_length = sum(length for op, length in aln.cigartuples if op in (0, 1, 4, 5, 7, 8)) # M/I/S/H/=/X
        else:
            query_length = aln.query_length

        if (aln.alen / query_length) < min_frac and (aln.alen / bam.get_reference_length(rname)) < min_frac:
            return False, "aln_quality"

    if min_alen > 0 and aln.alen < min_alen:
        return False, "aln_quality"

    return True, None

def _is_pair_owner(aln):
    """Return True for exactly one alignment of a mapped pair."""
    if not aln.is_paired:
        return True

    if aln.mate_is_unmapped:
        return True

    # Canonical location for current alignment
    this_key = (
        aln.reference_id,
        aln.reference_start,
        0 if aln.is_read1 else 1,
    )

    # Canonical location for mate
    mate_key = (
        aln.next_reference_id,
        aln.next_reference_start,
        1 if aln.is_read1 else 0,
    )

    return this_key < mate_key


def _process_chunk(task: Tuple[str, int, int]) -> List:
    """
    Process one (rname, start0, end0) chunk.

    Returns:
      (rname, start0, end0, numreads, covbases, mismatches_total,
       consensus_diff, mean_depth)
    """
    global _BAM, _CFG
    assert _BAM is not None, "Worker BAM handle not initialized"

    rname, start0, end0 = task
    L = end0 - start0
    if L <= 0:
        return [rname, start0, end0, 0,0,0,0,0,0,0,0]

    # Difference arrays (signed) so we can do O(segments) updates and O(L) cumsums.
    # depth[pos] = #reads with an aligned base at that position (from CIGAR ops M/= /X)
    # mm[pos] = #reads with CIGAR X at that position
    depth_diff = np.zeros(L + 1, dtype=np.int32)
    mm_diff = np.zeros(L + 1, dtype=np.int32)
    split_read_flag = _CFG["split_read_flag"]

    numreads = 0
    readlength = 0
    indels = 0
    invalid_alns = 0
    bam = _BAM

    # Iterate reads overlapping this region.
    for aln in bam.fetch(rname, start0, end0):
        valid, reason = _is_aln_valid(aln, rname)

        if aln.reference_start < start0:
            valid, reason = False, "aln_position"
        
        if not valid:
            invalid_alns += 1
            continue

        # Count valid reads
        # If split_read_flag is set, only count reads with ZC tag (the first chunked reads) towards numreads.
        if split_read_flag:
            if aln.has_tag('ZC'):
                numreads += 1
        else:
            if aln.is_secondary or aln.is_supplementary:
                pass
            elif _is_pair_owner(aln):
                numreads += 1
        
        # count total read length (including softclips) for mean depth calculation
        readlength += aln.alen

        cig = aln.cigartuples
        if not cig:
            continue

        ref_pos = aln.reference_start
        block_start: Optional[int] = None

        # CIGAR operation codes in pysam:
        # 0=M, 1=I, 2=D, 3=N, 4=S, 5=H, 6=P, 7==, 8=X
        for op, length in cig:
            if length <= 0:
                continue

            if op in (0, 7, 8):  # aligned query bases consuming reference
                if block_start is None:
                    block_start = ref_pos

                if op == 8:  # X mismatches
                    seg_s = ref_pos
                    seg_e = ref_pos + length  # exclusive
                    if seg_e > start0 and seg_s < end0:
                        if seg_s < start0:
                            seg_s = start0
                        if seg_e > end0:
                            seg_e = end0
                        mm_diff[seg_s - start0] += 1
                        mm_diff[seg_e - start0] -= 1

                ref_pos += length

            elif op in (2, 3):  # D or N: consumes reference but not query -> breaks aligned block
                if block_start is not None:
                    seg_s = block_start
                    seg_e = ref_pos  # exclusive
                    if seg_e > start0 and seg_s < end0:
                        if seg_s < start0:
                            seg_s = start0
                        if seg_e > end0:
                            seg_e = end0
                        depth_diff[seg_s - start0] += 1
                        depth_diff[seg_e - start0] -= 1
                        indels += length
                    block_start = None

                ref_pos += length

            else:
                # I/S/H/P: does not consume reference; does not affect ref_pos.
                # We do not break the block on insertions/softclips because reference
                # positions remain contiguous.
                continue

        # Close any remaining aligned block
        if block_start is not None:
            seg_s = block_start
            seg_e = ref_pos
            if seg_e > start0 and seg_s < end0:
                if seg_s < start0:
                    seg_s = start0
                if seg_e > end0:
                    seg_e = end0
                depth_diff[seg_s - start0] += 1
                depth_diff[seg_e - start0] -= 1

    # Build per-base depth and mismatch arrays
    depth = np.cumsum(depth_diff[:-1], dtype=np.int32)  # length L
    mm = np.cumsum(mm_diff[:-1], dtype=np.int32)        # length L

    covbases = int(np.count_nonzero(depth))
    mismatches_total = int(mm.sum())
    mapped_bases = int(depth.sum())+indels # total aligned bases (including matches and mismatches)

    # Positions where mismatch fraction > 0.5 among reads with aligned bases
    # i.e. mm / depth > 0.5  ->  2*mm > depth
    consensus_diff = int(np.count_nonzero((depth > 0) & (mm * 2 > depth)))

    logging.debug(f"Processed {rname}: {numreads} reads, {covbases} covbases, {mismatches_total} mismatches, {consensus_diff} consensus_diff, {indels} indels, {mapped_bases} mapped_bases, {invalid_alns} invalid_alns")

    return [rname,
            start0,
            end0,
            numreads,
            covbases,
            mismatches_total,
            indels,
            consensus_diff,
            mapped_bases,
            invalid_alns,
            readlength]


def _iter_tasks(references: List[str], lengths: List[int], chunk_size: int) -> Iterable[Tuple[str, int, int]]:
    for rname, rlen in zip(references, lengths):
        if rlen <= 0:
            continue
        cs = rlen if chunk_size <= 0 else chunk_size
        for start0 in range(0, rlen, cs):
            end0 = min(start0 + cs, rlen)
            yield (rname, start0, end0)
def _taxid_from_rname(rname: str) -> str:
    """Extract TAXID from reference name: ...|...|TAXID|... ."""
    parts = rname.rsplit("|", 2)
    return parts[-2] if len(parts) >= 3 else rname


def _read_key(aln) -> Tuple[str, int]:
    """Distinguish mates while keeping ordinary/long reads simple."""
    if aln.is_paired:
        if aln.is_read1:
            return (aln.query_name, 1)
        if aln.is_read2:
            return (aln.query_name, 2)
    return (aln.query_name, 0)


def write_taxid_network(bam_path: str, node_path: str, edge_path: str) -> None:
    """
    Write TAXID node/edge statistics from primary and supplementary alignments.

    Edges are generated from:

        1. Primary -> supplementary alignment of the same read/mate.

        2. For paired-end reads:
               read1 primary TAXID -> read2 primary TAXID

           when the two mates map to different TAXIDs.

    Secondary alignments are ignored.
    """
    

    node_reads = defaultdict(set)
    node_counts = defaultdict(lambda: [0, 0, 0])

    # Primary alignment for each individual read/mate:
    #
    #   (QNAME, 0) -> TAXID     single/long read
    #   (QNAME, 1) -> TAXID     mate 1
    #   (QNAME, 2) -> TAXID     mate 2
    primary_by_read = {}

    # Supplementary TAXIDs associated with each read/mate.
    supp_by_read = defaultdict(list)

    # Primary mappings for paired-end reads:
    #
    #   QNAME -> {1: taxid1, 2: taxid2}
    pair_primary = defaultdict(dict)

    global _BAM, _CFG

    if not _CFG:
        logfile_prev = Path(bam_path).with_suffix(".log")
        if logfile_prev.is_file():
            (mi, mf, mg, sni_argv) = extract_reads.load_criteria_from_log(logfile_prev)
        
        _CFG = {
            "min_mapq": 0,
            "min_idt": mi if mi is not None else 0.85,
            "min_frac": mf if mf is not None else 0,
            "min_alen": mg if mg is not None else 100,
            "include_secondary": False,
            "include_supplementary": True,
            "include_duplicates": False,
            "include_qcfail": False,
            "split_read_flag": False,
        }

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        _BAM = bam

        for aln in bam.fetch(until_eof=True):

            # Avoid trying to obtain an RNAME for an unmapped alignment.
            if aln.is_unmapped or aln.reference_id < 0:
                continue

            # Do not use secondary alignments in this network.
            if aln.is_secondary:
                continue

            rname = bam.get_reference_name(aln.reference_id)
            taxid = _taxid_from_rname(rname)

            valid, reason = _is_aln_valid(aln, rname)
            if not valid:
                continue

            key = _read_key(aln)

            # Node statistics.
            #
            # This preserves the original meaning of READ_COUNT as
            # unique QNAMEs/fragments.
            node_reads[taxid].add(aln.query_name)
            node_counts[taxid][0] += 1  # all valid alignments

            if aln.is_supplementary:
                node_counts[taxid][2] += 1
                supp_by_read[key].append(taxid)

            else:
                # Since secondary alignments have already been excluded,
                # this is a true primary alignment.
                node_counts[taxid][1] += 1

                primary_by_read[key] = taxid

                # Remember mate-1 and mate-2 primary mappings separately.
                if aln.is_paired:
                    if aln.is_read1:
                        pair_primary[aln.query_name][1] = taxid
                    elif aln.is_read2:
                        pair_primary[aln.query_name][2] = taxid

    # ------------------------------------------------------------------
    # Build edges
    # ------------------------------------------------------------------

    edge_alignment_count = Counter()
    edge_reads = defaultdict(set)

    def add_edge(source_taxid, target_taxid, qname):
        if source_taxid == target_taxid:
            return

        edge = (source_taxid, target_taxid)
        edge_alignment_count[edge] += 1
        edge_reads[edge].add(qname)

    # 1. Primary -> supplementary edges.
    #
    # Supplementary alignments remain associated with their own mate:
    #
    #   read1 primary -> read1 supplementary
    #   read2 primary -> read2 supplementary
    #
    for key, targets in supp_by_read.items():
        source_taxid = primary_by_read.get(key)

        if source_taxid is None:
            continue

        qname = key[0]

        for target_taxid in targets:
            add_edge(source_taxid, target_taxid, qname)

    # 2. Paired-end read1 -> read2 edges.
    #
    # Example:
    #
    #   readX/1 -> taxid1
    #   readX/2 -> taxid2
    #
    # produces:
    #
    #   taxid1 -> taxid2
    #
    for qname, mates in pair_primary.items():
        taxid1 = mates.get(1)
        taxid2 = mates.get(2)

        if taxid1 is None or taxid2 is None:
            continue

        add_edge(taxid1, taxid2, qname)

    # ------------------------------------------------------------------
    # Write node table
    # ------------------------------------------------------------------

    with open(node_path, "w", encoding="utf-8") as out:
        out.write(
            "TAXID\tREAD_COUNT\tALIGNMENT_COUNT\t"
            "PRIMARY_ALIGNMENT_COUNT\tSUPP_ALIGNMENT_COUNT\n"
        )

        for taxid in sorted(node_counts):
            aln_count, primary_count, supp_count = node_counts[taxid]

            out.write(
                f"{taxid}\t"
                f"{len(node_reads[taxid])}\t"
                f"{aln_count}\t"
                f"{primary_count}\t"
                f"{supp_count}\n"
            )

    logging.info(
        f"{len(node_counts)} nodes. Node file written to: {node_path}",
    )

    # ------------------------------------------------------------------
    # Write edge table
    # ------------------------------------------------------------------

    with open(edge_path, "w", encoding="utf-8") as out:
        out.write(
            "SOURCE_TAXID\tTARGET_TAXID\t"
            "READ_COUNT\tALIGNMENT_COUNT\n"
        )

        for edge, aln_count in sorted(
            edge_alignment_count.items(),
            key=lambda x: (
                -len(edge_reads[x[0]]),
                -x[1],
                x[0][0],
                x[0][1],
            ),
        ):
            source_taxid, target_taxid = edge

            out.write(
                f"{source_taxid}\t"
                f"{target_taxid}\t"
                f"{len(edge_reads[edge])}\t"
                f"{aln_count}\n"
            )

    logging.info(
        f"{len(edge_alignment_count)} edges. Edge file written to: {edge_path}",
    )

    return node_path, edge_path


def parse_aln_from_bam(bam_path: str,
                       processes: int, 
                       min_frac: float,
                       min_idt: float,
                       min_alen: int,
                       min_mapq: Optional[int] = 0,
                       htslib_threads: Optional[int] = 1,
                       chunk_size: Optional[int] = 10_000_000,
                       imap_chunksize: Optional[int] = 1,
                       include_secondary: Optional[bool] = False,
                       include_supplementary: Optional[bool] = False,
                       include_duplicates: Optional[bool] = False,
                       include_qcfail: Optional[bool] = False,
                       split_read_flag: Optional[bool] = False,
                       ) -> int:
    if not os.path.exists(bam_path):
        print(f"ERROR: BAM not found: {bam_path}", file=sys.stderr)
        return 2

    # Open BAM in main process to validate index and obtain reference lengths
    try:
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            if not bam.has_index():
                print("ERROR: BAM index (.bai) not found or not readable. Pysam requires an index.", file=sys.stderr)
                return 2
            references = list(bam.references)
            lengths = list(bam.lengths)
    except Exception as e:
        print(f"ERROR: Failed to open BAM: {e}", file=sys.stderr)
        return 2

    logging.debug(f"Parsing {len(references)} references with {processes} processes...")
    logging.debug(f"Parameters: min_mapq={min_mapq}, min_frac={min_frac}, min_idt={min_idt}, min_alen={min_alen}, include_secondary={include_secondary}, include_supplementary={include_supplementary}, include_duplicates={include_duplicates}, include_qcfail={include_qcfail}, split_read_flag={split_read_flag}")

    tasks = _iter_tasks(references, lengths, chunk_size)

    # (default is one-based if neither flag used; argparse sets one_based True by default)
    # endpos will be end0 in both conventions; interpretation differs.

    pool = mp.Pool(
        processes=processes,
        initializer=_init_worker,
        initargs=(
            bam_path,
            htslib_threads,
            min_mapq,
            min_frac,
            min_idt,
            min_alen,
            include_secondary,
            include_supplementary,
            include_duplicates,
            include_qcfail,
            split_read_flag,
        ),
    )

    try:
        ref_chunk_results = []
        header = [
            "RNAME",             # rname,
            "STARTPOS",          # start0,
            "ENDPOS",            # end0,
            "NUMREADS",          # numreads,
            "COVBASES",          # covbases,
            "MISMATCHES",        # mismatches_total,
            "INDELS",            # indels,
            "CONSENSUS_DIFF",    # consensus_diff,
            "MAPPED_BASES",      # mapped_bases,
            "INVALID_ALNS",      # invalid_alns,
            "READLENGTH"         # readlength
        ]

        ref_chunk_results.append(header)
        mapper = pool.imap_unordered
        for result in mapper(_process_chunk, tasks, chunksize=imap_chunksize):
            result[1] +=  1 # start0 to 1-based startpos for output (endpos remains end0, which is exclusive in both conventions)
            ref_chunk_results.append(result)

    finally:
        pool.close()
        pool.join()

    logging.debug(f"Total signature fragments processed: {len(ref_chunk_results)-1}")

    return ref_chunk_results


def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(
        description="Compute coverage and consensus mismatch metrics from a BAM"
    )
    p.add_argument("bam", help="Input BAM path (requires .bai index).")
    p.add_argument("-o", "--out", required=True, help="Output TSV path.")
    p.add_argument(
        "-c",
        "--chunk-size",
        type=int,
        default=1_000_000,
        help="Chunk size in reference bases for parallel tasks (default: 1,000,000). Use 0 for whole-contig.",
    )
    p.add_argument(
        "-p",
        "--processes",
        type=int,
        default=max(1, mp.cpu_count() - 1),
        help="Worker processes (default: cpu_count-1).",
    )
    p.add_argument(
        "-t",
        "--htslib-threads",
        type=int,
        default=1,
        help="HTSlib threads per worker for BAM decompression (default: 1).",
    )

    # Filters
    p.add_argument("--min-mapq", type=int, default=0, help="Minimum MAPQ to keep an alignment (default: 0).")
    p.add_argument("--min-frac", type=float, default=0.0, help="Minimum fraction to keep an alignment (default: 0.0).")
    p.add_argument("--min-idt", type=float, default=0.0, help="Minimum identity to keep an alignment (default: 0.0).")
    p.add_argument("--min-alen", type=int, default=0, help="Minimum alignment length to keep an alignment (default: 0).")
    p.add_argument("--include-secondary", action="store_true", help="Include secondary alignments (default: off).")
    p.add_argument("--include-supplementary", action="store_true", help="Include supplementary alignments (default: off).")
    p.add_argument("--include-duplicates", action="store_true", help="Include duplicate-marked reads (default: off).")
    p.add_argument("--include-qcfail", action="store_true", help="Include QC-failed reads (default: off).")
    p.add_argument(
        "--taxid-network",
        action="store_true",
        help="Also write <out>.nodes.tsv and <out>.edges.tsv from primary/supplementary TAXID links (default: off).",
    )

    # Coordinate output style
    p.add_argument(
        "--imap-chunksize",
        type=int,
        default=1,
        help="chunksize passed to multiprocessing imap/imap_unordered (default: 1).",
    )

    args = p.parse_args(argv)

    ref_results = parse_aln_from_bam(
        bam_path=args.bam,
        processes=args.processes,
        min_frac=args.min_frac,
        min_idt=args.min_idt,
        min_alen=args.min_alen,
        min_mapq=args.min_mapq,
        htslib_threads=args.htslib_threads,
        chunk_size=args.chunk_size,
        imap_chunksize=args.imap_chunksize,
        include_secondary=args.include_secondary,
        include_supplementary=args.include_supplementary,
        include_duplicates=args.include_duplicates,
        include_qcfail=args.include_qcfail,
    )
    
    out_path=args.out
    with open(out_path, "w", encoding="utf-8") as out:
        for res in ref_results:
            out.write("\t".join(map(str, res)) + "\n")

    if args.taxid_network:
        base = out_path[:-4] if out_path.lower().endswith(".tsv") else out_path
        node_path = base + ".nodes.tsv"
        edge_path = base + ".edges.tsv"
        logging.debug(f"Writing TAXID network: {node_path}, {edge_path}")
        write_taxid_network(args.bam, node_path, edge_path)

    return 0

if __name__ == "__main__":
    raise SystemExit(main())
