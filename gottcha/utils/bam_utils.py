#!/usr/bin/env python3
"""
bam_cov_mismatch.py

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
import math
import multiprocessing as mp
import os
import sys
import logging
import pandas as pd
from typing import Iterable, List, Optional, Tuple

import numpy as np
import pysam

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


def _process_chunk(task: Tuple[str, int, int]) -> Tuple[str, int, int, int, int, int, int, float]:
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

    min_mapq = _CFG["min_mapq"]
    min_frac = _CFG["min_frac"]
    min_idt = _CFG["min_idt"]
    inc_sec = _CFG["include_secondary"]
    inc_sup = _CFG["include_supplementary"]
    inc_dup = _CFG["include_duplicates"]
    inc_qcf = _CFG["include_qcfail"]
    min_alen = _CFG["min_alen"]
    split_read_flag = _CFG["split_read_flag"]

    aln_starts_in_chunk_flag = False
    numreads = 0
    readlength = 0
    indels = 0
    invalid_alns = 0
    bam = _BAM

    # Iterate reads overlapping this region.
    for aln in bam.fetch(rname, start0, end0):
        # Basic filters
        if aln.is_unmapped:
            continue
        if (not inc_sec) and aln.is_secondary:
            continue
        if (not inc_sup) and aln.is_supplementary:
            continue
        if (not inc_dup) and aln.is_duplicate:
            continue
        if (not inc_qcf) and aln.is_qcfail:
            continue
        if aln.mapping_quality < min_mapq:
            continue

        if aln.reference_start >= start0:
            aln_starts_in_chunk_flag = True

        if min_idt > 0.0 and aln.has_tag('NM'):
            mm_idt = aln.get_tag('NM') / aln.alen
            if min_idt > (1-mm_idt):
                if aln_starts_in_chunk_flag: invalid_alns += 1
                continue

        if min_frac > 0.0:
            if (aln.alen / aln.query_length) < min_frac or (aln.alen / aln.reference_length) < min_frac:
                if aln_starts_in_chunk_flag: invalid_alns += 1
                continue

        if min_alen > 0 and aln.alen < min_alen:
            if aln_starts_in_chunk_flag: invalid_alns += 1
            continue

        # Note: aln.reference_start is 0-based leftmost coordinate of the alignment on the reference.
        # Only count reads that have their aligned portion starting within the chunk towards numreads, to avoid double-counting reads that span multiple chunks.
        if aln_starts_in_chunk_flag:
            # If split_read_flag is set, only count reads with ZC tag (the first chunked reads) towards numreads.
            if split_read_flag:
                if aln.has_tag('ZC'):
                    numreads += 1
            else:
                numreads += 1
            
            # count total read length (including softclips) for mean depth calculation
            readlength += aln.query_length

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
    mapped_bases = int(depth.sum()) # total aligned bases (including matches and mismatches)

    # Positions where mismatch fraction > 0.5 among reads with aligned bases
    # i.e. mm / depth > 0.5  ->  2*mm > depth
    consensus_diff = int(np.count_nonzero((depth > 0) & (mm * 2 > depth)))

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

def comp_cov_mismatches(bam_path: str,
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
            result[1] +=  1 # start0
            ref_chunk_results.append(result)

    finally:
        pool.close()
        pool.join()

    return ref_chunk_results


def pile_lvl_zscore(tol_bp, tol_sig_len, linear_len):
    """
    Calculate Z-score for the depth of coverage of mapped regions.
    
    This determines how unusual the coverage depth is compared to expected depth
    based on a statistical model. Higher Z-scores may indicate biased mapping.
    
    Parameters:
        tol_bp (int): Total number of mapped bases
        tol_sig_len (int): Total length of the signature
        linear_len (int): Linear length (de-duplicated) covered by mappings
        
    Returns:
        float: Z-score for the depth distribution (or 0 if calculation fails)
    """
    try:
        avg_doc = tol_bp/tol_sig_len
        lin_doc = tol_bp/linear_len
        v = (linear_len*(lin_doc-avg_doc)**2 + (tol_sig_len-linear_len)*(avg_doc)**2)/tol_sig_len
        sd = math.sqrt(v)
        if sd == 0.0:
            return 0
        else:
            return (lin_doc-avg_doc)/sd
    except:
        return 0


def group_refs_to_strains(ref_chunk_results, acc_list, acc_list_action, df_stats):
    """
    Group reference mapping results by strains and calculate strain-level statistics.
    
    Converts the mapping results dictionary to a pandas DataFrame and groups by
    taxonomic identifier. Calculates various statistics including total mapped bases,
    read counts, coverage, and depth of coverage.
    
    Parameters:
        ref_chunk_results (list): List of mapping statistics for each reference fragment chunks (output from comp_cov_mismatches)
        acc_list (list, optional): List of accessions of interest
        acc_list_action (str, optional): Action to take with the accession list (e.g., "exclude")
        df_stats (pandas.DataFrame): DataFrame containing genome signature statistics
    Returns:
        pandas.DataFrame: DataFrame with strain-level statistics
    """
    # covert mapping info to df
    r_chunk_df = pd.DataFrame(ref_chunk_results[1:], columns=ref_chunk_results[0])

    # retrieve sig fragment info
    r_chunk_df['RNAME'] = r_chunk_df['RNAME'].str.rstrip('|')

    r_df = r_chunk_df.groupby('RNAME').agg({
        'COVBASES':'sum', # of covered signature bases
        'NUMREADS':'sum', # of mapped reads
        'MISMATCHES':'sum', # of mismatches
        'INDELS':'sum', # of indels
        'MAPPED_BASES':'sum', # total length of mapped bases (including matches and mismatches)
        'CONSENSUS_DIFF':'sum', # number of positions with >50% mismatches among aligned reads
        'INVALID_ALNS':'sum', # total invalid alignments (after filters) for this reference
        'READLENGTH':'sum', # total length of reads
    }).reset_index()

    # add reportable read count
    r_df['AOI_READ_COUNT'] = 0
    aoi_read_count = 0

    r_df[['ACC','RSTART','REND','TAXID']] = r_df['RNAME'].str.split('|', expand=True)

    if acc_list:
        idx = (r_df['ACC'].isin(acc_list) | r_df['RNAME'].isin(acc_list))
        r_df.loc[idx, 'AOI_READ_COUNT'] = r_df.loc[idx, 'NUMREADS'] # report the read count for the accession#s of interest
        aoi_read_count = r_df.loc[idx, 'NUMREADS'].sum()

        if acc_list_action == 'filter_out':
            r_df = r_df.loc[~idx] # set mapped bases, read count, mismatch and covered sig len to 0 for the accession#s of interest
        elif acc_list_action == 'filter_in':
            r_df = r_df[idx].reset_index(drop=True)

        # if after applying the accession list filter, there is no valid mapping left, exit the program
        if len(r_df) == 0:
            logging.info(f"No valid mappings after applying accession list filter. Exiting.")
            sys.exit(0)

    r_df['RSTART'] = r_df['RSTART'].astype(int)
    r_df['REND'] = r_df['REND'].astype(int)
    r_df['SLEN'] = r_df['REND']-r_df['RSTART']+1 # length of the signature fragment

    # group by strain
    str_df = r_df.groupby(['TAXID']).agg({
        'COVBASES':'sum', # of covered signature bases
        'NUMREADS':'sum', # of mapped reads
        'MISMATCHES':'sum', # of mismatches
        'INDELS':'sum', # of indels
        'MAPPED_BASES':'sum', # total length of mapped bases (including matches and mismatches)
        'CONSENSUS_DIFF':'sum', # number of positions with >50% mismatches among aligned reads
        'INVALID_ALNS':'sum', # total invalid alignments (after filters) for this reference
        'READLENGTH':'sum', # total length of reads
        'SLEN':'sum', # length of this signature fragments (mapped)
        'AOI_READ_COUNT':'sum'  # reportable read count
    }).reset_index()
    # total length of signatures
    str_df['TOTAL_SIG_LEN'] = str_df['TAXID'].map(df_stats['TotalLength'])
    str_df['BEST_SIG_COV'] = str_df['COVBASES']/str_df['TOTAL_SIG_LEN'] # bLC:  best linear coverage of a strain
    str_df['DEPTH'] = str_df['MAPPED_BASES']/str_df['TOTAL_SIG_LEN'] # roll-up DoC
    str_df['NOTE'] = str_df['TAXID'].map(df_stats['Note']).fillna('') # note for the strain
    
    # rename columns
    str_df.rename(columns={
        "MAPPED_BASES": "TOTAL_BP_MAPPED",
        "NUMREADS":     "READ_COUNT",
        "MISMATCHES":   "TOTAL_BP_MISMATCH",
        "INDELS":       "TOTAL_BP_INDEL",
        "READLENGTH":   "TOTAL_READ_LEN",
        "COVBASES":     "COVERED_SIG_LEN",
        "SLEN":         "MAPPED_SIG_LEN",
    }, inplace=True)

    # check if TOTAL_SIG_LEN is 0, report the TAXID and exit
    # this should not happen if the database and corresponding stats file are correct
    if str_df['TOTAL_SIG_LEN'].eq(0).any():
        logging.fatal(f"Error: total signature length is ZERO for some mapped strains. Please check your database.")
        sys.exit(1)

    # get genome size
    str_df['SIG_LEVEL'] = str_df['TAXID'].map(df_stats['DB_level'])
    str_df['GENOME_SIZE'] = str_df['TAXID'].map(df_stats['GenomeSize'])
    str_df['GENOME_COUNT'] = 1

    # infer total genome contents
    str_df['GENOMIC_CONTENT_EST'] = str_df['TOTAL_BP_MAPPED']/str_df['TOTAL_SIG_LEN']*str_df['GENOME_SIZE']

    # estimate z-score
    str_df['ZSCORE'] = str_df.apply(lambda x: pile_lvl_zscore(x.TOTAL_BP_MAPPED, x.TOTAL_SIG_LEN, x.COVERED_SIG_LEN), axis=1)

    return str_df, aoi_read_count


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

    # Coordinate output style
    p.add_argument(
        "--imap-chunksize",
        type=int,
        default=1,
        help="chunksize passed to multiprocessing imap/imap_unordered (default: 1).",
    )

    args = p.parse_args(argv)

    ref_results = comp_cov_mismatches(
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

    return 0

if __name__ == "__main__":
    raise SystemExit(main())
