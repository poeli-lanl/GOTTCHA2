#!/usr/bin/env python3
from __future__ import annotations

import multiprocessing as mp
import os
import sys
import logging
import pandas as pd
from typing import Iterable, List, Optional, Tuple

import pysam
from . import taxonomy

# Global BAM handle and config for worker processes
taxa_dict = {}
lineage_cache = {}  # Cache for reference taxid to qualified taxids mapping

def parse_taxids(taxid_arg: str, res_df: pd.DataFrame, full_tsv_fn: str) -> Tuple[dict, list]:
    """Parse taxids from command line arg or file"""

    taxa_list = []
    qualified_taxa = pd.DataFrame()
    taxa_df = pd.DataFrame()

    if res_df.shape[1] > 0:
        taxa_df = res_df
        logging.info(f"Successfully loaded taxonomy profile with {len(taxa_df)} entries")
    else:
        try:
            logging.info(f"Reading taxonomy file {full_tsv_fn}...")
            taxa_df = pd.read_csv(full_tsv_fn,
                                sep='\t',
                                engine='python',
                                quoting=3,
                                on_bad_lines='skip',
                                dtype={'NOTE': str})

            logging.info(f"Successfully loaded taxonomy profile with {len(taxa_df)} entries")
        except Exception as e:
            logging.error(f"Error reading taxonomy file: {e}")
            sys.exit(1)

    # Filter in entries specified by taxid_arg
    filtered_idx = None

    if 'NOTE' in taxa_df.columns:
        filtered_idx = ~taxa_df['NOTE'].str.contains('Filtered out', na=False)

    if taxid_arg and taxid_arg != 'all':
        if taxid_arg.startswith('@'):
            # Read taxids from file
            filename = taxid_arg[1:]  # Remove @ prefix
            try:
                with open(filename) as f:
                    taxa_list = [x.strip() for x in f.readlines() if x.strip() and not x.startswith('#')]
            except IOError as e:
                logging.error(f"Error reading taxid file {filename}: {e}")
                sys.exit(1)
        else:
            # Parse comma-separated list
            taxa_list = [x.strip() for x in taxid_arg.split(',')]

        if taxa_list:
            filtered_idx &= (taxa_df['TAXID'].isin(taxa_list) | taxa_df['NAME'].isin(taxa_list))

    # Filter out entries with "Filtered out" notes
    if filtered_idx is not None:
        qualified_taxa = taxa_df[filtered_idx]
        logging.info(f"Found {len(qualified_taxa)} qualified taxa after filtering")
    else:
        qualified_taxa = taxa_df

    # Ensure these columns exist
    if not all(col in qualified_taxa.columns for col in ['LEVEL', 'NAME', 'TAXID']):
        logging.error(f"Required columns missing in taxonomy file. Available columns: {qualified_taxa.columns.tolist()}")
        sys.exit(1)

    # Pre-compute a mapping from reference taxids to qualified taxids
    # This avoids expensive lineage lookups during processing
    logging.info("Building taxonomy lookup index...")

    taxa_dict = {}

    # Gather all qualified taxids
    qualified_taxids = []
    for _, row in qualified_taxa[['LEVEL', 'NAME', 'TAXID']].iterrows():
        if pd.notna(row['TAXID']):
            taxid = str(row['TAXID']).strip()
            qualified_taxids.append(taxid)
            taxa_dict[taxid] = {
                'level': str(row['LEVEL']).replace(' ', '_') if pd.notna(row['LEVEL']) else 'unknown',
                'name': str(row['NAME']).replace(' ', '_') if pd.notna(row['NAME']) else 'unknown'
            }

    logging.debug(f"Qualified taxa:\n{qualified_taxa}")
    logging.debug(f"Qualified taxids: {qualified_taxids}")

    return taxa_dict, qualified_taxids


def _init_worker(bam_path: str,
                 taxa_dict: dict,
                 format: str,
                 min_frac: float,
                 min_idt: float,
                 min_alen: int,
                 min_mapq: int = 0,
                 htslib_threads: int = 1,
                 include_secondary: bool = False,
                 include_supplementary: bool = False,
                 include_duplicates: bool = False,
                 include_qcfail: bool = False) -> None:
    """Initializer for each worker process: open BAM once and stash filters."""
    global _BAM, _CFG
    _BAM = pysam.AlignmentFile(bam_path, "rb", threads=htslib_threads)
    _CFG = {
        "min_mapq": min_mapq,
        "taxa_dict": taxa_dict,
        "min_frac": min_frac,
        "min_idt": min_idt,
        "min_alen": min_alen,
        "include_secondary": include_secondary,
        "include_supplementary": include_supplementary,
        "include_duplicates": include_duplicates,
        "include_qcfail": include_qcfail,
        "format": format
    }


def _iter_tasks(references: List[str],
                max_per_taxon: int,
                acc_list: list, 
                acc_list_action: str
                ) -> Iterable[Tuple[str, int, int]]:
    
    global lineage_cache

    for ref in references:
        # Extract taxid from reference
        try:
            acc, rstart, rend, ref_taxid, _ = ref.split('|')
        except ValueError:
            logging.debug(f"Malformed reference: {ref}")
            continue  # Skip malformed references

        # Skip if accession is in the exclusion list (if applicable)
        aoi_flag = False
        if acc_list:
            if (acc in acc_list) or (ref in acc_list):
                aoi_flag = True
                if acc_list_action == 'filter_out':
                    continue
            else:
                if acc_list_action == 'filter_in':
                    continue

        if len(lineage_cache[ref_taxid]) > 0:
            yield (ref, max_per_taxon, aoi_flag)


def _extract_worker(task: Tuple[str, int, bool]) -> dict:
    """
    Process one (rname, start0, end0) chunk.

    Returns:
      (rname, start0, end0, numreads, covbases, mismatches_total,
       consensus_diff, mean_depth)
    """
    global _BAM, _CFG, lineage_cache
    assert _BAM is not None, "Worker BAM handle not initialized"
    
    taxon_seqs = {}  # Dictionary to hold sequences for each taxid
    ref, max_per_taxon, aoi_flag = task
    
    try:
        acc, rstart, rend, ref_taxid, _ = ref.split('|')
    except ValueError:
        logging.debug(f"Malformed reference: {ref}")

    bam = _BAM

    min_mapq = _CFG["min_mapq"]
    taxa_dict = _CFG["taxa_dict"]
    format = _CFG["format"]
    min_frac = _CFG["min_frac"]
    min_idt = _CFG["min_idt"]
    min_alen = _CFG["min_alen"]
    inc_sec = _CFG["include_secondary"]
    inc_sup = _CFG["include_supplementary"]
    inc_dup = _CFG["include_duplicates"]
    inc_qcf = _CFG["include_qcfail"]

    # Iterate reads overlapping this region.
    for aln in bam.fetch(ref):
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

        if min_idt > 0.0 and aln.has_tag('NM'):
            mm_idt = aln.get_tag('NM') / aln.alen
            if min_idt > (1-mm_idt):
                continue

        if min_frac > 0.0:
            if (aln.alen / aln.query_length) < min_frac and (aln.alen / bam.get_reference_length(ref)) < min_frac:
                continue

        if min_alen > 0 and aln.alen < min_alen:
            continue

        matching_taxids = lineage_cache[ref_taxid]

        # Process the matching taxa
        for taxid in matching_taxids:
            # Initialize list for this taxid if needed
            if taxid not in taxon_seqs:
                taxon_seqs[taxid] = []

            # Only collect up to max_per_taxon sequences per taxon
            if (max_per_taxon==0) or (len(taxon_seqs[taxid]) < max_per_taxon):
                # Create FASTA entry with taxonomy information
                level = taxa_dict[taxid]['level']
                name = taxa_dict[taxid]['name']
                rname = aln.query_name
                region = (aln.reference_start+1, aln.reference_end)
                mapping_idt = aln.get_tag('NM') / aln.alen if aln.has_tag('NM') else 0
                mapping_frac = max((aln.alen / aln.query_length), (aln.alen / bam.get_reference_length(ref))) if bam.get_reference_length(ref) > 0 else 0

                # determine if the read is the first or second mate
                mate = ''
                if aln.is_paired:
                    mate = '.1' if aln.is_read1 else '.2'

                if format == 'fasta':
                    fasta_entry = f">{rname}{mate}|{ref}:{region[0]}..{region[1]} LEVEL={level} NAME={name} TAXID={taxid} AOI={aoi_flag} MG={aln.alen} MI={mapping_idt:.2f} MF={mapping_frac:.2f}\n{aln.query_sequence}\n"
                else:
                    fasta_entry = f"@{rname}{mate}|{ref}:{region[0]}..{region[1]} LEVEL={level} NAME={name} TAXID={taxid} AOI={aoi_flag} MG={aln.alen} MI={mapping_idt:.2f} MF={mapping_frac:.2f}\n{aln.query_sequence}\n+\n{aln.query_qualities_str}\n"
                
                taxon_seqs[taxid].append(fasta_entry)

    return taxon_seqs


def extract_sequences_by_taxonomy(bam_path: str,
                                  taxa_dict: dict,
                                  qualified_taxids: list,
                                  o,
                                  numthreads: int,
                                  matchFraction: float,
                                  matchIdentity: float,
                                  matchLength: int,
                                  max_per_taxon: int,
                                  acc_list: list,
                                  acc_list_action: str,
                                  format: str = 'fasta'):
    """
    Extract sequences mapping to taxa from the full taxonomy report.

    For each taxon in the full report, extract up to max_per_taxon sequences.

    Parameters:
        bam_path (str): Path to the BAM file
        taxa_dict (dict): Dictionary containing taxonomy information
        qualified_taxids (list): List of qualified taxonomic IDs to extract
        o (file): Output file handle for the extracted sequences
        numthreads (int): Number of threads to use for processing
        matchFraction (float): Minimum fraction required for a valid match
        matchIdentity (float): Minimum identity required for a valid match
        matchLength (int): Minimum length required for a valid match
        max_per_taxon (int): Maximum number of sequences to extract per taxon; 0 is unlimited.
        format (str): Output format ('fasta' or 'fastq')

    Returns:
        tuple: (taxon_count, seq_count) - Number of taxa and total sequences extracted
    """

    global lineage_cache

    logging.debug(f"taxa_dict: {taxa_dict}...")
    logging.debug(f"qualified_taxids: {qualified_taxids}...")

    if not os.path.exists(bam_path):
        logging.fatal(f"ERROR: BAM not found: {bam_path}")
        return 2

    # Open BAM in main process to validate index and obtain reference lengths
    try:
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            if not bam.has_index():
                logging.fatal("ERROR: BAM index (.bai) not found or not readable. Pysam requires an index.")
                return 2
            references = list(bam.references)
    except Exception as e:
        logging.fatal(f"ERROR: Failed to open BAM: {e}")
        return 2

    # Pre-compute a mapping from reference taxids to qualified taxids to avoid expensive lineage lookups during processing
    for ref in references:
        # Extract taxid from reference
        try:
            ref_taxid = ref.split('|')[3]
        except ValueError:
            logging.debug(f"Malformed reference: {ref}")
            continue  # Skip malformed references

        # Check if we already know what qualified taxa this reference belongs to
        if ref_taxid in lineage_cache:
            matching_taxids = lineage_cache[ref_taxid]
        else:
            # If not, find all qualified taxa this reference belongs to
            matching_taxids = []
            ref_lineage = None

            for q_taxid in qualified_taxids:
                # Avoid recomputing the lineage for each taxid check
                if ref_lineage is None:
                    ref_lineage = taxonomy.taxid2fullLineage(ref_taxid, space2underscore=False)

                if f"|{q_taxid}|" in ref_lineage:
                    matching_taxids.append(q_taxid)

            # Cache the result
            lineage_cache[ref_taxid] = matching_taxids

    # Generate tasks for worker processes
    tasks = _iter_tasks(references, max_per_taxon, acc_list, acc_list_action)
    chunk_results = {}
    # Merge results from this chunk
    all_taxon_seqs = {}

    # (default is one-based if neither flag used; argparse sets one_based True by default)
    # endpos will be end0 in both conventions; interpretation differs.

    pool = mp.Pool(
        processes=numthreads,
        initializer=_init_worker,
        initargs=(
            bam_path,
            taxa_dict,
            format,
            matchFraction,
            matchIdentity,
            matchLength,
        ),
    )

    logging.info(f"Starting extraction for {len(qualified_taxids)} qualified taxa...")

    try:
        mapper = pool.imap_unordered
        for result in mapper(_extract_worker, tasks):
            for key, value in result.items():
                chunk_results.setdefault(key, []).extend(value)
    
        for taxid, seqs in chunk_results.items():
            logging.info(f"Extracted {len(seqs)} sequences for taxid {taxid}")

            if taxid not in all_taxon_seqs:
                all_taxon_seqs[taxid] = []

            # Add sequences, respecting the max_per_taxon limit
            if max_per_taxon > 0:
                remaining = max_per_taxon - len(all_taxon_seqs[taxid])

                if remaining > 0:
                    all_taxon_seqs[taxid].extend(seqs[:remaining])
            else:
                all_taxon_seqs[taxid].extend(seqs)
    finally:
        pool.close()
        pool.join()


    # Write sequences to output file
    logging.info("Writing sequences to output file...")

    total_seqs = 0
    taxon_count = 0

    for taxid, seqs in all_taxon_seqs.items():
        if seqs:  # If we got any sequences for this taxon
            taxon_count += 1
            total_seqs += len(seqs)
            
            output_seqs = "".join(seqs)
            o.write(output_seqs)

    return taxon_count, total_seqs