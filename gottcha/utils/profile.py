#!/usr/bin/env python3

import re
import sys, os, time, subprocess
import pandas as pd
from pathlib import Path
import gc
from multiprocessing import set_start_method
import logging

if __package__ in (None, ''):
    # Keep direct execution routed through the package CLI.
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from gottcha.utils import (
    taxonomy,
    report,
    sam_to_bam,
    process_bam,
    ont_utils,
    read_mapping,
    aggregate_results,
    reciprocal_graph,
    extract_reads,
    prefilter,
    sig_archive,
)


def dependency_check(cmd: str) -> None:
    """
    Verify that external dependencies are available in the system.

    Attempts to execute the specified command with --help and checks if it runs
    successfully. Exits the program if the command is not found or fails.

    Parameters:
        cmd (str): Command to check

    Returns:
        None

    Raises:
        SystemExit: If the command is not found or fails
    """
    try:
        subprocess.check_call([cmd, "--help"], stdout=subprocess.DEVNULL)
    except Exception as e:
        sys.stderr.write(f"[ERROR] {cmd}: {e}\n")
        sys.exit(1)


def time_spend(start: float) -> str:
    """
    Calculate and format elapsed time since a given start time.

    Parameters:
        start (float): Starting time in seconds (as returned by time.time())

    Returns:
        str: Formatted time string in HH:MM:SS format
    """
    done = time.time()
    elapsed = done - start
    return time.strftime("%H:%M:%S", time.gmtime(elapsed))


def load_acc_list(filepath: str) -> set[str]:
    """
    Load a list of accession numbers to exclude from processing.

    Reads a file containing accession numbers (one per line) and returns them
    as a set for efficient lookup during processing. Empty lines are ignored.

    Parameters:
        filepath (str): Path to the file containing accession numbers to exclude

    Returns:
        set: Set of accession numbers to exclude. Returns empty set if input file is empty.

    Example:
        exclude_list = load_acc_list('exclude.txt')
    """
    with open(filepath) as f:
        acc_list = f.read().splitlines()

    if len(acc_list) == 0:
        logging.warning(f"Exclude accession list is empty.")
        return set()
    else:
        return set(acc_list)


def load_database_stats(db_stats_file: str) -> pd.DataFrame:
    """
    Load database signature statistics from a stats file.

    Reads a tab-delimited stats file containing information about
    taxonomic signatures and their lengths.

    Parameters:
        db_stats_file (str): Path to the database stats file

    Returns:
        pd.DataFrame: df indexed with taxid, contains signature lengths and genome sizes

    Note:
        The input stats file is an 9-column tab-delimited file with:
        1. Rank
        2. Name
        3. Taxid
        4. Superkingdom
        5. NumOfSeq
        6. Max
        7. Min
        8. TotalLength
        9. GenomeSize
       10. Note (optional)
    """

    # Determine the number of columns in the stats file
    header = pd.read_csv(db_stats_file, nrows=0, sep='\t', header=None)
    valid_col_count = len(header.columns)

    usecols = [0, 2, 7, 8]
    names = ['DB_level', 'Taxid', 'TotalLength', 'Note']

    if valid_col_count == 10:
        usecols=[0, 2, 7, 8, 9]
        names=['DB_level', 'Taxid', 'TotalLength', 'GenomeSize', 'Note']

    # Set header to None to support files without headers
    df_stats = pd.read_csv(db_stats_file,
                           low_memory=False,
                           sep='\t',
                           header=None,
                           usecols=usecols,
                           names=names,
                           dtype={'DB_level': str, 'Taxid': str},
                           index_col='Taxid')

    # If 'Note' column is not present, create it with empty strings
    if not 'Note' in df_stats:
        df_stats['Note'] = ''
    if not 'GenomeSize' in df_stats:
        df_stats['GenomeSize'] = 0

    # Remove the row with index 'Taxid' if it exists
    # This is to handle the case when the stats file having the header
    if 'Taxid' in df_stats.index:
        df_stats = df_stats.drop('Taxid')

    # Make sure the format is consistent
    try:
        df_stats['TotalLength'] = pd.to_numeric(df_stats['TotalLength'], errors='raise')
    except ValueError:
        logging.error(f"Error processing stats file. Please check the format of the file.")
        sys.exit(1)

    # This is to handle the case when the stats file does not have GenomeSize column
    # In that case, the 'Note' column will be loaded as 'GenomeSize' and filled with 0
    df_stats['GenomeSize'] = pd.to_numeric(df_stats['GenomeSize'], errors='coerce')
    df_stats['GenomeSize'] = df_stats['GenomeSize'].fillna(0).astype(int)

    return df_stats


def print_message(msg: str, silent: bool, start: float, logfile: Path, errorout: int = 0):
    """
    Print and log a timestamped message.

    Writes a message to the log file and optionally to stderr. Can also
    terminate the program with an error message.

    Parameters:
        msg (str): Message to print
        silent (bool): If True, suppress output to stderr
        start (float): Start time for timestamp calculation
        logfile (Path): Path to the log file
        errorout (int): If non-zero, exit with error after printing

    Returns:
        None

    Raises:
        SystemExit: If errorout is non-zero
    """
    message = "[%s] %s\n" % (time_spend(start), msg)

    with logfile.open("a", encoding="utf-8") as f:
        f.write(message)

    if errorout:
        sys.exit(message)
    elif not silent:
        sys.stderr.write(message)


def main(argvs):
    """
    Run the profiling/extraction workflow with a validated argparse.Namespace.
    """
    from gottcha.gottcha2 import __version__

    global logfile
    global begin_t
    global df_stats
    global acc_list

    direct_ont_flag = bool(argvs.nanopore and not argvs.ont_chunk)
    begin_t  = time.time()
    bamfile  = Path(argvs.bam) if argvs.bam else Path(argvs.outdir) / f"{argvs.prefix}.gottcha_{argvs.dbLevel}.bam"
    samfile  = Path(argvs.outdir) / f"{argvs.prefix}.gottcha_{argvs.dbLevel}.sam"
    logfile  = Path(argvs.outdir) / f"{argvs.prefix}.gottcha_{argvs.dbLevel}.log"
    set_start_method("fork") # for default multiprocessing method
    acc_list = set()
    split_read_flag = False
    multi_part_index_flag = False
    res_df = pd.DataFrame() # aggregated restuls
    logfile_prev = ""
    reciprocal_groups = {}

    logging_level = logging.WARNING

    if argvs.debug:
        logging_level = logging.DEBUG
    elif argvs.silent:
        logging_level = logging.FATAL
    elif argvs.verbose:
        logging_level = logging.INFO

    logging.basicConfig(
        level=logging_level,
        format='[%(asctime)s] [%(levelname)s] [%(module)s] %(message)s',
        datefmt='%Y%m%d %H:%M:%S',
   )

    #dependency check
    if sys.version_info < (3,9):
        sys.exit("[ERROR] Python 3.9 or above is required.")

    dependency_check("minimap2")
    dependency_check("samtools")
    if argvs.fast:
        dependency_check("sylph")

    #prepare output object
    argvs.relAbu = argvs.relAbu.upper()
    outfile_full = Path(argvs.outdir) / f"{argvs.prefix}.full.tsv"
    outfile_lineage = Path(argvs.outdir) / f"{argvs.prefix}.lineage.tsv"
    outfile_mpa = Path(argvs.outdir) / f"{argvs.prefix}.mpa.tsv"

    #create output directory if not exists
    Path(argvs.outdir).mkdir(parents=True, exist_ok=True)

    outfile = Path(argvs.outdir) / f"{argvs.prefix}.tsv"
    if argvs.format == "csv":
        outfile = Path(argvs.outdir) / f"{argvs.prefix}.csv"
    elif argvs.format == "biom":
        outfile = Path(argvs.outdir) / f"{argvs.prefix}.biom"

    out_fp = outfile.open("w", encoding="utf-8")

    if argvs.extractOnly:
        # repalce bamfile name from ".gottcha_\w+.bam" to ".log"
        logfile_prev = bamfile.with_suffix(".log")
        (mi, mf, mg, sni_argv) = (None, None, None, None)

        # if match criteria (mi/mf/mg) are not provided, load them from the log file
        if logfile_prev.is_file():
            (mi, mf, mg, sni_argv) = extract_reads.load_criteria_from_log(logfile_prev)
        else:
            logfile_prev = None

        if (argvs.sniScore is None) and sni_argv is not None:
            argvs.sniScore = sni_argv

        if (argvs.matchIdentity is None) and mi and mi >= 0:
            argvs.matchIdentity = mi
        else:
            argvs.matchIdentity = 0

        if (argvs.matchFraction is None) and mf and mf >= 0:
            argvs.matchFraction = mf
        else:
            argvs.matchFraction = 0

        if (argvs.matchLength is None) and mg and mg >= 0:
            argvs.matchLength = mg
        else:
            argvs.matchLength = 0

    (sni_score_cutoff, sni_score_species, sni_score_strain) = [float(x) for x in argvs.sniScore.split(',')]

    # display the command line
    logging.info(' '.join(sys.argv))

    print_message(f"GOTTCHA (v{__version__})", argvs.silent, begin_t, logfile)
    print_message(f"Arguments and dependencies checked:", argvs.silent, begin_t, logfile)
    print_message(f" - Database           : {argvs.database}",    argvs.silent, begin_t, logfile)
    print_message(f" - Database level     : {argvs.dbLevel}",     argvs.silent, begin_t, logfile)
    print_message(f" - Abundance          : {argvs.relAbu}",      argvs.silent, begin_t, logfile)
    print_message(f" - Output directory   : {argvs.outdir}",      argvs.silent, begin_t, logfile)
    print_message(f" - Output prefix      : {argvs.prefix}",      argvs.silent, begin_t, logfile)
    print_message(f" - Threads            : {argvs.threads}",     argvs.silent, begin_t, logfile)
    print_message(f" - Fast mode          : {argvs.fast}",        argvs.silent, begin_t, logfile)
    if argvs.input:
        print_message(f" - Input Reads        : {argvs.input}",     argvs.silent, begin_t, logfile)
    if argvs.bam:
        print_message(f" - Input BAM File     : {bamfile}",           argvs.silent, begin_t, logfile)
    if argvs.nanopore:
        print_message(f" - Nanopore Mode      : Enabled",              argvs.silent, begin_t, logfile)
    if argvs.errorRate:
        print_message(f" - Read Error Rate    : {argvs.errorRate}", argvs.silent, begin_t, logfile)
    if argvs.sigList:
        print_message(f" - Sig-of-Int List    : {argvs.sigList}", argvs.silent, begin_t, logfile)
    if argvs.sigList:
        print_message(f" - Sig-of-Int Action  : {argvs.sigListAction}", argvs.silent, begin_t, logfile)
    if argvs.minCov > 0:
        print_message(f" - Minimal SIG Cov    : {argvs.minCov}",      argvs.silent, begin_t, logfile)
    if argvs.minLen > 0:
        print_message(f" - Minimal SIG Length : {argvs.minLen}",      argvs.silent, begin_t, logfile)
    if argvs.minReads > 0:
        print_message(f" - Minimal Reads      : {argvs.minReads}",    argvs.silent, begin_t, logfile)
    if argvs.extract:
        print_message(f" - Extract Taxa       : {argvs.extract}",     argvs.silent, begin_t, logfile)
    if argvs.extractOnly:
        print_message(f" - Extract Only       : {argvs.extractOnly}", argvs.silent, begin_t, logfile)
    if argvs.maxZscore > 0:
        print_message(f" - Maximal zScore     : {argvs.maxZscore}",   argvs.silent, begin_t, logfile)
    if logfile_prev:
        print_message(f" - Load criteria from : {logfile_prev}",      argvs.silent, begin_t, logfile)
    if argvs.matchIdentity != None:
        print_message(f" - Min Match Identity : {argvs.matchIdentity}", argvs.silent, begin_t, logfile)
    if argvs.matchFraction != None:
        print_message(f" - Min Match Fraction : {argvs.matchFraction}", argvs.silent, begin_t, logfile)
    if argvs.matchLength != None:
        print_message(f" - Min Match Length   : {argvs.matchLength}", argvs.silent, begin_t, logfile)
    if argvs.sniScore != None:
        print_message(f" - SNI-score (g,s,n)  : {argvs.sniScore}",    argvs.silent, begin_t, logfile)

    #load taxonomy for taxonomic aggregation and annotation
    if not argvs.extractOnly:
        print_message("Loading taxonomy information...", argvs.silent, begin_t, logfile)

        if Path(argvs.database + ".tax.tsv").exists():
            custom_taxa_tsv = Path(argvs.database + ".tax.tsv")
        elif Path(argvs.database + ".taxa").exists():
            custom_taxa_tsv = Path(argvs.database + ".taxa")

        taxonomy.loadTaxonomy(cus_taxonomy_file=custom_taxa_tsv, auto_download=False)
        print_message(f" - {len(taxonomy.taxNames):,} taxa loaded.", argvs.silent, begin_t, logfile)

        #load database stats
        print_message("Loading database stats...", argvs.silent, begin_t, logfile)
        if Path(argvs.database + ".stats").is_file():
            df_stats = load_database_stats(argvs.database+".stats")
        else:
            print_message(f"ERROR: {argvs.database+'.stats'} not found.", argvs.silent, begin_t, logfile, errorout=1)

        print_message(f" - {df_stats.shape[0]:,} entries loaded.", argvs.silent, begin_t, logfile)
        print_message(f" - signatures at {df_stats['DB_level'].unique().tolist()} levels loaded.", argvs.silent, begin_t, logfile)

    if argvs.sigList:
        print_message("Loading accession#s of interest list...", argvs.silent, begin_t, logfile)
        acc_list = load_acc_list(argvs.sigList)
        print_message(f" - {len(acc_list):,} accession/signature of interest loaded.", argvs.silent, begin_t, logfile)

    # Summary of the Main Process:
    #
    # Input Reads
    #     ↓
    # [Nanopore Preprocessing] (optional)
    #     ↓
    # [Run fast query] (sylph; optionall; if fast mode is on)
    #     ↓
    # [Extract queried signatures] (optionall; if fast mode is on)
    #     ↓
    # Read Mapping (minimap2)
    #     ↓
    # Alignments (SAM File)
    #     ↓
    # [Remove Multiple Hits] (if multi-part index)
    #     ↓
    # [Remove Inconsistent Chunks] (if nanopore)
    #     ↓
    # BAM Conversion + Indexing
    #     ↓
    # Parse & Filter Alignments
    #     ↓
    # Group to Strains
    #     ↓
    # Aggregate Taxonomy
    #     ↓
    # Generate Reports
    #     ↓
    # [Extract Reads] (optional)

    if argvs.input:
        # if fast mode is on, run Sylph query to prefilter the reference genomes and create a smaller reference for read mapping; 
        # otherwise, use the full database index for read mapping
        minimap2_index = "" if argvs.fast else f"{argvs.database}.mmi"
        
        # The original input reads for Sylph sketch and query; Not the split reads for minimap2 mapping if nanopore option is on
        sylph_input = argvs.input

        # if nanopore option is on, preprocessing reads
        if argvs.nanopore:
            print_message("Checking nanopore read files...", argvs.silent, begin_t, logfile)
            if direct_ont_flag:
                print_message(" - Direct mapping ONT reads", argvs.silent, begin_t, logfile)
            else:
                print_message(" - Splitting ONT reads to chunks", argvs.silent, begin_t, logfile)
                argvs.input = ont_utils.preprocess_nanopore_reads(argvs.input, argvs.outdir, argvs.prefix, argvs.silent)
                split_read_flag = True

        if argvs.fast:
            print_message("Prefiltering reference genomes...", argvs.silent, begin_t, logfile)
            sylph_db = f"{argvs.database}.syldb"
            g2_archive = f"{argvs.database}.zip"
            sylph_query_tsv = Path(argvs.outdir) / f"{argvs.prefix}.sylph_query.tsv"
            queried_signatures_file = Path(argvs.outdir) / f"{argvs.prefix}.sylph_queried_signatures.txt"
            extracted_reference = Path(argvs.outdir) / f"{argvs.prefix}.sylph_extracted.fa.gz"
            
            # extract subsample (cXXX) rate from sylph_db string, default set to 100
            subsampling_rate = 100
            match = re.search(r'c(\d+)\.', sylph_db)
            if match:
                subsampling_rate = int(match.group(1))

            fast_min_kmer = argvs.fast_min_kmer
            fast_min_ani = argvs.fast_min_ani

            # Run Sylph sketch if the input file is in FASTA format
            if Path(sylph_input[0]).name.endswith(('.fa', '.fasta', '.fa.gz', '.fna', '.fna.gz', '.fasta.gz')):
                print_message("Generating sketchs for FASTA input reads...", argvs.silent, begin_t, logfile)
                try:
                    sylph_result = prefilter.run_sylph_sketch(
                        read_file=sylph_input[0],
                        outdir=str(argvs.outdir),
                        threads=argvs.threads,
                        subsampling_rate=subsampling_rate,
                    )
                except (FileNotFoundError, subprocess.CalledProcessError) as e:
                    print_message(f"ERROR: prefiltering failed: {e}", argvs.silent, begin_t, logfile, errorout=1)

                with logfile.open("a", encoding="utf-8") as f:
                    if sylph_result.stdout:
                        f.write(sylph_result.stdout)
                    if sylph_result.stderr:
                        f.write(sylph_result.stderr)
                
                sylph_input = [str(Path(argvs.outdir) / f"{Path(sylph_input[0]).name}.sylsp")]

            # Run Sylph query to get the list of signatures that are likely present in the input reads
            try:
                sylph_result = prefilter.run_sylph_query(
                    database=sylph_db,
                    reads=sylph_input,
                    output=str(sylph_query_tsv),
                    threads=argvs.threads,
                    subsampling_rate=subsampling_rate,
                    minimum_kmer=fast_min_kmer,
                    minimum_ani=fast_min_ani,
                    read_seq_id=float(100-argvs.errorRate*100)
                )
            except (FileNotFoundError, subprocess.CalledProcessError) as e:
                print_message(f"ERROR: prefiltering failed: {e}", argvs.silent, begin_t, logfile, errorout=1)

            with logfile.open("a", encoding="utf-8") as f:
                if sylph_result.stdout:
                    f.write(sylph_result.stdout)
                if sylph_result.stderr:
                    f.write(sylph_result.stderr)

            try:
                pd.read_csv(sylph_query_tsv, 
                            sep='\t', 
                            usecols=['Genome_file'], 
                            dtype={'Genome_file': str})['Genome_file'].dropna().str.strip().to_csv(queried_signatures_file, index=False, header=False)
            except pd.errors.EmptyDataError:
                queried_signatures = pd.Series(dtype=str)
            except (FileNotFoundError, ValueError) as e:
                print_message(f"ERROR: unable to parse Sylph query output {sylph_query_tsv}: {e}", argvs.silent, begin_t, logfile, errorout=1)

            filenames = sig_archive.read_file_list(queried_signatures_file, filename_only=True)
            print_message(f" - Identified {len(filenames):,} reference genomes.", argvs.silent, begin_t, logfile)
            
            if len(filenames) == 0:
                print_message("No references identified. GOTTCHA2 stopped.", argvs.silent, begin_t, logfile)
                sys.exit(0)

            # Extract those signatures from the archive to create a smaller reference for read mapping
            extracted_content, processed_files, skipped_files = sig_archive.quick_concat(g2_archive,
                                                                                         separator=str('\n').encode('utf-8'),
                                                                                         skip_missing=False, 
                                                                                         filenames=filenames)

            extracted_reference.write_bytes(extracted_content)
            if extracted_reference.stat().st_size == 0:
                print_message(
                    f"ERROR: no queried signatures could be extracted from {sylph_db}; {extracted_reference} is empty.",
                    argvs.silent,
                    begin_t,
                    logfile,
                    errorout=1
                )

            if skipped_files:
                print_message(f" - {len(skipped_files):,} queried signatures were not found in the archive.", argvs.silent, begin_t, logfile)
            print_message(f" - {len(processed_files):,} reference genomes extracted.", argvs.silent, begin_t, logfile)
            minimap2_index = str(extracted_reference)

        print_message("Running read-mapping...", argvs.silent, begin_t, logfile)
        exitcode, cmd, input_read_count, multi_part_index_flag = read_mapping.minimap2(
            argvs.input,
            minimap2_index,
            argvs.threads,
            argvs.m2_options,
            argvs.presetx,
            samfile,
            logfile,
            allow_secondary=(argvs.secondary == 'yes'),
            max_secondary=argvs.max_secondary,
            secondary_ratio=argvs.secondary_ratio,
        )
        logging.info(f"COMMAND: {cmd}")

        if exitcode != 0:
            # if size of the samfile is zero
            print_message(f"Logfile saved to {logfile}.", argvs.silent, begin_t, logfile)
            sys.exit("[%s] ERROR: error occurred while running read mapping (exit: %s).\n" % (time_spend(begin_t), exitcode))
        else:
            print_message(f" - {input_read_count:,} input reads processed.", argvs.silent, begin_t, logfile)
            print_message(f"Mapped SAM file saved to {samfile}.", argvs.silent, begin_t, logfile)
        gc.collect()

    # remove multiple hits
    if multi_part_index_flag:
        # remove multiple hits from the SAM file
        print_message("Removing multiple hits from SAM file...", argvs.silent, begin_t, logfile)
        samfile_temp = Path(argvs.outdir) / f"{argvs.prefix}.gottcha_{argvs.dbLevel}.sam.temp"
        flag, aln_count, top_hits_count = read_mapping.post_processing_sam(samfile, samfile_temp)
        if flag:
            samfile_temp.rename(samfile)
            # Note:
            # When input of the gottcha2 is a SAM file and new outdir/prefix is provided, the output will be saved to that location.
            # If not, the output will overwrite the original SAM file.
            print_message(f" - {aln_count:,} total alignments", argvs.silent, begin_t, logfile)
            print_message(f" - {top_hits_count:,} best hits among index-partitions", argvs.silent, begin_t, logfile)

        gc.collect()

    if Path(samfile).is_file() and argvs.reciprocal_groups == 'yes':
        print_message("Resolving reciprocal relationships from SAM file...", argvs.silent, begin_t, logfile)
        reciprocal_groups = reciprocal_graph.reciprocal_relationships_from_sam(samfile, min_alen=argvs.matchLength)

        for group in reciprocal_groups:
            logging.debug(f"{group}:")
            for species_taxid in reciprocal_groups[group]:
                logging.debug(f" - {taxonomy.taxid2name(species_taxid)} ({species_taxid}); g__{taxonomy.taxid2nameOnRank(species_taxid, 'genus')}")

        tol_reciprocal_groups = len(reciprocal_groups)
        print_message(f" - {tol_reciprocal_groups:,} reciprocal groups identified", argvs.silent, begin_t, logfile)
        gc.collect()

    # processing alignments and generate results
    if not argvs.extractOnly:
        if os.path.isfile(os.path.abspath(samfile)):
            print_message("Converting to BAM file...", argvs.silent, begin_t, logfile)
            sam_to_bam.convert_sam_to_bam(input_sam=os.path.abspath(samfile),
                                          output_bam=os.path.abspath(bamfile),
                                          threads=argvs.threads,
                                          quiet=argvs.silent)
            print_message(f"BAM file saved to {bamfile}...", argvs.silent, begin_t, logfile)
            
            file_path = Path(samfile)
            if file_path.exists():
                file_path.unlink()
            gc.collect()

        if Path(bamfile).exists() and Path(f"{bamfile}.bai").exists():
            print_message("Processing alignments...", argvs.silent, begin_t, logfile)
            ref_chunk_results = process_bam.parse_aln_from_bam(
                bam_path=bamfile,
                processes=argvs.threads,
                min_frac=argvs.matchFraction,
                min_idt=argvs.matchIdentity,
                min_alen=argvs.matchLength,
                include_secondary=(argvs.secondary == 'yes'),
                include_supplementary=(argvs.secondary == 'yes'),
                split_read_flag=split_read_flag
            )

            str_df, soi_read_count = aggregate_results.group_refs_to_strains(ref_chunk_results, acc_list, argvs.sigListAction, df_stats)

            tol_read_count = str_df['READ_COUNT'].sum()
            tol_invalid_match_count = str_df['INVALID_ALNS'].sum()

            print_message(f" - {tol_invalid_match_count:,} alignments did not meet matching criteria", argvs.silent, begin_t, logfile)
            print_message(f" - {tol_read_count:,} qualified reads processed", argvs.silent, begin_t, logfile)

            if not tol_read_count:
                print_message("No qualified alignments found. Stopping.", argvs.silent, begin_t, logfile)
                sys.exit(0)

            # aggregate the results
            _args = (str_df,
                     argvs.relAbu,
                     argvs.dbLevel,
                     argvs.minCov,
                     argvs.minReads,
                     argvs.minLen,
                     argvs.maxZscore,
                     sni_score_species,
                     sni_score_strain,
                     sni_score_cutoff,
                     argvs.errorRate,
                     df_stats,
                     reciprocal_groups)
            res_df, soi_read_count = aggregate_results.aggregate_taxonomy(*_args)

            if acc_list:
                print_message(f" - {soi_read_count:,} reads mapped to accession-of-interest", argvs.silent, begin_t, logfile)
                read_count_after_soi = tol_read_count
                if argvs.sigListAction == 'filter_out':
                    read_count_after_soi = tol_read_count - soi_read_count
                elif argvs.sigListAction == 'filter_in':
                    read_count_after_soi = soi_read_count
                print_message(f" - {read_count_after_soi:,} reads after applying accession-of-interest action ({argvs.sigListAction})", argvs.silent, begin_t, logfile)

            print_message("Done taxonomy aggregation.", argvs.silent, begin_t, logfile)

            if not len(res_df):
                print_message("No qualified taxonomy profiled.", argvs.silent, begin_t, logfile)
            else:
                # generate output results
                if argvs.format == "biom":
                    report.generate_biom_file(res_df, out_fp, argvs.dbLevel, argvs.prefix)
                else:
                    report.generate_report_file(res_df, out_fp, outfile_full, argvs.format)
                # generate lineage file
                target_idx = (res_df['LEVEL']==argvs.dbLevel) & \
                                (res_df['NOTE'].str.contains('Filtered out', na=False) == False) & \
                                (res_df['NOTE'].str.contains('Not shown', na=False) == False)
                target_df = res_df.loc[target_idx, ['ABUNDANCE','TAXID']]
                tax_num = len(target_df)

                print_message(f"{tax_num} qualified {argvs.dbLevel} profiled.", argvs.silent, begin_t, logfile)

                if tax_num:
                    report.generate_lineage_file(target_df, outfile_lineage)

                    if argvs.mpa:
                        target_df = res_df.loc[target_idx, ['TAXID', 'REL_ABUNDANCE', 'REL_ABUNDANCE_GC','READ_COUNT', 'SIG_COV']]
                        report.generate_mpa_file(target_df, outfile_mpa)
                        print_message(f"MPA format file saved to {outfile_mpa}.", argvs.silent, begin_t, logfile)

                print_message(f"Results saved to {outfile}.", argvs.silent, begin_t, logfile)
        else:
            print_message(f"ERROR: BAM file {bamfile} or its index not found.", argvs.silent, begin_t, logfile, errorout=1)
            print_message("GOTTCHA2 stopped.", argvs.silent, begin_t, logfile)
            sys.exit(0)

    # extracting reads
    if argvs.extract:
        (taxa_arg, max_per_taxon, out_format) = (argvs.extract.split(':', maxsplit=2) + ['all', 'fasta'])[:3]

        print_message(f"Extracting {max_per_taxon} sequences per taxa in {out_format} format...", argvs.silent, begin_t, logfile)

        full_report_file = ""

        if max_per_taxon.isdigit() or max_per_taxon == 'all':
            max_per_taxon = int(max_per_taxon) if max_per_taxon != 'all' else 0

        if argvs.extractOnly:
            # repalce bamfile name to replace from ".gottcha_\w+.bam" to ".full.tsv"
            full_report_file = re.sub(r"\.gottcha_\w+\.bam$", ".full.tsv", str(bamfile))

        taxa_dict, ref_to_extract_taxid = extract_reads.parse_taxids(taxa_arg, 
                                                                     res_df, 
                                                                     full_report_file, 
                                                                     sni_score_cutoff,
                                                                     sni_score_species, 
                                                                     sni_score_strain
                                                                     )

        if not len(ref_to_extract_taxid):
            print_message("No qualified taxonomy profiled.", argvs.silent, begin_t, logfile)

        outfile = Path(argvs.outdir) / f"{argvs.prefix}.extract.{out_format.lower()}"

        _args = (os.path.abspath(bamfile),
                    taxa_dict,
                    ref_to_extract_taxid,
                    outfile.open("w", encoding="utf-8"),
                    argvs.threads,
                    argvs.matchFraction,
                    argvs.matchIdentity,
                    argvs.matchLength,
                    max_per_taxon,
                    acc_list,
                    argvs.sigListAction,
                    out_format)
        taxon_count, seq_count = extract_reads.extract_sequences_by_taxonomy(*_args)
        print_message(f"Done extracting {seq_count} sequences from {taxon_count} taxa to '{outfile}'.",
                        argvs.silent, begin_t, logfile)

if __name__ == '__main__':
    from gottcha.gottcha2 import cli
    args = sys.argv[1:]
    if not args or args[0].startswith('-'):
        args = ['profile', *args]
    cli(args)
