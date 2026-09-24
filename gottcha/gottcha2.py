#!/usr/bin/env python3
__version__   = "2.5.1"
__author__    = "Po-E (Paul) Li, B-GEN, Bioscience Division, Los Alamos National Laboratory"

import argparse as ap
import os
import sys
from pathlib import Path
from re import search

if __package__ in (None, ''):
    # Use this checkout when running the CLI as a script.
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from gottcha.utils import profile, download, sam_to_bam, coverage_browser

def parse_args(ver, args):
    """
    Parse and validate command line arguments for GOTTCHA2.

    This function sets up the argument parser, defines all possible command-line
    options, parses the provided arguments, and performs validation to ensure
    the configuration is valid and complete.

    Parameters:
        ver (str): Version string to display in help messages
        args (list): Command line arguments to parse

    Returns:
        argparse.Namespace: Object containing all validated arguments

    Raises:
        SystemExit: If validation fails or --version is specified
    """
    command = args[0] if args and not args[0].startswith('-') else None
    parser_args = args[1:] if command else args
    taxonomic_levels = [
        'superkingdom', 'phylum', 'class', 'order',
        'family', 'genus', 'species', 'strain',
    ]

    p = ap.ArgumentParser(
        prog=f'gottcha2 {command}' if command else 'gottcha2 profile',
        formatter_class=ap.RawTextHelpFormatter,
        description=(
            "Genomic Origin Through Taxonomic CHAllenge (GOTTCHA) is an annotation-independent,\n"
            "signature-based metagenomic taxonomic profiling tool with substantially low false\n"
            "discovery rates. This command maps input reads to precomputed signature databases\n"
            f"using minimap2 and profiles the organisms present in a sample. (Version: {ver})"
        ),
    )

    input_group = p.add_argument_group('Input source')
    eg = input_group.add_mutually_exclusive_group(required=True)
    eg.add_argument(
        '-i', '--input',
        metavar='FASTQ',
        nargs='+',
        type=str,
        help='Input FASTQ/FASTA file(s). Separate multiple files with spaces.',
    )
    eg.add_argument(
        '-b', '--bam',
        metavar='BAMFILE',
        type=str,
        help='Input sorted and indexed BAM file.',
    )

    database_group = p.add_argument_group('Database')
    database_group.add_argument(
        '-d', '--database',
        metavar='GOTTCHA2_DB',
        type=str,
        default=None,
        help='Path to the GOTTCHA2 database prefix, index file, or database directory.',
    )
    database_group.add_argument(
        '-l', '--dbLevel',
        metavar='LEVEL',
        type=str,
        default='',
        choices=taxonomic_levels,
        help=(
            'Taxonomic level of the input database.\n'
            'Choices: superkingdom, phylum, class, order, family, genus, species, strain.\n'
            'Auto-detected when the database prefix contains a rank, e.g. GOTTCHA_db.species.'
        ),
    )

    platform_group = p.add_argument_group('Sequencing and mapping')
    platform_group.add_argument(
        '-np', '--nanopore',
        action='store_true',
        help=(
            'Treat input reads as Oxford Nanopore (ONT) reads.\n'
            'Enables read preprocessing and ONT defaults: -er 0.03 -mi 0.85 -mf 0.85 -mg 100.'
        ),
    )
    platform_group.add_argument(
        '--ont-chunk',
        action='store_true',
        help=(
            'Split ONT reads to 150-bp chunks and map them to GOTTCHA2 signature fragments.\n'
            'If not specified, the default direct mapping workflow is used.'
        ),
    )
    platform_group.add_argument(
        '-xm', '--presetx',
        metavar='STR',
        type=str,
        default=None,
        choices=['sr', 'map-pb', 'map-ont', 'lr:hq'],
        help='minimap2 preset passed with -x. [default: lr:hq for --nanopore; sr otherwise]',
    )
    platform_group.add_argument(
        '--secondary', choices=['yes', 'no'], default='no',
        help='Allow minimap2 secondary candidates. [default: no]',
    )
    platform_group.add_argument(
        '--max-secondary', type=int, default=10,
        help='Maximum minimap2 secondary candidates per primary alignment. [default: 10]',
    )
    platform_group.add_argument(
        '--secondary-ratio', type=float, default=0.9,
        help='Minimum minimap2 secondary/primary chaining-score ratio. [default: 0.9]',
    )
    platform_group.add_argument(
        '--m2-options',
        metavar='STR',
        type=str,
        default='auto',
        help='Additional minimap2 options. Use with care. [default: platform-specific auto settings]',
    )
    platform_group.add_argument(
        '-er', '--errorRate',
        metavar='FLOAT',
        type=float,
        help='Estimated sequencing error rate. [default: 0.005, or 0.03 with --nanopore]',
    )
    platform_group.add_argument(
        '-t', '--threads',
        metavar='INT',
        type=int,
        default=1,
        help='Number of threads. [default: 1]',
    )

    alignment_group = p.add_argument_group('Alignment thresholds')
    alignment_group.add_argument(
        '-mi', '--matchIdentity',
        metavar='FLOAT',
        type=float,
        help='Minimum identity (0.0-1.0) required for a valid match. [default: 0.95, or 0.85 with --nanopore]',
    )
    alignment_group.add_argument(
        '-mf', '--matchFraction',
        metavar='FLOAT',
        type=float,
        help='Minimum aligned fraction (0.0-1.0) of the read or signature fragment. [default: 0.95, or 0.85 with --nanopore]',
    )
    alignment_group.add_argument(
        '-mg', '--matchLength',
        metavar='INT',
        type=int,
        help='Minimum alignment length in bp required for a valid match. [default: 100]',
    )

    profiling_group = p.add_argument_group('Profiling thresholds')
    profiling_group.add_argument(
        '-ss', '--sniScore',
        metavar='FLOAT[,FLOAT,FLOAT]',
        type=str,
        help=(
            'Signature nucleotide identity (SNI) thresholds for taxonomic aggregation.\n'
            'One value applies to all ranks; two values append strain default 0.99;\n'
            'three values mean other ranks, species, and strain. [default: 0.8,0.95,0.99]'
        ),
    )
    profiling_group.add_argument(
        '-Mc', '-mc', '--minCov',
        metavar='FLOAT',
        type=float,
        default=0,
        help='Minimum signature coverage used in abundance calculation. [default: 0]',
    )
    profiling_group.add_argument(
        '-Mr', '-mr', '--minReads',
        metavar='INT',
        type=int,
        default=0,
        help='Minimum read count used in abundance calculation. [default: 0]',
    )
    profiling_group.add_argument(
        '-Ml', '-ml', '--minLen',
        metavar='INT',
        type=int,
        default=0,
        help='Minimum signature length used in abundance calculation. [default: 0]',
    )
    profiling_group.add_argument(
        '-Mz', '-mz', '--maxZscore',
        metavar='FLOAT',
        type=float,
        default=0,
        help='Maximum estimated z-score for mapped-region depths; 0 disables the filter. [default: 0]',
    )
    profiling_group.add_argument(
        '-nc', '--noCutoff',
        action='store_true',
        help='Disable profiling-stage cutoffs. Equivalent to -Mc 0 -Mr 0 -Ml 0 -Mz 0 -ss 0,0,0.',
    )
    profiling_group.add_argument(
        '-r', '--relAbu',
        metavar='FIELD',
        type=str,
        default='DEPTH',
        choices=['DEPTH', 'READ_COUNT', 'GENOMIC_CONTENT_EST'],
        help='Field used to calculate relative abundance. [default: DEPTH]',
    )
    profiling_group.add_argument(
        '--reciprocal-groups',
        choices=['yes', 'no'], 
        default='no',
        help='(EXPERIMENTAL) Enable or disable reciprocal groups. [default: no]',
    )

    signature_group = p.add_argument_group('Signature-of-interest filtering')
    signature_group.add_argument(
        '-sl', '--sigList',
        metavar='FILE',
        type=str,
        help='File containing accessions/signatures of interest, one per line.',
    )
    signature_group.add_argument(
        '-sa', '--sigListAction',
        choices=['filter_out', 'filter_in', 'report_only'],
        default='report_only',
        type=str,
        help=(
            'Action for aligned reads mapped to signatures of interest:\n'
            '  filter_out  discard matching reads\n'
            '  filter_in   keep only matching reads\n'
            '  report_only report matching reads as SOI_READ_COUNT without filtering [default]'
        ),
    )

    extraction_group = p.add_argument_group('Read extraction')
    extraction_group.add_argument(
        '-e', '--extract',
        metavar='TAXON[,TAXON2,...]',
        type=str,
        default=None,
        help=(
            'Extract mapped reads for specific taxa to FASTA or FASTQ.\n'
            'Accepted forms:\n'
            "  1234,5678              comma-separated taxon IDs\n"
            "  @taxids.txt            one taxon ID per line\n"
            "  @taxids.txt:1000:fasta limit reads per taxon and choose fasta/fastq\n"
            "  all                    extract all matching taxa/reads"
        ),
    )
    extraction_group.add_argument(
        '-ef', '--extractFullRef',
        action='store_true',
        help="Extract up to 20 sequences per reference as FASTA. Equivalent to -e 'all:20:fasta'.",
    )
    extraction_group.add_argument(
        '-eo', '--extractOnly',
        action='store_true',
        help='Only extract reads from the alignment file; skip profiling.',
    )

    output_group = p.add_argument_group('Output')
    output_group.add_argument(
        '-fm', '--format',
        metavar='STR',
        type=str,
        default='tsv',
        choices=['tsv', 'csv', 'biom'],
        help='Output format. [default: tsv]',
    )
    output_group.add_argument(
        '-o', '--outdir',
        metavar='DIR',
        type=str,
        default='.',
        help='Output directory. [default: .]',
    )
    output_group.add_argument(
        '-p', '--prefix',
        metavar='STR',
        type=str,
        help='Output file prefix. [default: input file prefix]',
    )
    output_group.add_argument(
        '--mpa',
        action='store_true',
        help='Generate output in MetaPhlAn format (*.mpa).',
    )

    fast_group = p.add_argument_group('Fast profile')
    fast_group.add_argument(
        '--fast',
        action='store_true',
        help='Enable fast-profile mode.',
    )
    fast_group.add_argument(
        '--fast-min-kmer',
        metavar='INT',
        type=int,
        default=5,
        help='Minimum k-mer size for fast-profile prefiltering. [default: 5]',
    )
    fast_group.add_argument(
        '--fast-min-ani',
        metavar='INT',
        type=int,
        default=80,
        help='Minimum ANI for fast-profile prefiltering. [default: 80]',
    )
    
    logging_group = p.add_argument_group('Logging')
    logging_group.add_argument(
        '--silent',
        action='store_true',
        help='Disable status messages.',
    )

    logging_group.add_argument(
        '--verbose',
        action='store_true',
        help='Show verbose status messages.',
    )
    logging_group.add_argument(
        '--debug',
        action='store_true',
        help='Enable debug logging and keep temporary files.',
    )
    
    eg.add_argument(
        '-v', '--version',
        action='store_true',
        help='Print version number and exit.',
    )
    p.set_defaults(fast=(command == 'fast-profile'), extractOnly=(command == 'extract'))
    args_parsed = p.parse_args(parser_args)

    """
    Checking options
    """
    if args_parsed.version:
        print(ver)
        sys.exit(0)

    if args_parsed.extract and args_parsed.extractFullRef:
        p.error('--extract and --extractFullRef are incompatible options.')

    if args_parsed.input and args_parsed.bam:
        p.error('--input / --bam are incompatible options.')

    if not args_parsed.extractOnly:
        if not args_parsed.database:
            p.error('--database option is missing.')
        if args_parsed.sniScore is None:
            args_parsed.sniScore='0.9,0.95,0.99'

    # Auto-detect database path and prefix, and check the existence of input files and database index
    if args_parsed.database:
        #assign default path for database name
        db_extfn = "syldb" if args_parsed.fast else "mmi"

        # find the database index file if a directory is provided, and set the database prefix accordingly;
        if Path(args_parsed.database).is_dir():
            dbs = list(Path(args_parsed.database).glob(f"*.{db_extfn}"))
            if len(dbs) > 1:
                p.error(f'Multiple .{db_extfn} files found in {args_parsed.database}. Please specify one with database prefix.')
            elif len(dbs) == 0:
                p.error(f'No .{db_extfn} file found in {args_parsed.database}. Please specify the database prefix or the path to the .{db_extfn} file.')
            else:
                args_parsed.database = str(dbs[0])

        # check if the database file is provided with or without the extension, and remove the extension if provided
        if args_parsed.database.endswith(f'.{db_extfn}'):
            args_parsed.database = args_parsed.database.replace(f'.{db_extfn}', '')

        # Only check the existence of the database index file if input reads are provided
        if args_parsed.input:
            if not Path(f'{args_parsed.database}.{db_extfn}').is_file():
                p.error(f'Database index {args_parsed.database}.{db_extfn} not found.')

        if args_parsed.fast:
            if not Path(f'{args_parsed.database}.{db_extfn}').is_file():
                p.error(f'Database index {args_parsed.database}.{db_extfn} not found.')
            if not Path(f'{args_parsed.database}.zip').is_file():
                p.error(f'Signature sequences file {args_parsed.database}.zip not found.')

        # check the existence of the taxonomic information file for the specified database
        if not Path(f'{args_parsed.database}.tax.tsv').is_file():
            p.error(f'Taxonomic file {args_parsed.database}.tax.tsv not found.')
        if not Path(f'{args_parsed.database}.stats').is_file():
            p.error(f'Database stats file {args_parsed.database}.stats not found.')

    if args_parsed.input:
        for path in args_parsed.input:
            if path == '-':
                p.error('--input does not support reading from stdin ("-"). Please provide a file path.')
            if not Path(path).is_file():
                p.error(f'Input file {path} not found.')

    if args_parsed.bam:
        if not Path(args_parsed.bam).is_file():
            p.error(f'BAM file {args_parsed.bam} not found.')

    if args_parsed.sigList:
        if not Path(args_parsed.sigList).is_file():
            p.error(f'Signature-of-interest list {args_parsed.sigList} not found.')
        args_parsed.sigList = os.path.abspath(args_parsed.sigList)

    if args_parsed.nanopore and args_parsed.input and len(args_parsed.input) != 1:
        p.error('--nanopore option requires a single input read file.')

    if not args_parsed.prefix:
        if args_parsed.input:
            name = search(r'([^\/\.]+)\..*$', args_parsed.input[0])
            args_parsed.prefix = name.group(1)
        elif args_parsed.bam:
            name = search(r'([^\/]+)\.\w+\.bam$', args_parsed.bam)
            args_parsed.prefix = name.group(1)
        else:
            args_parsed.prefix = "GOTTCHA_"

    if not args_parsed.dbLevel:
        if args_parsed.database:
            major_ranks = {"superkingdom":1,"phylum":2,"class":3,"order":4,"family":5,"genus":6,"species":7, "strain":8}
            parts = args_parsed.database.split('.')
            for part in parts:
                if part in major_ranks:
                    args_parsed.dbLevel = part
                    break
        elif args_parsed.bam:
            name = search(r'\.gottcha_(\w+).bam$', args_parsed.bam)
            try:
                args_parsed.dbLevel = name.group(1)
            except:
                pass

        if not args_parsed.dbLevel:
            p.error('--dbLevel is missing and cannot be auto-detected.')

    # Set SNI-SCORE default to 0.9, species 0.95, strain 0.99
    if args_parsed.sniScore:
        if args_parsed.sniScore.count(',') == 0:
            args_parsed.sniScore = ','.join([args_parsed.sniScore]*3)
        elif args_parsed.sniScore.count(',') == 1:
            args_parsed.sniScore = args_parsed.sniScore + ',0.99'
    else:
        args_parsed.sniScore = '0.9,0.95,0.99'

    # If mi/mf/mg are not specified, set them to default values based on whether --nanopore is specified
    # But if --extractOnly is specified, do not set default values for matchIdentity and matchFraction, the value will be load from the log file if not provided
    if args_parsed.matchIdentity is None:
        if not args_parsed.extractOnly:
            if args_parsed.nanopore:
                args_parsed.matchIdentity = 0.85
            else:
                args_parsed.matchIdentity = 0.95
    else:
        if args_parsed.matchIdentity < 0 or args_parsed.matchIdentity > 1:
            p.error('--matchIdentity must be between 0 and 1.')

    if args_parsed.matchFraction is None:
        if not args_parsed.extractOnly:
            if args_parsed.ont_chunk:
                args_parsed.matchFraction = 0.85
            elif args_parsed.nanopore:
                args_parsed.matchFraction = 0
            else:
                args_parsed.matchFraction = 0.95
    else:
        if args_parsed.matchFraction < 0 or args_parsed.matchFraction > 1:
            p.error('--matchFraction must be between 0 and 1.')

    if args_parsed.matchLength is None:
        if not args_parsed.extractOnly:
            args_parsed.matchLength = 100
    else:
        if args_parsed.matchLength < 0:
            p.error('--matchLength must be a non-negative integer.')

    if args_parsed.extractFullRef:
        args_parsed.extract = 'all:20:fasta'

    if args_parsed.extractOnly:
        error_message = ""
        if not args_parsed.extract:
            error_message += "--extract must be specified. "

        if not args_parsed.bam:
            error_message += "--bam must be specified. "

        if error_message:
            p.error(error_message)

    if args_parsed.ont_chunk and not args_parsed.nanopore:
        p.error('--ont-chunk requires --nanopore.')

    if args_parsed.max_secondary < 0:
        p.error('--max-secondary must be >= 0.')

    if not 0 <= args_parsed.secondary_ratio <= 1:
        p.error('--secondary-ratio must be between 0 and 1.')

    if args_parsed.presetx is None:
        args_parsed.presetx = 'lr:hq' if (args_parsed.nanopore and not args_parsed.ont_chunk) else 'sr'

    if args_parsed.m2_options == 'auto':
        if args_parsed.nanopore and not args_parsed.ont_chunk:
            # Fast mode builds the signature index on the fly and adds k24/w12
            # below, so two minimizers are a useful short-fragment safeguard.
            # Prebuilt GOTTCHA2 .mmi indexes retain k28/w24; allow one seed to
            # initiate DP for ~100-bp signatures.
            seed_chain = '-n2' if args_parsed.fast else '-n1'
            args_parsed.m2_options = f'{seed_chain} -m25 -s120 --no-long-join'
        else:
            args_parsed.m2_options = '-s120'

    if args_parsed.noCutoff:
        args_parsed.sniScore = '0,0,0'

    if not args_parsed.errorRate:
        if args_parsed.ont_chunk:
            args_parsed.errorRate = 0.03
        elif args_parsed.nanopore:
            args_parsed.errorRate = 0.01
        else:
            args_parsed.errorRate = 0.005

    return args_parsed


def usage():
    """Display usage information for GOTTCHA2 command-line tool."""
    print(f"""
GOTTCHA2 - Genomic Origin Through Taxonomic CHAllenge v{__version__}

Usage:
    gottcha2 <command> [options]

Commands:
    profile       Use GOTTCHA2 to profile metagenomic reads against a signature database

    fast-profile  Faster version of profile that uses a more aggressive prefiltering strategy to speed up the read-mapping process

    extract       Extract reads of a specific taxon from profiled results

    sam2bam       Convert GOTTCHA2 SAM to sorted/indexed BAM

    coverage-browser  Generate a coverage browser HTML from profiling results or a BAM

    download      Download the latest GOTTCHA2 database

    version       Display version information
    
Examples:
    gottcha2 profile -i reads.fastq -d database/db_prefix

    gottcha2 fast-profile -i reads.fastq -d database/db_prefix

    gottcha2 extract -d prefix.bam -d database/db_prefix -e 666

    gottcha2 download -d fast

    gottcha2 sam2bam -i prefix.sam -o prefix.bam

    gottcha2 coverage-browser -r results/ -o sample.coverage.html

For detailed help on a specific command:
    gottcha2 <command> --help
""")
    sys.exit(1)

def cli(args=None):
    """Parse command-line options and dispatch to the selected workflow."""
    args = sys.argv[1:] if args is None else list(args)
    if len(args) < 1:
        usage()
    elif args[0] in ("profile", "fast-profile", "extract"):
        profile.main(parse_args(__version__, args))
    elif args[0] == "download":
        download.main(args[1:])
    elif args[0] == "sam2bam":
        sam_to_bam.main(args[1:])
    elif args[0] == "coverage-browser":
        coverage_browser.main(args[1:])
    elif args[0] == "version":
        print(f"{__version__}")

    else:
        print(f"Error: '{args[0]}' is not a valid command")
        usage()


if __name__ == '__main__':
    cli()
