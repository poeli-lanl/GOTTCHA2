#!/usr/bin/env python3
try:
    # Package usage (installed)
    from . import pull_database
    from .. import gottcha2
    from . import gottcha_sam_to_bam
    from . import gottcha2_bam
except ImportError:  # pragma: no cover
    # Script usage (running from source directory)
    import pull_database
    import gottcha.gottcha2 as gottcha2
    import gottcha_sam_to_bam
    import gottcha2_bam
import sys

def usage():
    """Display usage information for GOTTCHA2 command-line tool."""
    version = gottcha2.__version__
    print(f"""
GOTTCHA2 - Genomic Origin Through Taxonomic CHAllenge v{version}

Usage:
    gottcha2 <command> [options]

Commands:
    profile    Taxonomic profiling of metagenomic reads
              (Map reads to signature database and classify)

    extract    Extract reads of a specific taxon from profiled results

    sam2bam    Convert SAM to sorted/indexed BAM

    pull       Download/update GOTTCHA2 databases

    version    Display version information
    
Examples:
    gottcha2 profile -i reads.fastq -d database/db_prefix

    gottcha2 extract -s prefix.sam -d database/db_prefix -e 666

    gottcha2 sam2bam -i alignments.sam -o alignments.bam

For detailed help on a specific command:
    gottcha2 <command> --help
""")
    sys.exit(1)

def gottcha2_command():
    args = sys.argv[1:]
    if len(args) < 1:
        usage()
    elif args[0] == "profile":
        gottcha2.main(args[1:])
    elif args[0] == "profile_bam":
        gottcha2_bam.main(args[1:])
    elif args[0] == "pull":
        pull_database.main(args[1:])
    elif args[0] == "sam2bam":
        gottcha_sam_to_bam.main(args[1:])
    elif args[0] == "version":
        print(f"{gottcha2.__version__}")
    elif args[0] == "extract":
        gottcha2.main(args[1:])

    else:
        print(f"Error: '{args[0]}' is not a valid command")
        usage()
