import subprocess
import logging
from typing import List, Tuple
import pandas as pd
import logging
import pandas as pd

def minimap2(reads: List, db: str, threads: int, mm_options: str, presetx: str, samfile: str, logfile: str) -> Tuple[int, str, str]:
    """
    Map reads to the reference database using minimap2.

    Builds and executes a command to run minimap2 for read mapping, with parameters
    adjusted based on input settings. Filters the SAM output to keep only relevant
    alignments.

    Parameters:
        reads (List): List of input read file paths
        db (str): Path to the minimap2 database (without .mmi extension)
        threads (int): Number of threads to use
        mm_options (str): Minimap2 options for read mapping
        presetx (str): Minimap2 preset mode ('sr', 'map-pb', or 'map-ont')
        samfile (str): Output SAM file path
        logfile (str): Log file path
        nanopore (bool): Whether to use Nanopore-specific settings

    Returns:
        Tuple[int, str, str]: (
            exitcode (int): Exit code from the mapping process,
            cmd (str): Command that was executed,
            errs (str): Error output from the command
        )
    """
    input_file = " ".join(reads)

    # Minimap2 options for short reads: the options here is essentailly the -x 'sr' equivalent with some modifications on scoring
    sr_opts = f"-x sr {mm_options} -a -N20 --eqx --secondary=no --sam-hit-only"

    if presetx != 'sr':
        sr_opts = f"-x {presetx} -N20 --secondary=no --sam-hit-only -a"

    bash_cmd   = f"set -o pipefail; set -x;"
    mm2_cmd    = f"minimap2 {sr_opts} -t{threads} {db}.mmi {input_file}"
    filter_cmd = f"sed '/^@/d'"  # filter out header lines
    cmd        = f"{bash_cmd} {mm2_cmd} 2>> {logfile} | {filter_cmd} > {samfile}"

    logging.info(f"Readmapping command: {mm2_cmd}")

    proc = subprocess.Popen(cmd, shell=True, executable='/bin/bash', stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    outs, errs = proc.communicate()
    exitcode = proc.poll()

    return exitcode, mm2_cmd, errs


def post_processing_sam(samfile: str, samfile_temp: str) -> bool:
    """
    Removing multiple hits from the SAM file by keeping only the best alignment for each read.

    Parameters:
        samfile (str): Path to the SAM file
        samfile_temp (str): Path to the temporary SAM file with only the best alignments

    Returns:
        bool: False if no multiple hits were found, True if multiple hits were removed
    """
    logging.info(f'Loading the SAM file...')

    df = pd.read_csv(samfile,
                     sep='\t',
                     header=None,
                     usecols=[0, 1, 13],
                     names=['QNAME', 'FLAG', 'AS'],
                     converters={
                         'AS': lambda x: x.replace('AS:i:', '')
                     },
                     dtype={'QNAME': 'str', 'FLAG': 'uint16'}
    )

    aln_count = len(df)
    logging.info(f'Total alignments in SAM file: {aln_count}')

    df[['AS']] = df[['AS']].astype('int16')

    logging.info(f'Filtering non-primary hits...')
    # for each row, if the flag bitwise AND with 256 (not primary alignment) or 2048 (supplementary), then remove them from the df
    df = df[~(df['FLAG'] & (256|2048)).astype(bool)]
    logging.info(f'After removing non-primary hits: {len(df)}')

    logging.info(f'Identifying top score hits...')
    # if FLAG bitwise AND with 128 (second in pair), append '/2' to the QNAME
    idx = (df['FLAG'] & 128).astype(bool)
    df.loc[idx, 'QNAME'] = df.loc[idx, 'QNAME'] + '/2'

    # get the index with the best alignment score for each read
    idxmax = df.groupby('QNAME')['AS'].idxmax()
    logging.info(f'Total top score hits: {len(idxmax)}')

    if len(idxmax) == aln_count:
        logging.info(f'No multiple hits found. Keeping the original SAM file.')
        return False, aln_count, aln_count
    else:
        # Create a set of indices for faster lookup
        idxmax_set = set(idxmax.values)
        del idxmax

        logging.info(f'Writing top score hits...')
        with open(samfile_temp, 'w') as fout, open(samfile, 'r') as fin:
            for idx, line in enumerate(fin):
                if idx > 0 and idx % 100000 == 0:
                    logging.debug(f'Processed {idx} lines...')

                if idx in idxmax_set:
                    fout.write(line)
        logging.info(f'{len(idxmax_set)} hits written to {samfile_temp}.')

        return True, aln_count, len(idxmax_set)
