import math
import logging
import numpy as np
import pandas as pd
import sys
from . import taxonomy

def pile_lvl_zscore(tol_bp: int, tol_sig_len: int, linear_len: int) -> float:
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


def infer_sni_score(df: pd.DataFrame, error_rate: float) -> pd.DataFrame:
    """
    Estimate the Average Nucleotide Identity (SNI-score) together with 95% confidence intervals:
    - widens the interval when only a fraction of the signature space is actually covered ( SIG_COV )
    - project mismatches onto those unique positions
    - automatically becomes narrower as more signature bases are covered
    """

    df = df.copy()

    # from scipy.stats import norm
    # z = norm.ppf(0.5 + conf/2)  # ≈1.96
    z = 1.9599639845

    # use only unique covered signature bases
    n = df["COVERED_SIG_LEN"]
    cov = df["SIG_COV"].clip(lower=1e-12)  # avoid n_eff = 0

    # remove the expected sequencing-error penalty
    m_rate = (df["CONSENSUS_DIFF"]/df["COVERED_SIG_LEN"])
    m_rate_adj = m_rate - error_rate
    m_rate_adj = m_rate_adj.clip(lower=1e-12)  # avoid negative values

    # observed SCORE
    p_hat = 1 - m_rate_adj
    p = 1 - m_rate

    # cov is the coverage of the signature space, used to widen the confidence interval when only a fraction of the signature is covered
    n_eff  = n * cov

    z2     = z ** 2
    denom  = 1 + z2 / n_eff
    center = (p_hat + z2 / (2 * n_eff)) / denom
    hw     = (z * np.sqrt(
                (p_hat * (1 - p_hat)) / n_eff + z2 / (4 * n_eff ** 2)
             ) / denom)

    # observed-identity CI
    id_low, id_high = center - hw, center + hw

    # convert to true SCORE by adding the sequencing-error penalty
    score_low  = np.clip(id_low, 0, 1)
    score_high = np.clip(id_high, 0, 1)

    score_ci95 = "[" + score_low.round(6).astype(str) + "-" + score_high.round(6).astype(str) + "]"

    df = df.assign(
        SNI_SCORE    = center.round(6),
        SNI_CI95_LH  = score_ci95
    )

    return df


def group_refs_to_strains(ref_chunk_results: list, acc_list: list, acc_list_action: str, df_stats: pd.DataFrame) -> pd.DataFrame:
    """
    Group reference mapping results by strains and calculate strain-level statistics.
    
    Converts the mapping results dictionary to a pandas DataFrame and groups by
    taxonomic identifier. Calculates various statistics including total mapped bases,
    read counts, coverage, and depth of coverage.
    
    Parameters:
        ref_chunk_results (list): List of mapping statistics for each reference fragment chunks (output from parse_aln_from_bam)
        acc_list (list, optional): List of accessions of interest
        acc_list_action (str, optional): Action to take with the accession list (e.g., "exclude")
        df_stats (pandas.DataFrame): DataFrame containing genome signature statistics
    Returns:
        pandas.DataFrame: DataFrame with strain-level statistics
        int: Number of reads mapped to accessions of interest
    """
    # covert mapping info to df
    r_chunk_df = pd.DataFrame(ref_chunk_results[1:], columns=ref_chunk_results[0])

    # retrieve sig fragment info
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

    r_df[['ACC','RSTART','REND','TAXID','MISC']] = r_df['RNAME'].str.split('|', expand=True)

    logging.debug(f"Initial number of mapped reference fragments: {len(r_df)}")

    if acc_list:
        idx = (r_df['ACC'].isin(acc_list) | r_df['RNAME'].isin(acc_list))

        logging.debug(f"AOI list has {len(acc_list)} records. {idx.sum()} out of {len(r_df)} references match the accession list.")

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


def aggregate_taxonomy(str_df: pd.DataFrame, 
                       abu_col: str, 
                       tg_rank: str, 
                       mc: float, 
                       mr: int, 
                       ml: int, 
                       mz: float, 
                       sni_score_species: float, 
                       sni_score_strain: float, 
                       sni_score_cutoff: float, 
                       error_rate: float,
                       groups: dict) -> pd.DataFrame:
    """
    Aggregate strain-level results to higher taxonomic ranks.

    Starting from strain-level mapping data, this function rolls up statistics to
    higher taxonomic ranks (species, genus, family, etc.). It applies the specified
    cutoff criteria to filter results and marks entries that fall below these thresholds.

    Additional behaviors:
        1. For a species where all strains are supported only by strain-level signatures
           (SIG_LEVEL == 8), only strains passing the strain SNI_SCORE threshold are
           aggregated into species and higher ranks.
        2. Qualified species can be merged according to the provided groups dictionary.
           Each group is collapsed into its most abundant qualified member, with the
           relative abundance contribution of each merged member recorded in NOTE.

    Parameters:
        str_df (pandas.DataFrame): DataFrame with genome-level mapping stats
        abu_col (str): Column name to use for abundance calculations
        tg_rank (str): Target taxonomic rank
        mc (float): Minimum linear coverage threshold
        mr (int): Minimum read count threshold
        ml (int): Minimum covered signature length threshold
        mz (float): Maximum Z-score threshold (0 to disable)
        sni_score_cutoff (float): SNI-score cutoff for all levels
        sni_score_species (float): SNI-score cutoff for species level
        sni_score_strain (float): SNI-score cutoff for strain level
        error_rate (float): Error rate for SNI-score inference
        groups (dict): Species groups to merge

    Returns:
        pandas.DataFrame: DataFrame with rolled-up taxonomy at all ranks
    """

    major_ranks = {"superkingdom":1,"phylum":2,"class":3,"order":4,"family":5,"genus":6,"species":7,"strain":8}

    groups = groups or {}
    str_df = str_df.copy()

    # total reads mapped to accession#s of interest
    total_aoi_read_count = str_df['AOI_READ_COUNT'].sum()

    # produce columns for the final report at each ranks
    rep_df = pd.DataFrame()

    # add taxonomic lineage info
    ranks = list(major_ranks.keys())[::-1]

    def get_taxid_lineage(taxid):
        """get taxid lineage with {rank}_names and {rank}_taxids"""
        lineage = taxonomy.taxid2lineageDICT(taxid).values()
        return [d['name'] for d in lineage]+[d['taxid'] for d in lineage]

    def join_notes(x):
        """combine notes and remove duplicated entries"""
        notes = []
        for n in x.dropna():
            for note in str(n).split(';'):
                note = note.strip()
                if note and note not in notes:
                    notes.append(note)
        return '; '.join(notes)

    def taxid_key(x):
        """normalize taxid to string for group matching"""
        if pd.isna(x):
            return None
        if isinstance(x, float) and x.is_integer():
            return str(int(x))
        return str(x)

    try:
        cols = [f'{r}_name' for r in ranks]+[f'{r}_taxid' for r in ranks]
        str_df[cols] = str_df['TAXID'].map(get_taxid_lineage).to_list()
    except Exception as e:
        logging.error(f"Error processing rank {ranks}: {e}. Please verify that your taxonomy file matches the expected database.")
        sys.exit(1)
    
    logging.debug(f"Taxonomic lineage info added to {len(str_df)} strains.")

    # decide top signature level, convert the rank to the corresponding number
    str_df['SIG_LEVEL'] = str_df['SIG_LEVEL'].map(major_ranks)

    # identify species where all strains have only strain-level signatures
    str_df['_ALL_STRAIN_SIG'] = str_df.groupby('species_name')['SIG_LEVEL'].transform(
        lambda x: (x == major_ranks['strain']).all()
    )

    # infer the SNI-score for each strain
    idx = str_df['COVERED_SIG_LEN'] > 0
    str_df = str_df[idx].reset_index(drop=True)
    str_df["SIG_COV"] = str_df["COVERED_SIG_LEN"]/str_df["TOTAL_SIG_LEN"]
    str_df = infer_sni_score(str_df, error_rate)

    logging.debug(f"SNI-score inferred for {len(str_df[str_df['SIG_COV']>0])} strains.")

    if 'NOTE' not in str_df.columns:
        str_df['NOTE'] = ''
    else:
        str_df['NOTE'] = str_df['NOTE'].fillna('')

    # For species where all strains use strain-level signatures, only strains
    # passing strain SNI are allowed to aggregate into species and higher ranks.
    str_df['_ROLLUP'] = (~str_df['_ALL_STRAIN_SIG']) | (str_df['SNI_SCORE'] >= sni_score_strain)

    filtered = str_df['_ALL_STRAIN_SIG'] & (str_df['SNI_SCORE'] < sni_score_strain)
    str_df.loc[filtered, 'NOTE'] += (
        f"Excluded from higher-rank aggregation "
        f"(strain SNI_SCORE threshold {sni_score_strain} > "
        + str_df.loc[filtered, 'SNI_SCORE'].astype(str)
        + "); "
    )

    # keep all strains for strain-level reporting, but use only qualified strains
    # when aggregating to species and higher ranks
    rollup_df = str_df[str_df['_ROLLUP']].copy()

    total_abundance_strain = str_df[abu_col].sum()
    total_abundance_rollup = rollup_df[abu_col].sum()

    # check taxids are not assigned to multiple merge groups
    seen_taxids = {}
    for group_name, taxids in groups.items():
        for taxid in taxids:
            key = taxid_key(taxid)
            if key in seen_taxids:
                raise ValueError(
                    f"Taxid {taxid} occurs in both {seen_taxids[key]} and {group_name}."
                )
            seen_taxids[key] = group_name

    # iterate through ranks to get index and value
    for idx, rank in enumerate(ranks):

        # strain-level report uses all strains; higher ranks use rollup-qualified strains
        if rank == 'strain':
            tmp_df = str_df.copy()
        else:
            tmp_df = rollup_df.copy()

        tmp_df['LEVEL'] = rank
        tmp_df[['LVL_NAME', 'LVL_TAXID']] = tmp_df[[f'{rank}_name', f'{rank}_taxid']]

        if rank == 'superkingdom':
            tmp_df[['PARENT_NAME', 'PARENT_TAXID']] = ['root', '1']
        else:
            tmp_df[['PARENT_NAME', 'PARENT_TAXID']] = tmp_df[[f'{ranks[idx+1]}_name', f'{ranks[idx+1]}_taxid']]

        # rollup strains that make cutoffs
        lvl_df = None

        if rank == 'strain':
            lvl_df = tmp_df.copy()

        else:
            lvl_df = tmp_df.groupby('LVL_NAME').agg({
                'LEVEL':'first',
                'LVL_TAXID':'first',
                'PARENT_NAME':'first',
                'PARENT_TAXID':'first',
                'TOTAL_BP_MAPPED': 'sum',
                'READ_COUNT': 'sum',
                'TOTAL_BP_MISMATCH': 'sum',
                'TOTAL_BP_INDEL': 'sum',
                'TOTAL_READ_LEN': 'sum',
                'COVERED_SIG_LEN': 'sum',
                'MAPPED_SIG_LEN': 'sum',
                'TOTAL_SIG_LEN': 'sum',
                'CONSENSUS_DIFF': 'sum',
                'DEPTH': 'sum',
                'AOI_READ_COUNT': 'sum',
                'BEST_SIG_COV': 'max',
                'ZSCORE': 'min',
                'GENOMIC_CONTENT_EST': 'sum',
                'SIG_LEVEL': 'max',
                'GENOME_COUNT': 'count',
                'GENOME_SIZE': 'sum',
                'SNI_SCORE': 'max',
                'NOTE': join_notes
            })

            # find the index of the row with max SNI_SCORE in each group
            # pull out the low/high bounds from those rows
            max_idx = tmp_df.groupby('LVL_NAME')['SNI_SCORE'].idxmax()
            score_bounds = (tmp_df
                            .loc[max_idx, ['LVL_NAME', 'SNI_CI95_LH']]
                            .set_index('LVL_NAME'))
            lvl_df = lvl_df.join(score_bounds).reset_index()

        # merge qualified species according to groups
        if rank == 'species' and groups:

            # qualified species must pass all species-level filtering criteria
            qualified = (
                (lvl_df['SIG_LEVEL'] >= major_ranks['species']) &
                (lvl_df['SNI_SCORE'] >= sni_score_species) &
                (lvl_df['BEST_SIG_COV'] >= mc) &
                (lvl_df['READ_COUNT'] >= mr) &
                (lvl_df['COVERED_SIG_LEN'] >= ml)
            )

            if mz > 0:
                qualified &= lvl_df['ZSCORE'] <= mz

            sum_cols = [
                'TOTAL_BP_MAPPED',
                'READ_COUNT',
                'TOTAL_BP_MISMATCH',
                'TOTAL_BP_INDEL',
                'TOTAL_READ_LEN',
                'COVERED_SIG_LEN',
                'MAPPED_SIG_LEN',
                'TOTAL_SIG_LEN',
                'CONSENSUS_DIFF',
                'DEPTH',
                'AOI_READ_COUNT',
                'GENOMIC_CONTENT_EST',
                'GENOME_COUNT',
                'GENOME_SIZE'
            ]

            for group_name, taxids in groups.items():

                group_taxids = {taxid_key(x) for x in taxids}

                idx_group = qualified & lvl_df['LVL_TAXID'].map(taxid_key).isin(group_taxids)
                members = lvl_df[idx_group].copy()

                # only merge if at least two qualified species are present
                if len(members) < 2:
                    continue

                # representative = most abundant qualified species
                rep_idx = members[abu_col].idxmax()
                rep = lvl_df.loc[rep_idx].copy()

                group_abundance = members[abu_col].sum()

                fractions = []
                for _, row in members.sort_values(abu_col, ascending=False).iterrows():
                    frac = row[abu_col]/group_abundance if group_abundance > 0 else 0
                    fractions.append(
                        f"{row['LVL_NAME']} (taxid={row['LVL_TAXID']}): {frac:.4f}"
                    )

                # sum statistics of all merged species
                for col in sum_cols:
                    if col in lvl_df.columns:
                        rep[col] = members[col].sum()

                rep['BEST_SIG_COV'] = members['BEST_SIG_COV'].max()
                rep['ZSCORE'] = members['ZSCORE'].min()
                rep['SIG_LEVEL'] = members['SIG_LEVEL'].max()

                # if abundance column is not one of the summed columns,
                # explicitly sum it as well
                if abu_col not in sum_cols and abu_col not in [
                    'BEST_SIG_COV', 'ZSCORE', 'SIG_LEVEL', 'SNI_SCORE'
                ]:
                    rep[abu_col] = members[abu_col].sum()

                # keep representative SNI_SCORE and SNI_CI95_LH
                note = join_notes(members['NOTE'])
                if note:
                    note += '; '

                note += (
                    f"Merged {group_name} into {rep['LVL_NAME']} "
                    f"(taxid={rep['LVL_TAXID']}); "
                    f"abundance fractions: {', '.join(fractions)}"
                )
                rep['NOTE'] = note

                # replace group members with representative
                lvl_df.loc[rep_idx] = rep
                lvl_df = lvl_df.drop(members.index.difference([rep_idx])).reset_index(drop=True)

                # update qualified mask because dataframe indices changed
                qualified = (
                    (lvl_df['SIG_LEVEL'] >= major_ranks['species']) &
                    (lvl_df['SNI_SCORE'] >= sni_score_species) &
                    (lvl_df['BEST_SIG_COV'] >= mc) &
                    (lvl_df['READ_COUNT'] >= mr) &
                    (lvl_df['COVERED_SIG_LEN'] >= ml)
                )

                if mz > 0:
                    qualified &= lvl_df['ZSCORE'] <= mz

        # calculate the relative abundance of each taxon
        # strain level uses all strain abundance; higher ranks use rollup abundance
        if rank == 'strain':
            total_abundance = total_abundance_strain
        else:
            total_abundance = total_abundance_rollup

        lvl_df['ABUNDANCE'] = lvl_df[abu_col]

        if total_abundance > 0:
            lvl_df['REL_ABUNDANCE'] = lvl_df[abu_col]/total_abundance
        else:
            lvl_df['REL_ABUNDANCE'] = 0

        lvl_df['ABUNDANCE_DEPTH'] = lvl_df['DEPTH']
        if lvl_df['DEPTH'].sum() > 0:
            lvl_df['REL_ABUNDANCE_DEPTH'] = lvl_df['ABUNDANCE_DEPTH']/lvl_df['ABUNDANCE_DEPTH'].sum()
        else:
            lvl_df['REL_ABUNDANCE_DEPTH'] = 0

        lvl_df['ABUNDANCE_GC'] = lvl_df['GENOMIC_CONTENT_EST']
        if lvl_df['GENOMIC_CONTENT_EST'].sum() > 0:
            lvl_df['REL_ABUNDANCE_GC'] = lvl_df['GENOMIC_CONTENT_EST']/lvl_df['GENOMIC_CONTENT_EST'].sum()
        else:
            lvl_df['REL_ABUNDANCE_GC'] = 0

        # if 'NOTE' is not empty, add '; ' to the end of the string
        lvl_df['NOTE'] = lvl_df['NOTE'].fillna('')
        lvl_df['NOTE'] = lvl_df['NOTE'].apply(lambda x: f'{x}; ' if x and not x.endswith('; ') else x)

        # A note is added if a taxa has a rank with higher resolution than the target rank
        idx_biased = lvl_df['SIG_LEVEL'] < major_ranks[rank]
        lvl_df.loc[idx_biased, 'NOTE'] += f"Not shown ({rank}-result could be biased); "

        # add SCORE reason
        if rank == 'strain':
            filtered = (lvl_df['SNI_SCORE'] < sni_score_strain)
            lvl_df.loc[filtered, 'NOTE'] += (
                f"Filtered out (strain SNI_SCORE threshold {sni_score_strain} > "
                + lvl_df.loc[filtered, 'SNI_SCORE'].astype(str) + "); "
            )

        elif rank == 'species':
            filtered = (lvl_df['SNI_SCORE'] < sni_score_species)
            lvl_df.loc[filtered, 'NOTE'] += (
                f"Filtered out (species SNI_SCORE threshold {sni_score_species} > "
                + lvl_df.loc[filtered, 'SNI_SCORE'].astype(str) + "); "
            )

        else:
            filtered = (lvl_df['SNI_SCORE'] < sni_score_cutoff)
            lvl_df.loc[filtered, 'NOTE'] += (
                f"Filtered out (SNI_SCORE threshold {sni_score_cutoff} > "
                + lvl_df.loc[filtered, 'SNI_SCORE'].astype(str) + "); "
            )

        # add filtered reason
        filtered = (lvl_df['BEST_SIG_COV'] < mc)
        lvl_df.loc[filtered, 'NOTE'] += (
            f"Filtered out (minCov threshold {mc} > "
            + lvl_df.loc[filtered, 'BEST_SIG_COV'].astype(str) + "); "
        )

        filtered = (lvl_df['READ_COUNT'] < mr)
        lvl_df.loc[filtered, 'NOTE'] += (
            f"Filtered out (minReads threshold {mr} > "
            + lvl_df.loc[filtered, 'READ_COUNT'].astype(str) + "); "
        )

        filtered = (lvl_df['COVERED_SIG_LEN'] < ml)
        lvl_df.loc[filtered, 'NOTE'] += (
            f"Filtered out (minLen threshold {ml} > "
            + lvl_df.loc[filtered, 'COVERED_SIG_LEN'].astype(str) + "); "
        )

        if mz > 0:
            filtered = (lvl_df['ZSCORE'] > mz)
            lvl_df.loc[filtered, 'NOTE'] += (
                f"Filtered out (maxZscore threshold {mz} < "
                + lvl_df.loc[filtered, 'ZSCORE'].astype(str) + "); "
            )

        # remove internal columns
        lvl_df.drop(columns=['_ALL_STRAIN_SIG', '_ROLLUP'], errors='ignore', inplace=True)

        # concart ranks-dataframe to the report-dataframe
        rep_df = pd.concat(
            [lvl_df.sort_values('ABUNDANCE', ascending=False), rep_df],
            ignore_index=True
        )

    # add additional columns
    rep_df = rep_df.assign(
        SIG_COV                = rep_df["COVERED_SIG_LEN"]/rep_df["TOTAL_SIG_LEN"],
        READ_WT_SNI            = 1-(rep_df["TOTAL_BP_MISMATCH"]/rep_df["TOTAL_READ_LEN"]),
        CONSENSUS_SEQ_SNI      = 1-(rep_df['CONSENSUS_DIFF']/rep_df["COVERED_SIG_LEN"]),
        COVERED_MAPPED_SIG_COV = rep_df["COVERED_SIG_LEN"]/rep_df["MAPPED_SIG_LEN"],
        COVERED_SIG_DEPTH      = rep_df["TOTAL_BP_MAPPED"]/rep_df["COVERED_SIG_LEN"],
    )

    rep_df.drop(columns=['TAXID'], inplace=True)
    rep_df.rename(columns={"LVL_NAME": "NAME", "LVL_TAXID": "TAXID"}, inplace=True)

    logging.debug(f'rep_df:\n{rep_df}')

    return rep_df, total_aoi_read_count