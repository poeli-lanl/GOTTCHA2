from pathlib import Path
import logging
from typing import List, Tuple
import pandas as pd
import numpy as np
import logging
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from . import taxonomy as t

def _parse_cigar_query_interval(cigar):
    """
    Convert a CIGAR into query coordinates.

    Returns
    -------
    qstart : int
        Start of aligned portion in the query.
    qend : int
        End of aligned portion in the query, half-open [qstart, qend).
    qlen : int
        Original read length inferred from CIGAR, including hard clips.

    Query-consuming:
        M, I, S, =, X

    Hard clipping H is included in the original read coordinate system.

    Reference-only operations:
        D, N

    Neither:
        P
    """
    if not cigar or cigar == "*":
        return -1, -1, -1

    ops = []
    n = 0
    have_number = False

    for ch in cigar:
        if "0" <= ch <= "9":
            n = n * 10 + (ord(ch) - 48)
            have_number = True
        else:
            if not have_number:
                raise ValueError(f"Invalid CIGAR: {cigar}")

            ops.append((n, ch))
            n = 0
            have_number = False

    if have_number:
        raise ValueError(f"Invalid CIGAR: {cigar}")

    # Original query length.
    # H does not appear in SEQ but is part of the original read.
    qlen = sum(
        length
        for length, op in ops
        if op in "MISH=X"
    )

    # Leading S/H = query margin before alignment
    left = 0
    i = 0
    while i < len(ops) and ops[i][1] in "SH":
        left += ops[i][0]
        i += 1

    # Trailing S/H = query margin after alignment
    right = 0
    i = len(ops) - 1
    while i >= 0 and ops[i][1] in "SH":
        right += ops[i][0]
        i -= 1

    return left, qlen - right, qlen


def _add_query_intervals(df, flag_col="FLAG"):
    """
    Add QSTART, QEND, and READ_LEN.

    If FLAG exists, reverse-strand alignments are converted into the
    original read coordinate system.
    """
    n = len(df)

    qstart = np.full(n, -1, dtype=np.int32)
    qend = np.full(n, -1, dtype=np.int32)
    qlen = np.full(n, -1, dtype=np.int32)

    cigars = df["CIGAR"].to_numpy()

    if flag_col is not None and flag_col in df.columns:
        flags = df[flag_col].to_numpy(dtype=np.int64, copy=False)
    else:
        flags = None

    for i, cigar in enumerate(cigars):
        start, end, length = _parse_cigar_query_interval(cigar)

        if length < 0:
            continue

        # SAM FLAG 0x10 = query mapped to reverse strand.
        #
        # CIGAR is oriented relative to the aligned sequence, so mirror
        # the query interval to recover coordinates in the original read.
        if flags is not None and (int(flags[i]) & 0x10):
            start, end = length - end, length - start

        qstart[i] = start
        qend[i] = end
        qlen[i] = length

    df = df.copy()

    df["QSTART"] = qstart
    df["QEND"] = qend
    # df["READ_LEN"] = qlen

    return df

def _filter_overlapping_alignments(
    df,
    min_overlap_bp=1,
    min_overlap_frac=0.0,
):
    """
    Mark each alignment if it overlaps any earlier/higher-scoring
    alignment for the same read.

    Assumptions
    -----------
    df is sorted by:
        1. QNAME
        2. alignment score descending

    Parameters
    ----------
    min_overlap_bp : int
        Minimum number of overlapping query bases.

    min_overlap_frac : float
        Minimum fraction of the current alignment's query span that
        must overlap higher-scoring alignments.

    Returns
    -------
    DataFrame with:
        OVERLAP_BP
        OVERLAP_FRAC
        OVERLAPS_HIGHER
    """
    names = df["QNAME"].to_numpy()
    starts = df["QSTART"].to_numpy(dtype=np.int64, copy=False)
    ends = df["QEND"].to_numpy(dtype=np.int64, copy=False)

    n = len(df)

    overlap_bp = np.zeros(n, dtype=np.int32)
    overlap_frac = np.zeros(n, dtype=np.float32)
    overlaps = np.zeros(n, dtype=bool)

    current_qname = None
    covered_mask = 0

    for i in range(n):
        qname = names[i]

        # New read: reset query coverage
        if i == 0 or qname != current_qname:
            current_qname = qname
            covered_mask = 0

        start = int(starts[i])
        end = int(ends[i])

        if start < 0 or end <= start:
            continue

        span = end - start

        # Example:
        # [20, 70) =>
        # bits 20..69 set to 1
        aln_mask = ((1 << span) - 1) << start

        # All previous records for this QNAME have higher/equal score.
        intersect = aln_mask & covered_mask

        if intersect:
            bp = intersect.bit_count()
            frac = bp / span

            overlap_bp[i] = bp
            overlap_frac[i] = frac

            if bp >= min_overlap_bp and frac >= min_overlap_frac:
                overlaps[i] = True

        # Add this alignment to higher-score coverage for later records.
        covered_mask |= aln_mask

    # return non overlaps
    df = df[overlaps == False].copy()

    # result["OVERLAP_BP"] = overlap_bp
    # result["OVERLAP_FRAC"] = overlap_frac
    # result["OVERLAPS_HIGHER"] = overlaps
    return df

def _get_species_hit_groups(df, include_singletons=True, same_genus_only=True):
    """
    Group species using reciprocal directed read-mapping relationships.

    For each read, the species of the FIRST alignment is treated as the
    source species and is linked to every other species observed for that
    read.

    Only reciprocal (bidirectional) species relationships are retained.

    Optionally, reciprocal relationships are only allowed between species
    belonging to the same genus.

    Parameters
    ----------
    df : pandas.DataFrame
        Must contain:
            QNAME
            SPECIES_TAXID

        If same_genus_only=True, must also contain `genus_col`.

        Rows must be in alignment order within each read. The first row
        for each QNAME is treated as that read's first alignment.

    include_singletons : bool, default=True
        If True, species without reciprocal relationships are returned
        as one-species groups.

    same_genus_only : bool, default=False
        If True, only species belonging to the same genus can be connected
        and therefore placed in the same group.

    genus_col : str, default="GENUS_TAXID"
        Column containing the genus identifier/name.

    Returns
    -------
    dict
        Example:
        {
            "group1": [species1, species2],
            "group2": [species3],
        }
    """

    # Validate input
    required = {"QNAME", "SPECIES_TAXID"}
    genus_col="GENUS_TAXID"

    if same_genus_only:
        required.add(genus_col)

    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"Missing required column(s): {sorted(missing)}"
        )

    # Encode reads and species while preserving dataframe order.
    qcodes, _ = pd.factorize(df["QNAME"], sort=False)
    species_codes, species = pd.factorize(df["SPECIES_TAXID"], sort=False)

    valid = (qcodes >= 0) & (species_codes >= 0)
    qcodes = qcodes[valid]
    species_codes = species_codes[valid]
    n_species = len(species)

    # Map each species code to its genus.
    species_to_genus = None

    if same_genus_only:
        tmp = df.loc[valid, ["SPECIES_TAXID", genus_col]].dropna()

        # Check that each species maps to only one genus.
        genus_counts = (
            tmp.groupby("SPECIES_TAXID")[genus_col]
            .nunique()
        )

        inconsistent = genus_counts[genus_counts > 1]

        if len(inconsistent):
            raise ValueError(
                "Some species are associated with multiple genera: "
                f"{inconsistent.index.tolist()}"
            )

        genus_lookup = (
            tmp.drop_duplicates("SPECIES_TAXID")
            .set_index("SPECIES_TAXID")[genus_col]
            .to_dict()
        )

        species_to_genus = np.asarray(
            [genus_lookup.get(sp) for sp in species],
            dtype=object,
        )

    # --------------------------------------------------------------
    # Find read boundaries.
    #
    # This assumes all rows belonging to a QNAME are contiguous.
    # --------------------------------------------------------------
    if len(qcodes) == 0:
        return {}

    starts = np.r_[0, np.flatnonzero(qcodes[1:] != qcodes[:-1]) + 1]

    ends = np.r_[starts[1:], len(qcodes)]

    # --------------------------------------------------------------
    # Build directed edges:
    #
    #     first species -> every other species
    #
    # Duplicate species within the same read are ignored.
    # --------------------------------------------------------------
    src_list = []
    dst_list = []

    for start, end in zip(starts, ends):

        hits = species_codes[start:end]

        if len(hits) < 2:
            continue

        first_species = hits[0]
        seen = {first_species}

        for target_species in hits[1:]:

            if target_species in seen:
                continue

            seen.add(target_species)

            # Ignore cross-genus edges if requested.
            if same_genus_only:
                source_genus = species_to_genus[first_species]
                target_genus = species_to_genus[target_species]

                # Species without genus information cannot form
                # genus-restricted edges.
                if (
                    source_genus is None
                    or target_genus is None
                    or source_genus != target_genus
                ):
                    continue

            src_list.append(first_species)
            dst_list.append(target_species)

    # If there are no cross-species edges, every species is isolated.
    if not src_list:

        if not include_singletons:
            return {}

        return {
            f"group{i + 1}": [sp]
            for i, sp in enumerate(species)
        }

    src = np.asarray(src_list, dtype=np.int64)
    dst = np.asarray(dst_list, dtype=np.int64)

    # Directed graph.
    directed_graph = coo_matrix(
        (
            np.ones(len(src), dtype=np.uint8),
            (src, dst),
        ),
        shape=(n_species, n_species),
    ).tocsr()

    # Collapse multiple supporting reads to edge presence/absence.
    directed_graph.data[:] = 1

    # --------------------------------------------------------------
    # Retain only reciprocal relationships:
    #
    #     A -> B
    #     B -> A
    #
    # becomes:
    #
    #     A <-> B
    # --------------------------------------------------------------
    reciprocal_graph = directed_graph.multiply(
        directed_graph.T
    )

    reciprocal_graph.eliminate_zeros()

    # --------------------------------------------------------------
    # Connected components of the reciprocal graph.
    #
    # Transitivity is allowed:
    #
    #     A <-> B <-> C
    #
    # gives:
    #
    #     {A, B, C}
    #
    # With same_genus_only=True, all members are guaranteed to belong
    # to the same genus because cross-genus edges were removed earlier.
    # --------------------------------------------------------------
    _, labels = connected_components(
        reciprocal_graph,
        directed=False,
        return_labels=True,
    )

    groups = {}

    for sp, label in zip(species, labels):
        groups.setdefault(label, []).append(sp)

    if not include_singletons:
        groups = {
            label: members
            for label, members in groups.items()
            if len(members) > 1
        }

    return {
        f"group{i + 1}": members
        for i, members in enumerate(groups.values())
    }


def reciprocal_relationships_from_sam(samfile: Path, min_alen: int) -> dict:
    """
    Extract reciprocal relationships from a SAM file.

    Parameters:
        samfile (Path): Path to the SAM file.

    Returns:
        dict: A dictionary of species-level reciprocal relationship groups.
    """
    logging.info(f'Loading the SAM file...')

    df = pd.read_csv(
        samfile,
        sep="\t",
        header=None,
        usecols=[0, 1, 2, 5, 13],
        names=["QNAME", "FLAG", "REF", "CIGAR", "AS"],
        dtype={
            "FLAG": "int16",
        },
        converters={
            'AS': lambda x: x.replace('AS:i:', '')
        }
    )
    df[['AS']] = df[['AS']].astype('int16')

    logging.debug(f'Loaded SAM file with {len(df)} alignments.')

    logging.info(f'Identifying reads having >1 distinct species....')
    taxids = df["REF"].str.rsplit("|", n=2).str[-2]

    taxid_to_spe_taxid = {
        taxid: t.taxid2taxidOnRank(
            taxid,
            target_rank="species"
        )
        for taxid in taxids.unique()
    }

    df["SPECIES_TAXID"] = taxids.map(taxid_to_spe_taxid)

    # Keep only QNAMEs having alignments to >1 distinct SPECIES_TAXID
    mask = (
        df.groupby("QNAME", sort=False)["SPECIES_TAXID"]
        .transform("nunique")
        .gt(1)
    )
    df = df.loc[mask].reset_index(drop=True)
    logging.debug(f'Filtered SAM file to {len(df)} alignments with reads having >1 distinct species.')

    # Filter alignments based on minimum identity and minimum alignment length.
    logging.info(f'Filtering out alignments not meeting alignment criteria...')
    df = df[df['AS'] > min_alen*2]
    logging.debug(f'Filtered SAM file to {len(df)} alignments meeting min_alen criteria.')

    # Mark overlapping alignments and filter out those that overlap higher-scoring alignments.
    logging.info(f'Filtering out alignments overlapping higher-scoring ones...')
    df = _add_query_intervals(df)
    df = _filter_overlapping_alignments(df, min_overlap_bp=0)
    logging.debug(f'Filtered SAM file to {len(df)} alignments overlapping higher-scoring ones.')

    logging.info(f'Grouping species using reciprocal directed read-mapping relationships...')
    # Extract taxids from the REF column and map them to species and genus taxids.
    taxid_to_gen_taxid = {
        taxid: t.taxid2taxidOnRank(
            taxid,
            target_rank="genus"
        )
        for taxid in taxids.unique()
    }
    df["GENUS_TAXID"] = taxids.map(taxid_to_gen_taxid)

    # Group the alignments by species taxid to identify species-level hits.
    species_groups = _get_species_hit_groups(df)

    return species_groups