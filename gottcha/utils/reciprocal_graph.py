from pathlib import Path
import logging
from typing import List, Tuple, Hashable, Optional
import pandas as pd
import numpy as np
import logging
from statistics import NormalDist
from scipy.sparse import coo_matrix, triu
from scipy.sparse.csgraph import connected_components
from . import taxonomy as t

def _wilson_lower_bound(
    successes: np.ndarray,
    trials: np.ndarray,
    confidence: float,
) -> np.ndarray:
    """Vectorized lower Wilson bound for binomial proportions."""
    if not 0.0 < confidence < 1.0:
        raise ValueError("confidence must be between 0 and 1")

    z = NormalDist().inv_cdf(0.5 + confidence / 2.0)
    successes = successes.astype(np.float64, copy=False)
    trials = trials.astype(np.float64, copy=False)
    p = successes / trials
    z2 = z * z
    denominator = 1.0 + z2 / trials

    return (
        p
        + z2 / (2.0 * trials)
        - z * np.sqrt((p * (1.0 - p) + z2 / (4.0 * trials)) / trials)
    ) / denominator


def get_species_hit_groups(
    df: pd.DataFrame,
    include_singletons: bool = True,
    same_genus_only: bool = True,
    genus_col: str = "GENUS_TAXID",
    min_reciprocal_reads: int = 3,
    min_reciprocal_fraction: float = 0.005,
    fraction_confidence: Optional[float] = 0.95,
    min_relative_strength: float = 0.0,
    return_edges: bool = False,
):
    """
    Group species connected by sufficiently supported reciprocal mappings.

    For each read, the species of its first alignment is linked to every
    other distinct species observed for that read. A species pair can form
    a grouping edge only when the relationship is observed in both
    directions and passes absolute and normalized support thresholds.

    Parameters
    ----------
    df
        DataFrame containing QNAME and SPECIES_TAXID. Row order within each
        QNAME defines alignment order. QNAME blocks need not be contiguous.

    include_singletons
        Return species without retained grouping edges as one-species groups.

    same_genus_only
        Permit edges only between species assigned to the same genus.

    genus_col
        Column containing genus names or taxids.

    min_reciprocal_reads
        Minimum number of independent reads required in EACH direction.
        For A--B, both count(A -> B) and count(B -> A) must meet this value.

    min_reciprocal_fraction
        Minimum directional association in EACH direction. For A--B, this
        is applied to count(A -> B)/number_of_A_source_reads and the reverse.
        When fraction_confidence is not None, the Wilson lower confidence
        bound is used instead of the raw fraction.

    fraction_confidence
        Confidence level for the Wilson lower bound. Use None to threshold
        raw directional fractions instead.

    min_relative_strength
        Optional local bridge pruning in [0, 1]. After absolute filtering,
        an edge must have strength at least this fraction of the strongest
        retained edge incident on EACH endpoint. Set to 0 to disable.

    return_edges
        If True, return (groups, edge_table). edge_table contains support,
        rates, confidence bounds, edge strengths, and filtering decisions.

    Returns
    -------
    dict, or (dict, pandas.DataFrame)
        Groups keyed as group1, group2, ...; optionally with diagnostics.
    """
    required = {"QNAME", "SPECIES_TAXID"}
    if same_genus_only:
        required.add(genus_col)

    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"Missing required column(s): {sorted(missing)}")
    if min_reciprocal_reads < 1:
        raise ValueError("min_reciprocal_reads must be >= 1")
    if not 0.0 <= min_reciprocal_fraction <= 1.0:
        raise ValueError("min_reciprocal_fraction must be in [0, 1]")
    if not 0.0 <= min_relative_strength <= 1.0:
        raise ValueError("min_relative_strength must be in [0, 1]")

    qcode, _ = pd.factorize(df["QNAME"], sort=False)
    scode, species = pd.factorize(df["SPECIES_TAXID"], sort=False)
    valid = (qcode >= 0) & (scode >= 0)

    qcode = qcode[valid].astype(np.int64, copy=False)
    scode = scode[valid].astype(np.int64, copy=False)
    n_species = len(species)

    empty_columns = [
        "species_a", "species_b", "a_to_b_reads", "b_to_a_reads",
        "a_source_reads", "b_source_reads", "a_to_b_fraction",
        "b_to_a_fraction", "a_to_b_lower", "b_to_a_lower",
        "edge_strength", "passes_read_threshold",
        "passes_fraction_threshold", "passes_relative_threshold", "kept",
    ]

    if n_species == 0:
        result = {}
        if return_edges:
            return result, pd.DataFrame(columns=empty_columns)
        return result

    # Map each species to a single genus code before reordering rows.
    species_genus = None
    if same_genus_only:
        gcode, _ = pd.factorize(df.loc[valid, genus_col], sort=False)
        genus_map_df = pd.DataFrame({"species": scode, "genus": gcode})
        observed = genus_map_df[genus_map_df["genus"] >= 0]

        genus_counts = observed.groupby("species", sort=False)["genus"].nunique()
        inconsistent = genus_counts[genus_counts > 1]
        if not inconsistent.empty:
            bad_species = species[inconsistent.index.to_numpy()].tolist()
            raise ValueError(
                "Each species must map to only one genus; inconsistent species: "
                f"{bad_species}"
            )

        species_genus = np.full(n_species, -1, dtype=np.int64)
        first_genus = observed.drop_duplicates("species", keep="first")
        species_genus[first_genus["species"].to_numpy()] = (
            first_genus["genus"].to_numpy()
        )

    # Stable grouping by QNAME preserves the original alignment order within
    # every read even when QNAME blocks are not contiguous in the input.
    order = np.argsort(qcode, kind="stable")
    qcode = qcode[order]
    scode = scode[order]

    # Keep the first occurrence of each read--species combination. This
    # prevents multiple reference alignments from the same species from
    # inflating an edge count.
    pair_key = qcode * n_species + scode
    unique_hit = ~pd.Series(pair_key).duplicated(keep="first").to_numpy()
    qcode = qcode[unique_hit]
    scode = scode[unique_hit]

    first_hit = np.r_[True, qcode[1:] != qcode[:-1]]
    first_idx = np.flatnonzero(first_hit)
    group_sizes = np.diff(np.r_[first_idx, len(qcode)])
    source_for_hit = np.repeat(scode[first_idx], group_sizes)

    source_counts = np.bincount(
        scode[first_idx],
        minlength=n_species,
    ).astype(np.int64, copy=False)

    edge_rows = ~first_hit
    src = source_for_hit[edge_rows]
    dst = scode[edge_rows]

    if same_genus_only and len(src):
        src_genus = species_genus[src]
        dst_genus = species_genus[dst]
        genus_ok = (src_genus >= 0) & (src_genus == dst_genus)
        src = src[genus_ok]
        dst = dst[genus_ok]

    directed = coo_matrix(
        (
            np.ones(len(src), dtype=np.int64),
            (src, dst),
        ),
        shape=(n_species, n_species),
        dtype=np.int64,
    ).tocsr()
    directed.sum_duplicates()
    directed.setdiag(0)
    directed.eliminate_zeros()

    # The minimum matrix is nonzero only for pairs observed in both directions.
    reciprocal = triu(
        directed.minimum(directed.T),
        k=1,
        format="coo",
    )

    a = reciprocal.row.astype(np.int64, copy=False)
    b = reciprocal.col.astype(np.int64, copy=False)

    if len(a):
        a_to_b = directed[a, b].A1.astype(np.int64, copy=False)
        b_to_a = directed[b, a].A1.astype(np.int64, copy=False)
        n_a = source_counts[a]
        n_b = source_counts[b]

        frac_ab = a_to_b / n_a
        frac_ba = b_to_a / n_b

        if fraction_confidence is None:
            lower_ab = frac_ab.copy()
            lower_ba = frac_ba.copy()
        else:
            lower_ab = _wilson_lower_bound(a_to_b, n_a, fraction_confidence)
            lower_ba = _wilson_lower_bound(b_to_a, n_b, fraction_confidence)

        strength = np.minimum(lower_ab, lower_ba)
        passes_reads = (
            (a_to_b >= min_reciprocal_reads)
            & (b_to_a >= min_reciprocal_reads)
        )
        passes_fraction = strength >= min_reciprocal_fraction
        base_keep = passes_reads & passes_fraction

        passes_relative = np.ones(len(a), dtype=bool)
        if min_relative_strength > 0 and base_keep.any():
            strongest = np.zeros(n_species, dtype=np.float64)
            np.maximum.at(strongest, a[base_keep], strength[base_keep])
            np.maximum.at(strongest, b[base_keep], strength[base_keep])
            passes_relative = (
                (strength >= min_relative_strength * strongest[a])
                & (strength >= min_relative_strength * strongest[b])
            )

        kept = base_keep & passes_relative

        edge_table = pd.DataFrame(
            {
                "species_a": species[a],
                "species_b": species[b],
                "a_to_b_reads": a_to_b,
                "b_to_a_reads": b_to_a,
                "a_source_reads": n_a,
                "b_source_reads": n_b,
                "a_to_b_fraction": frac_ab,
                "b_to_a_fraction": frac_ba,
                "a_to_b_lower": lower_ab,
                "b_to_a_lower": lower_ba,
                "edge_strength": strength,
                "passes_read_threshold": passes_reads,
                "passes_fraction_threshold": passes_fraction,
                "passes_relative_threshold": passes_relative,
                "kept": kept,
            }
        ).sort_values(
            ["kept", "edge_strength", "a_to_b_reads", "b_to_a_reads"],
            ascending=[False, False, False, False],
            ignore_index=True,
        )

        kept_a = a[kept]
        kept_b = b[kept]
    else:
        edge_table = pd.DataFrame(columns=empty_columns)
        kept_a = np.empty(0, dtype=np.int64)
        kept_b = np.empty(0, dtype=np.int64)

    grouping_graph = coo_matrix(
        (
            np.ones(2 * len(kept_a), dtype=np.uint8),
            (
                np.r_[kept_a, kept_b],
                np.r_[kept_b, kept_a],
            ),
        ),
        shape=(n_species, n_species),
    ).tocsr()

    _, labels = connected_components(
        grouping_graph,
        directed=False,
        return_labels=True,
    )

    raw_groups: dict[int, list[Hashable]] = {}
    for sp, label in zip(species, labels):
        raw_groups.setdefault(int(label), []).append(sp)

    if not include_singletons:
        raw_groups = {
            label: members
            for label, members in raw_groups.items()
            if len(members) > 1
        }

    groups = {
        f"group{i + 1}": members
        for i, members in enumerate(raw_groups.values())
    }

    if return_edges:
        return groups, edge_table
    return groups

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

    logging.info(f'Grouping species using reciprocal directed read-mapping relationships...')
    # Extract taxids from the REF column and map them to species and genus taxids.
    taxids = df["REF"].str.rsplit("|", n=2).str[-2]
    taxid_to_gen_taxid = {
        taxid: t.taxid2taxidOnRank(
            taxid,
            target_rank="genus"
        )
        for taxid in taxids.unique()
    }
    df["GENUS_TAXID"] = taxids.map(taxid_to_gen_taxid)

    # Group the alignments by species taxid to identify species-level hits.
    species_groups = get_species_hit_groups(df)

    return species_groups