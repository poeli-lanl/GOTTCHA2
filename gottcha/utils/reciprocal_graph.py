from pathlib import Path
import gzip
import logging
from statistics import NormalDist
from typing import Hashable, Optional

import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

from . import taxonomy as t


SAM_UNMAPPED = 0x4
SAM_FIRST_SEGMENT = 0x40
SAM_LAST_SEGMENT = 0x80
SAM_SECONDARY = 0x100
SAM_SUPPLEMENTARY = 0x800


def _wilson_lower_bound(
    successes: np.ndarray,
    trials: np.ndarray,
    confidence: float,
) -> np.ndarray:
    """Vectorized lower Wilson bound for binomial proportions."""
    if not 0.0 < confidence < 1.0:
        raise ValueError("confidence must be between 0 and 1")

    successes = np.asarray(successes, dtype=np.float64)
    trials = np.asarray(trials, dtype=np.float64)
    if np.any(trials <= 0):
        raise ValueError("trials must be > 0")

    z = NormalDist().inv_cdf(0.5 + confidence / 2.0)
    p = successes / trials
    z2 = z * z
    denominator = 1.0 + z2 / trials

    return (
        p
        + z2 / (2.0 * trials)
        - z * np.sqrt((p * (1.0 - p) + z2 / (4.0 * trials)) / trials)
    ) / denominator


def _mate_number(flags: pd.Series) -> np.ndarray:
    """Return 1/2 for paired segments and 0 for unpaired/unspecified reads."""
    values = flags.to_numpy(dtype=np.int64, copy=False)
    return np.where(
        (values & SAM_FIRST_SEGMENT) != 0,
        1,
        np.where((values & SAM_LAST_SEGMENT) != 0, 2, 0),
    ).astype(np.int8)


def _parse_ref_metadata(ref: str) -> tuple[str, Optional[int], Optional[int], str, str]:
    """
    Parse a GOTTCHA2 signature reference name.

    The current GOTTCHA2 convention places the taxonomy id in the penultimate
    pipe-delimited field.  When the reference also has the common layout

        sequence_accession|signature_start|signature_end|taxid|genome_id

    the signature interval is returned as well.  The interval fields are
    optional here so this function remains compatible with older databases.
    """
    text = str(ref)
    parts = text.rsplit("|", 4)

    taxid = ""
    genome_id = ""
    seq_id = text
    start = None
    end = None

    if len(parts) >= 2:
        taxid = parts[-2]
        genome_id = parts[-1]

    if len(parts) == 5:
        seq_id = parts[0]
        try:
            start = int(parts[1])
            end = int(parts[2])
        except (TypeError, ValueError):
            start = None
            end = None

    return seq_id, start, end, taxid, genome_id


def _signature_length_from_ref(ref: str) -> float:
    _, start, end, _, _ = _parse_ref_metadata(ref)
    if start is None or end is None:
        return np.nan
    return float(abs(end - start) + 1)


def _extract_as_tag(optional_fields: list[str]) -> Optional[int]:
    """Return AS:i from SAM optional fields without assuming a tag column."""
    for field in optional_fields:
        if field.startswith("AS:i:"):
            try:
                return int(field[5:])
            except ValueError:
                return None
    return None


def _open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def _load_sam_without_secondary(samfile: Path, min_alen: int) -> pd.DataFrame:
    """
    Load mapped primary and supplementary SAM records only.

    Secondary alignments (FLAG 0x100) are always discarded.  This preserves
    the behavior requested for this graph: relationships are derived only
    from discordant primary mates and primary/supplementary split alignments.

    ``min_alen`` intentionally retains the score filter used by the previous
    implementation: AS > 2 * min_alen.  AS is now located by tag name rather
    than by a fixed SAM column.
    """
    records = []
    secondary_skipped = 0
    unmapped_skipped = 0
    missing_as_skipped = 0
    low_score_skipped = 0

    with _open_text(Path(samfile)) as handle:
        for line in handle:
            if not line or line.startswith("@"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                continue

            try:
                flag = int(fields[1])
            except ValueError:
                continue

            # Secondary mappings are deliberately not part of this model.
            if flag & SAM_SECONDARY:
                secondary_skipped += 1
                continue

            if flag & SAM_UNMAPPED or fields[2] == "*":
                unmapped_skipped += 1
                continue

            as_score = _extract_as_tag(fields[11:])
            if as_score is None:
                missing_as_skipped += 1
                continue
            if as_score <= min_alen * 2:
                low_score_skipped += 1
                continue

            try:
                pos = int(fields[3])
                mapq = int(fields[4])
            except ValueError:
                continue

            records.append(
                (
                    fields[0],
                    flag,
                    fields[2],
                    pos,
                    mapq,
                    fields[5],
                    as_score,
                )
            )

    logging.debug(
        "SAM load: kept=%d secondary_skipped=%d unmapped_skipped=%d "
        "missing_AS_skipped=%d low_score_skipped=%d",
        len(records),
        secondary_skipped,
        unmapped_skipped,
        missing_as_skipped,
        low_score_skipped,
    )

    return pd.DataFrame(
        records,
        columns=["QNAME", "FLAG", "REF", "POS", "MAPQ", "CIGAR", "AS"],
    )


def _prepare_alignment_df(df: pd.DataFrame) -> pd.DataFrame:
    """Validate flags, remove secondary/unmapped rows, and annotate row type."""
    required = {"QNAME", "FLAG", "REF", "SPECIES_TAXID"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"Missing required column(s): {sorted(missing)}")

    out = df.copy()
    flags = pd.to_numeric(out["FLAG"], errors="coerce")
    out = out.loc[flags.notna()].copy()
    out["FLAG"] = flags.loc[out.index].astype(np.int64)

    keep = (
        ((out["FLAG"].to_numpy() & SAM_UNMAPPED) == 0)
        & ((out["FLAG"].to_numpy() & SAM_SECONDARY) == 0)
        & out["SPECIES_TAXID"].notna().to_numpy()
    )
    out = out.loc[keep].copy()

    out["IS_SUPPLEMENTARY"] = (
        out["FLAG"].to_numpy(dtype=np.int64, copy=False) & SAM_SUPPLEMENTARY
    ) != 0
    out["IS_PRIMARY"] = ~out["IS_SUPPLEMENTARY"]
    out["MATE"] = _mate_number(out["FLAG"])

    return out


def _best_primary_per_segment(df: pd.DataFrame) -> pd.DataFrame:
    """Keep one representative primary record for each QNAME/read segment."""
    primary = df.loc[df["IS_PRIMARY"]].copy()
    if primary.empty:
        return primary

    if "AS" in primary.columns:
        primary["AS"] = pd.to_numeric(primary["AS"], errors="coerce")
        primary = primary.sort_values(
            ["QNAME", "MATE", "AS"],
            ascending=[True, True, False],
            kind="stable",
        )
    else:
        primary = primary.sort_values(["QNAME", "MATE"], kind="stable")

    # SAM should contain a single primary line for each segment.  Keeping the
    # best AS makes the function robust to malformed or merged input files.
    return primary.drop_duplicates(["QNAME", "MATE"], keep="first")


def _canonicalize_events(
    frame: pd.DataFrame,
    left_code: str,
    right_code: str,
    left_ref: str,
    right_ref: str,
    evidence: str,
    left_is_primary: bool,
    right_is_primary: bool,
) -> pd.DataFrame:
    """Canonicalize species pairs while preserving which endpoint is primary."""
    if frame.empty:
        return pd.DataFrame(
            columns=[
                "QNAME",
                "A_CODE",
                "B_CODE",
                "A_REF",
                "B_REF",
                "HAS_MATE",
                "HAS_SPLIT",
                "A_PRIMARY_LINK",
                "B_PRIMARY_LINK",
                "SPLIT_A_TO_B",
                "SPLIT_B_TO_A",
            ]
        )

    lcode = frame[left_code].to_numpy(dtype=np.int64, copy=False)
    rcode = frame[right_code].to_numpy(dtype=np.int64, copy=False)
    swap = lcode > rcode

    a_code = np.where(swap, rcode, lcode)
    b_code = np.where(swap, lcode, rcode)
    left_refs = frame[left_ref].astype(str).to_numpy()
    right_refs = frame[right_ref].astype(str).to_numpy()
    a_ref = np.where(swap, right_refs, left_refs)
    b_ref = np.where(swap, left_refs, right_refs)

    if evidence == "mate":
        a_primary = np.ones(len(frame), dtype=bool)
        b_primary = np.ones(len(frame), dtype=bool)
        split_a_to_b = np.zeros(len(frame), dtype=bool)
        split_b_to_a = np.zeros(len(frame), dtype=bool)
        has_mate = np.ones(len(frame), dtype=bool)
        has_split = np.zeros(len(frame), dtype=bool)
    elif evidence == "split":
        # The left endpoint is the primary alignment and the right endpoint is
        # the supplementary alignment in the caller below.
        if not left_is_primary or right_is_primary:
            raise ValueError("split event must be primary -> supplementary")
        a_primary = ~swap
        b_primary = swap
        split_a_to_b = ~swap
        split_b_to_a = swap
        has_mate = np.zeros(len(frame), dtype=bool)
        has_split = np.ones(len(frame), dtype=bool)
    else:
        raise ValueError(f"Unknown evidence type: {evidence}")

    return pd.DataFrame(
        {
            "QNAME": frame["QNAME"].to_numpy(),
            "A_CODE": a_code,
            "B_CODE": b_code,
            "A_REF": a_ref,
            "B_REF": b_ref,
            "HAS_MATE": has_mate,
            "HAS_SPLIT": has_split,
            "A_PRIMARY_LINK": a_primary,
            "B_PRIMARY_LINK": b_primary,
            "SPLIT_A_TO_B": split_a_to_b,
            "SPLIT_B_TO_A": split_b_to_a,
        }
    )


def _build_relation_events(
    df: pd.DataFrame,
    species_to_code: dict[Hashable, int],
    include_mate_links: bool,
    include_supplementary_links: bool,
) -> pd.DataFrame:
    """
    Build direct inter-species relationship events.

    Evidence comes from only two sources:
      1. discordant primary mates from the same template; and
      2. primary -> supplementary split alignments for the same read segment.

    Secondary alignments never enter this function.
    """
    primary = _best_primary_per_segment(df)
    primary = primary.loc[primary["SPECIES_TAXID"].isin(species_to_code)].copy()
    primary["SCODE"] = primary["SPECIES_TAXID"].map(species_to_code).astype(np.int64)

    events = []

    if include_mate_links:
        p1 = primary.loc[primary["MATE"] == 1, ["QNAME", "SCODE", "REF"]].rename(
            columns={"SCODE": "LEFT_CODE", "REF": "LEFT_REF"}
        )
        p2 = primary.loc[primary["MATE"] == 2, ["QNAME", "SCODE", "REF"]].rename(
            columns={"SCODE": "RIGHT_CODE", "REF": "RIGHT_REF"}
        )
        if not p1.empty and not p2.empty:
            mate = p1.merge(p2, on="QNAME", how="inner", sort=False)
            mate = mate.loc[mate["LEFT_CODE"] != mate["RIGHT_CODE"]]
            if not mate.empty:
                events.append(
                    _canonicalize_events(
                        mate,
                        "LEFT_CODE",
                        "RIGHT_CODE",
                        "LEFT_REF",
                        "RIGHT_REF",
                        evidence="mate",
                        left_is_primary=True,
                        right_is_primary=True,
                    )
                )

    if include_supplementary_links:
        supplementary = df.loc[df["IS_SUPPLEMENTARY"]].copy()
        supplementary = supplementary.loc[
            supplementary["SPECIES_TAXID"].isin(species_to_code)
        ].copy()
        if not supplementary.empty and not primary.empty:
            supplementary["SCODE"] = (
                supplementary["SPECIES_TAXID"].map(species_to_code).astype(np.int64)
            )

            p = primary[["QNAME", "MATE", "SCODE", "REF"]].rename(
                columns={"SCODE": "LEFT_CODE", "REF": "LEFT_REF"}
            )
            s = supplementary[["QNAME", "MATE", "SCODE", "REF"]].rename(
                columns={"SCODE": "RIGHT_CODE", "REF": "RIGHT_REF"}
            )
            split = p.merge(s, on=["QNAME", "MATE"], how="inner", sort=False)
            split = split.loc[split["LEFT_CODE"] != split["RIGHT_CODE"]]
            if not split.empty:
                events.append(
                    _canonicalize_events(
                        split,
                        "LEFT_CODE",
                        "RIGHT_CODE",
                        "LEFT_REF",
                        "RIGHT_REF",
                        evidence="split",
                        left_is_primary=True,
                        right_is_primary=False,
                    )
                )

    if not events:
        return _canonicalize_events(
            pd.DataFrame(),
            "LEFT_CODE",
            "RIGHT_CODE",
            "LEFT_REF",
            "RIGHT_REF",
            evidence="mate",
            left_is_primary=True,
            right_is_primary=True,
        )

    return pd.concat(events, ignore_index=True)


def _locus_metrics(
    events: pd.DataFrame,
    shared_templates: pd.Series,
    side: str,
    ref_lengths: dict[str, float],
) -> pd.DataFrame:
    """Calculate locus diversity/concentration for one endpoint of each edge."""
    ref_col = f"{side}_REF"
    prefix = side.lower()
    pair_cols = ["A_CODE", "B_CODE"]

    unique_refs = events[pair_cols + [ref_col]].drop_duplicates()
    distinct = unique_refs.groupby(pair_cols, sort=False).size().rename(
        f"{prefix}_distinct_loci"
    )

    unique_refs = unique_refs.copy()
    unique_refs["SIG_BP"] = unique_refs[ref_col].map(ref_lengths)
    sig_bp = unique_refs.groupby(pair_cols, sort=False)["SIG_BP"].sum(
        min_count=1
    ).rename(f"{prefix}_signature_bp")

    qref = events[["QNAME"] + pair_cols + [ref_col]].drop_duplicates()
    per_locus = qref.groupby(pair_cols + [ref_col], sort=False).size()
    top = per_locus.groupby(level=[0, 1], sort=False).max().rename("TOP_COUNT")
    top_fraction = (top / shared_templates).rename(f"{prefix}_top_locus_fraction")

    return pd.concat([distinct, sig_bp, top_fraction], axis=1)


def get_species_hit_groups(
    df: pd.DataFrame,
    include_singletons: bool = True,
    same_genus_only: bool = True,
    genus_col: str = "GENUS_TAXID",
    min_shared_templates: int = 3,
    min_dependency_fraction: float = 0.005,
    fraction_confidence: Optional[float] = 0.95,
    min_relative_strength: float = 0.0,
    include_mate_links: bool = True,
    include_supplementary_links: bool = True,
    return_edges: bool = False,
    return_nodes: bool = False,
    *,
    min_reciprocal_reads: Optional[int] = None,
    min_reciprocal_fraction: Optional[float] = None,
):
    """
    Group primary-supported species using template-level mapping dependency.

    This replaces the old first-alignment/reciprocal-direction model.  No
    secondary alignments are used.  A relationship is created only when:

      * the two primary mates of a paired-end template map to different
        species; or
      * a primary alignment and a supplementary alignment from the same read
        segment map to different species.

    For species pair A--B, the method measures the fraction of A's primary
    templates that participate in an A/B relationship and the corresponding
    fraction for B.  The edge strength is the larger of the two fractions (or
    Wilson lower bounds).  This is intentional: a low-abundance cross-mapping
    shadow can be almost completely dependent on an abundant source species,
    while only a small fraction of the source species' templates point back to
    the shadow.

    Parameters
    ----------
    df
        Alignment DataFrame. Required columns are QNAME, FLAG, REF, and
        SPECIES_TAXID. GENUS_TAXID is required when ``same_genus_only=True``.
        AS is optional for this function (it is used upstream when loading SAM).
    include_singletons
        Include primary-supported species without retained edges as one-member
        groups.
    same_genus_only
        Restrict graph edges to species assigned to the same genus.
    genus_col
        Column containing genus taxids/names.
    min_shared_templates
        Minimum number of independent QNAME templates supporting an A--B
        relationship.
    min_dependency_fraction
        Minimum dependency strength.  The dependency strength is max(f_A, f_B),
        where f_A is the fraction of A primary templates linked to B and f_B is
        the fraction of B primary templates linked to A.  Wilson lower bounds
        are used when ``fraction_confidence`` is not None.
    fraction_confidence
        Confidence level for Wilson lower bounds. Use None for raw fractions.
    min_relative_strength
        Optional local pruning in [0, 1].  An edge must be at least this
        fraction of the strongest edge for its more-dependent endpoint.
    include_mate_links
        Use discordant primary paired-end mates.
    include_supplementary_links
        Use primary/supplementary split relationships.
    return_edges
        Also return an edge diagnostics table.
    return_nodes
        Also return a node diagnostics table.
    min_reciprocal_reads, min_reciprocal_fraction
        Backward-compatible aliases for the old argument names. They now map
        to ``min_shared_templates`` and ``min_dependency_fraction``.

    Returns
    -------
    dict or tuple
        Groups keyed as group1, group2, ... . Depending on ``return_edges`` and
        ``return_nodes``, edge and node diagnostic tables are also returned.
    """
    if min_reciprocal_reads is not None:
        min_shared_templates = min_reciprocal_reads
    if min_reciprocal_fraction is not None:
        min_dependency_fraction = min_reciprocal_fraction

    if min_shared_templates < 1:
        raise ValueError("min_shared_templates must be >= 1")
    if not 0.0 <= min_dependency_fraction <= 1.0:
        raise ValueError("min_dependency_fraction must be in [0, 1]")
    if not 0.0 <= min_relative_strength <= 1.0:
        raise ValueError("min_relative_strength must be in [0, 1]")
    if not include_mate_links and not include_supplementary_links:
        raise ValueError("At least one relationship evidence type must be enabled")

    required = {"QNAME", "FLAG", "REF", "SPECIES_TAXID"}
    if same_genus_only:
        required.add(genus_col)
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"Missing required column(s): {sorted(missing)}")

    work = _prepare_alignment_df(df)
    primary = _best_primary_per_segment(work)

    # Graph nodes are species that have primary support somewhere in the sample.
    species = pd.Index(primary["SPECIES_TAXID"].dropna().unique())
    n_species = len(species)
    species_to_code = {sp: i for i, sp in enumerate(species)}

    edge_columns = [
        "species_a",
        "species_b",
        "shared_templates",
        "mate_templates",
        "split_templates",
        "mixed_templates",
        "a_primary_templates",
        "b_primary_templates",
        "a_linked_primary_templates",
        "b_linked_primary_templates",
        "a_dependency_fraction",
        "b_dependency_fraction",
        "a_dependency_lower",
        "b_dependency_lower",
        "edge_strength",
        "mutual_strength",
        "more_dependent_species",
        "candidate_anchor_species",
        "split_a_primary_to_b_supp",
        "split_b_primary_to_a_supp",
        "split_direction_balance",
        "a_distinct_loci",
        "b_distinct_loci",
        "a_signature_bp",
        "b_signature_bp",
        "a_top_locus_fraction",
        "b_top_locus_fraction",
        "passes_template_threshold",
        "passes_dependency_threshold",
        "passes_relative_threshold",
        "kept",
    ]
    node_columns = [
        "species",
        "primary_templates",
        "linked_primary_templates",
        "independent_primary_templates",
        "linked_fraction",
        "supplementary_target_templates",
        "kept_degree",
    ]

    if n_species == 0:
        groups = {}
        edge_table = pd.DataFrame(columns=edge_columns)
        node_table = pd.DataFrame(columns=node_columns)
        if return_edges and return_nodes:
            return groups, edge_table, node_table
        if return_edges:
            return groups, edge_table
        if return_nodes:
            return groups, node_table
        return groups

    # Species -> genus mapping is determined from all retained non-secondary
    # records but must be unique for each primary-supported species.
    species_genus = None
    if same_genus_only:
        genus_map = work.loc[
            work["SPECIES_TAXID"].isin(species),
            ["SPECIES_TAXID", genus_col],
        ].dropna()
        genus_counts = genus_map.groupby("SPECIES_TAXID", sort=False)[genus_col].nunique()
        inconsistent = genus_counts[genus_counts > 1]
        if not inconsistent.empty:
            raise ValueError(
                "Each species must map to only one genus; inconsistent species: "
                f"{inconsistent.index.tolist()}"
            )
        genus_first = genus_map.drop_duplicates("SPECIES_TAXID", keep="first").set_index(
            "SPECIES_TAXID"
        )[genus_col]
        species_genus = np.array([genus_first.get(sp, None) for sp in species], dtype=object)

    primary = primary.loc[primary["SPECIES_TAXID"].isin(species)].copy()
    primary["SCODE"] = primary["SPECIES_TAXID"].map(species_to_code).astype(np.int64)
    primary_template_counts = (
        primary[["QNAME", "SCODE"]]
        .drop_duplicates()
        .groupby("SCODE", sort=False)
        .size()
        .reindex(range(n_species), fill_value=0)
        .astype(np.int64)
    )

    events = _build_relation_events(
        work,
        species_to_code,
        include_mate_links=include_mate_links,
        include_supplementary_links=include_supplementary_links,
    )

    if same_genus_only and not events.empty:
        a_genus = species_genus[events["A_CODE"].to_numpy(dtype=np.int64)]
        b_genus = species_genus[events["B_CODE"].to_numpy(dtype=np.int64)]
        genus_ok = pd.notna(a_genus) & pd.notna(b_genus) & (a_genus == b_genus)
        events = events.loc[genus_ok].copy()

    if events.empty:
        kept_a = np.empty(0, dtype=np.int64)
        kept_b = np.empty(0, dtype=np.int64)
        edge_table = pd.DataFrame(columns=edge_columns)
        linked_counts = pd.Series(dtype=np.int64)
        supplementary_target_counts = pd.Series(dtype=np.int64)
    else:
        # Collapse multiple split segments or mixed evidence from the same
        # template into a single template/species-pair relationship.
        template_pair = (
            events.groupby(["QNAME", "A_CODE", "B_CODE"], sort=False)
            .agg(
                HAS_MATE=("HAS_MATE", "max"),
                HAS_SPLIT=("HAS_SPLIT", "max"),
                A_PRIMARY_LINK=("A_PRIMARY_LINK", "max"),
                B_PRIMARY_LINK=("B_PRIMARY_LINK", "max"),
                SPLIT_A_TO_B=("SPLIT_A_TO_B", "max"),
                SPLIT_B_TO_A=("SPLIT_B_TO_A", "max"),
            )
            .reset_index()
        )
        template_pair["MIXED"] = template_pair["HAS_MATE"] & template_pair["HAS_SPLIT"]

        grouped = template_pair.groupby(["A_CODE", "B_CODE"], sort=False)
        edge = grouped.agg(
            shared_templates=("QNAME", "size"),
            mate_templates=("HAS_MATE", "sum"),
            split_templates=("HAS_SPLIT", "sum"),
            mixed_templates=("MIXED", "sum"),
            a_linked_primary_templates=("A_PRIMARY_LINK", "sum"),
            b_linked_primary_templates=("B_PRIMARY_LINK", "sum"),
            split_a_primary_to_b_supp=("SPLIT_A_TO_B", "sum"),
            split_b_primary_to_a_supp=("SPLIT_B_TO_A", "sum"),
        )

        a_codes = edge.index.get_level_values(0).to_numpy(dtype=np.int64)
        b_codes = edge.index.get_level_values(1).to_numpy(dtype=np.int64)
        n_a = primary_template_counts.to_numpy()[a_codes]
        n_b = primary_template_counts.to_numpy()[b_codes]
        a_linked = edge["a_linked_primary_templates"].to_numpy(dtype=np.int64)
        b_linked = edge["b_linked_primary_templates"].to_numpy(dtype=np.int64)

        frac_a = a_linked / n_a
        frac_b = b_linked / n_b
        if fraction_confidence is None:
            lower_a = frac_a.copy()
            lower_b = frac_b.copy()
        else:
            lower_a = _wilson_lower_bound(a_linked, n_a, fraction_confidence)
            lower_b = _wilson_lower_bound(b_linked, n_b, fraction_confidence)

        strength = np.maximum(lower_a, lower_b)
        mutual_strength = np.minimum(lower_a, lower_b)
        a_more_dependent = lower_a >= lower_b
        dependent_code = np.where(a_more_dependent, a_codes, b_codes)
        anchor_code = np.where(a_more_dependent, b_codes, a_codes)

        split_ab = edge["split_a_primary_to_b_supp"].to_numpy(dtype=np.float64)
        split_ba = edge["split_b_primary_to_a_supp"].to_numpy(dtype=np.float64)
        split_max = np.maximum(split_ab, split_ba)
        split_min = np.minimum(split_ab, split_ba)
        split_balance = np.divide(
            split_min,
            split_max,
            out=np.full_like(split_min, np.nan, dtype=np.float64),
            where=split_max > 0,
        )

        shared = edge["shared_templates"].to_numpy(dtype=np.int64)
        passes_templates = shared >= min_shared_templates
        passes_dependency = strength >= min_dependency_fraction
        base_keep = passes_templates & passes_dependency

        passes_relative = np.ones(len(edge), dtype=bool)
        if min_relative_strength > 0 and base_keep.any():
            strongest_for_dependent = np.zeros(n_species, dtype=np.float64)
            np.maximum.at(
                strongest_for_dependent,
                dependent_code[base_keep],
                strength[base_keep],
            )
            passes_relative = strength >= (
                min_relative_strength * strongest_for_dependent[dependent_code]
            )

        kept = base_keep & passes_relative

        # Reference-position diagnostics.  Exact RNAME is treated as a locus;
        # when start/end are encoded in RNAME, total implicated signature bp is
        # also reported.  These metrics are intentionally diagnostic rather
        # than hard filters because localized cross-mapping is often precisely
        # the artifact we want the graph to reveal.
        ref_lengths = {
            ref: _signature_length_from_ref(ref)
            for ref in pd.unique(pd.concat([events["A_REF"], events["B_REF"]]))
        }
        shared_series = edge["shared_templates"]
        locus_a = _locus_metrics(events, shared_series, "A", ref_lengths)
        locus_b = _locus_metrics(events, shared_series, "B", ref_lengths)
        edge = edge.join(locus_a).join(locus_b)

        edge_table = edge.reset_index()
        edge_table["species_a"] = species.take(edge_table["A_CODE"].to_numpy()).to_numpy()
        edge_table["species_b"] = species.take(edge_table["B_CODE"].to_numpy()).to_numpy()
        edge_table["a_primary_templates"] = n_a
        edge_table["b_primary_templates"] = n_b
        edge_table["a_dependency_fraction"] = frac_a
        edge_table["b_dependency_fraction"] = frac_b
        edge_table["a_dependency_lower"] = lower_a
        edge_table["b_dependency_lower"] = lower_b
        edge_table["edge_strength"] = strength
        edge_table["mutual_strength"] = mutual_strength
        edge_table["more_dependent_species"] = species.take(dependent_code).to_numpy()
        edge_table["candidate_anchor_species"] = species.take(anchor_code).to_numpy()
        edge_table["split_direction_balance"] = split_balance
        edge_table["passes_template_threshold"] = passes_templates
        edge_table["passes_dependency_threshold"] = passes_dependency
        edge_table["passes_relative_threshold"] = passes_relative
        edge_table["kept"] = kept

        edge_table = edge_table[edge_columns].sort_values(
            ["kept", "edge_strength", "shared_templates", "mutual_strength"],
            ascending=[False, False, False, False],
            ignore_index=True,
        )

        kept_a = a_codes[kept]
        kept_b = b_codes[kept]

        # Node-level dependency: only QNAMEs in which the species itself has a
        # primary alignment count as linked primary templates.
        a_link_rows = template_pair.loc[
            template_pair["A_PRIMARY_LINK"], ["QNAME", "A_CODE"]
        ].rename(columns={"A_CODE": "SCODE"})
        b_link_rows = template_pair.loc[
            template_pair["B_PRIMARY_LINK"], ["QNAME", "B_CODE"]
        ].rename(columns={"B_CODE": "SCODE"})
        linked_counts = (
            pd.concat([a_link_rows, b_link_rows], ignore_index=True)
            .drop_duplicates(["QNAME", "SCODE"])
            .groupby("SCODE", sort=False)
            .size()
        )

        # Count templates where a species appears as the supplementary target
        # of a split relation, regardless of whether it was primary elsewhere.
        split_events = events.loc[events["HAS_SPLIT"]].copy()
        supp_rows = []
        a_to_b = split_events.loc[
            split_events["SPLIT_A_TO_B"], ["QNAME", "B_CODE"]
        ].rename(columns={"B_CODE": "SCODE"})
        if not a_to_b.empty:
            supp_rows.append(a_to_b)
        b_to_a = split_events.loc[
            split_events["SPLIT_B_TO_A"], ["QNAME", "A_CODE"]
        ].rename(columns={"A_CODE": "SCODE"})
        if not b_to_a.empty:
            supp_rows.append(b_to_a)
        if supp_rows:
            supplementary_target_counts = (
                pd.concat(supp_rows, ignore_index=True)
                .drop_duplicates(["QNAME", "SCODE"])
                .groupby("SCODE", sort=False)
                .size()
            )
        else:
            supplementary_target_counts = pd.Series(dtype=np.int64)

    # Connected components are built from retained dependency edges.  This
    # groups likely source/shadow species without claiming which one is true.
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

    # Node diagnostics.
    primary_counts_arr = primary_template_counts.to_numpy(dtype=np.int64)
    linked_arr = np.zeros(n_species, dtype=np.int64)
    if not linked_counts.empty:
        idx = linked_counts.index.to_numpy(dtype=np.int64)
        linked_arr[idx] = linked_counts.to_numpy(dtype=np.int64)
    supplementary_arr = np.zeros(n_species, dtype=np.int64)
    if not supplementary_target_counts.empty:
        idx = supplementary_target_counts.index.to_numpy(dtype=np.int64)
        supplementary_arr[idx] = supplementary_target_counts.to_numpy(dtype=np.int64)

    degree = np.zeros(n_species, dtype=np.int64)
    if len(kept_a):
        np.add.at(degree, kept_a, 1)
        np.add.at(degree, kept_b, 1)

    node_table = pd.DataFrame(
        {
            "species": species,
            "primary_templates": primary_counts_arr,
            "linked_primary_templates": linked_arr,
            "independent_primary_templates": primary_counts_arr - linked_arr,
            "linked_fraction": np.divide(
                linked_arr,
                primary_counts_arr,
                out=np.zeros(n_species, dtype=np.float64),
                where=primary_counts_arr > 0,
            ),
            "supplementary_target_templates": supplementary_arr,
            "kept_degree": degree,
        }
    ).sort_values(
        ["linked_fraction", "primary_templates"],
        ascending=[False, False],
        ignore_index=True,
    )

    if return_edges and return_nodes:
        return groups, edge_table, node_table
    if return_edges:
        return groups, edge_table
    if return_nodes:
        return groups, node_table
    return groups


def reciprocal_relationships_from_sam(
    samfile: Path,
    min_alen: int,
    *,
    include_singletons: bool = True,
    same_genus_only: bool = True,
    min_shared_templates: int = 3,
    min_dependency_fraction: float = 0.005,
    fraction_confidence: Optional[float] = 0.95,
    min_relative_strength: float = 0.0,
    include_mate_links: bool = True,
    include_supplementary_links: bool = True,
    return_edges: bool = False,
    return_nodes: bool = False,
):
    """
    Extract primary-mate/supplementary species relationships from a SAM file.

    Secondary alignments are explicitly discarded.  Graph nodes are species
    with primary support.  Edges summarize discordant primary mates and/or
    primary-to-supplementary split mappings.
    """
    logging.info("Loading primary and supplementary SAM alignments...")
    df = _load_sam_without_secondary(Path(samfile), min_alen=min_alen)
    logging.debug("Loaded %d retained SAM alignments", len(df))

    if df.empty:
        empty = get_species_hit_groups(
            pd.DataFrame(columns=["QNAME", "FLAG", "REF", "SPECIES_TAXID", "GENUS_TAXID"]),
            include_singletons=include_singletons,
            same_genus_only=same_genus_only,
            min_shared_templates=min_shared_templates,
            min_dependency_fraction=min_dependency_fraction,
            fraction_confidence=fraction_confidence,
            min_relative_strength=min_relative_strength,
            include_mate_links=include_mate_links,
            include_supplementary_links=include_supplementary_links,
            return_edges=return_edges,
            return_nodes=return_nodes,
        )
        return empty

    logging.info("Mapping reference taxids to species/genus taxids...")
    ref_meta = {ref: _parse_ref_metadata(ref) for ref in df["REF"].unique()}
    ref_taxid = {ref: meta[3] for ref, meta in ref_meta.items()}
    df["TAXID"] = df["REF"].map(ref_taxid)

    taxids = pd.Index(df["TAXID"].dropna().unique())
    taxid_to_species = {
        taxid: t.taxid2taxidOnRank(taxid, target_rank="species")
        for taxid in taxids
    }
    taxid_to_genus = {
        taxid: t.taxid2taxidOnRank(taxid, target_rank="genus")
        for taxid in taxids
    }
    df["SPECIES_TAXID"] = df["TAXID"].map(taxid_to_species)
    df["GENUS_TAXID"] = df["TAXID"].map(taxid_to_genus)

    # Retain only QNAMEs that contain a real cross-species relationship under
    # the primary-mate/supplementary model.  get_species_hit_groups performs
    # the definitive event construction; pre-filtering here is intentionally
    # omitted so primary-template denominators remain correct.
    logging.info(
        "Grouping species from discordant primary mates and supplementary split mappings..."
    )
    return get_species_hit_groups(
        df,
        include_singletons=include_singletons,
        same_genus_only=same_genus_only,
        min_shared_templates=min_shared_templates,
        min_dependency_fraction=min_dependency_fraction,
        fraction_confidence=fraction_confidence,
        min_relative_strength=min_relative_strength,
        include_mate_links=include_mate_links,
        include_supplementary_links=include_supplementary_links,
        return_edges=return_edges,
        return_nodes=return_nodes,
    )
