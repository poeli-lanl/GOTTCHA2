from __future__ import annotations

import importlib
import math
import pathlib
import sys
import types
from collections import Counter

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
UTILS = ROOT / "gottcha" / "utils"
if str(UTILS) not in sys.path:
    sys.path.insert(0, str(UTILS))

from species_linkage import (  # noqa: E402
    EvidenceConfig,
    EvidenceResult,
    GroupConfig,
    LinkageObservation,
    LinkageObservationStore,
    PairEvidence,
    SpeciesEvidence,
    build_linkage_evidence,
    cluster_read_observations,
    finalize_species_groups,
    make_read_key,
    query_interval_from_cigar,
    write_species_evidence,
    write_species_groups,
    write_species_links,
)


def obs(
    read: str,
    qstart: int,
    qend: int,
    species: int,
    reference: int,
    score: float,
    *,
    identity: float = 0.99,
    alignment_class: int = 0,
    has_as: int = 1,
) -> LinkageObservation:
    return LinkageObservation(
        read,
        qstart,
        qend,
        species,
        reference,
        score,
        identity,
        qend - qstart,
        alignment_class,
        has_as,
    )


def permissive_evidence_config(**overrides) -> EvidenceConfig:
    values = dict(
        min_query_overlap=0.80,
        max_species_per_segment=5,
        score_tau=0.05,
        length_saturation=100,
        min_alignment_weight=0.0,
        independent_score_margin=0.05,
    )
    values.update(overrides)
    return EvidenceConfig(**values)


def permissive_group_config(**overrides) -> GroupConfig:
    values = dict(
        min_link_reads=1,
        min_link_segments=1,
        min_link_loci=1,
        min_link_weight=0.0,
        same_genus_only=False,
        min_shadow_containment=0.50,
        min_mutual_ambiguity=0.50,
        max_shadow_independent_fraction=0.20,
        min_independent_loci_to_retain=2,
        min_anchor_independent_loci_for_shadow=2,
    )
    values.update(overrides)
    return GroupConfig(**values)


def test_query_interval_reconstructs_forward_and_reverse_coordinates() -> None:
    # 10H + 5S + 100M + 7I + 3D + 20S + 8H
    cigar = [(5, 10), (4, 5), (0, 100), (1, 7), (2, 3), (4, 20), (5, 8)]
    # Original span is 150 query bases, including hard clips.
    assert query_interval_from_cigar(cigar, False) == (15, 122)
    assert query_interval_from_cigar(cigar, True) == (28, 135)


def test_nonoverlapping_ont_segments_do_not_create_species_link() -> None:
    read_groups = [
        (
            "read-1",
            [
                obs("read-1", 0, 500, 0, 10, 1000),
                obs("read-1", 4000, 4500, 1, 20, 990),
            ],
        )
    ]
    result = build_linkage_evidence(
        read_groups,
        n_species=2,
        config=permissive_evidence_config(),
    )
    assert result.pair_evidence == {}
    assert result.species_evidence[0].independent_weight == 1.0
    assert result.species_evidence[1].independent_weight == 1.0
    assert result.segment_degree_histogram == Counter({1: 2})


def test_overlapping_primary_and_secondary_create_weighted_link() -> None:
    read_groups = [
        (
            "read-1",
            [
                obs("read-1", 0, 100, 0, 10, 200, alignment_class=0),
                obs("read-1", 0, 100, 1, 20, 198, alignment_class=1),
            ],
        )
    ]
    result = build_linkage_evidence(
        read_groups,
        n_species=2,
        config=permissive_evidence_config(),
    )
    pair = result.pair_evidence[(0, 1)]
    expected = math.exp(-((200 - 198) / 100) / 0.05)
    assert pair.shared_weight == pytest.approx(expected)
    assert pair.shared_reads == 1
    assert pair.shared_segments == 1
    assert pair.locus_pairs == {(10, 20)}


def test_one_query_segment_has_at_most_one_unit_cross_species_mass() -> None:
    read_groups = [
        (
            "read-1",
            [
                obs("read-1", 0, 100, 0, 10, 200),
                obs("read-1", 0, 100, 1, 20, 200, alignment_class=1),
                obs("read-1", 0, 100, 2, 30, 200, alignment_class=1),
            ],
        )
    ]
    result = build_linkage_evidence(
        read_groups,
        n_species=3,
        config=permissive_evidence_config(),
    )
    assert result.pair_evidence[(0, 1)].shared_weight == pytest.approx(0.5)
    assert result.pair_evidence[(0, 2)].shared_weight == pytest.approx(0.5)
    assert sum(pair.shared_weight for pair in result.pair_evidence.values()) == pytest.approx(1.0)


def test_pair_shared_reads_is_deduplicated_across_segments_on_same_read() -> None:
    read_groups = [
        (
            "read-1",
            [
                obs("read-1", 0, 100, 0, 10, 200),
                obs("read-1", 0, 100, 1, 20, 200, alignment_class=1),
                obs("read-1", 500, 600, 0, 11, 200),
                obs("read-1", 500, 600, 1, 21, 200, alignment_class=1),
            ],
        )
    ]
    result = build_linkage_evidence(
        read_groups,
        n_species=2,
        config=permissive_evidence_config(),
    )
    pair = result.pair_evidence[(0, 1)]
    assert pair.shared_reads == 1
    assert pair.shared_segments == 2
    assert pair.locus_pairs == {(10, 20), (11, 21)}


def test_cluster_is_anchor_centered_not_overlap_transitive() -> None:
    observations = [
        obs("r", 0, 100, 0, 10, 300),
        obs("r", 20, 120, 1, 20, 200),  # 80% overlap with first
        obs("r", 40, 140, 2, 30, 100),  # 80% with second, 60% with first
    ]
    clusters = cluster_read_observations(observations, min_query_overlap=0.80)
    assert len(clusters) == 2
    assert {item.species_idx for item in clusters[0]} == {0, 1}
    assert {item.species_idx for item in clusters[1]} == {2}


def test_anchor_centered_species_groups_do_not_follow_transitive_chain() -> None:
    species = [
        SpeciesEvidence(total_weight=10.0, independent_weight=5.0, read_count=10),
        SpeciesEvidence(total_weight=6.0, independent_weight=2.0, read_count=6),
        SpeciesEvidence(total_weight=4.0, independent_weight=1.0, read_count=4),
    ]
    pairs = {
        (0, 1): PairEvidence(
            shared_weight=4.0,
            shared_segments=4,
            shared_reads=4,
            locus_pairs={(1, 2)},
        ),
        (1, 2): PairEvidence(
            shared_weight=3.0,
            shared_segments=3,
            shared_reads=3,
            locus_pairs={(2, 3)},
        ),
    }
    evidence = EvidenceResult(species, pairs, Counter(), Counter(), Counter())
    groups = finalize_species_groups(
        ["A", "B", "C"],
        ["G", "G", "G"],
        evidence,
        permissive_group_config(
            min_shadow_containment=0.60,
            min_mutual_ambiguity=0.40,
        ),
    )
    assert groups.groups == {1: ["A", "B"]}
    assert "C" not in groups.species_to_group


def test_same_genus_filter_blocks_cross_genus_group() -> None:
    species = [
        SpeciesEvidence(total_weight=2.0, read_count=2),
        SpeciesEvidence(total_weight=2.0, read_count=2),
    ]
    pairs = {
        (0, 1): PairEvidence(
            shared_weight=2.0,
            shared_segments=2,
            shared_reads=2,
            locus_pairs={(1, 2)},
        )
    }
    evidence = EvidenceResult(species, pairs, Counter(), Counter(), Counter())
    groups = finalize_species_groups(
        ["A", "B"],
        ["G1", "G2"],
        evidence,
        permissive_group_config(same_genus_only=True),
    )
    assert groups.groups == {}
    assert groups.edge_rows[0]["PASS_GENUS"] == 0


def test_shadow_requires_independent_anchor_loci() -> None:
    anchor = SpeciesEvidence(
        total_weight=10.0,
        independent_weight=6.0,
        read_count=10,
        independent_signature_ids={1, 2},
    )
    candidate = SpeciesEvidence(
        total_weight=3.0,
        independent_weight=0.0,
        read_count=3,
        independent_signature_ids=set(),
    )
    pair = PairEvidence(
        shared_weight=3.0,
        shared_segments=3,
        shared_reads=3,
        locus_pairs={(1, 3), (2, 4)},
    )
    evidence = EvidenceResult(
        [anchor, candidate], {(0, 1): pair}, Counter(), Counter(), Counter()
    )
    config = permissive_group_config(
        min_shadow_containment=0.80,
        min_mutual_ambiguity=0.80,
        min_anchor_independent_loci_for_shadow=2,
    )
    result = finalize_species_groups(["A", "B"], ["G", "G"], evidence, config)
    assert result.shadows == {"B": "A"}

    anchor.independent_signature_ids = {1}
    result = finalize_species_groups(["A", "B"], ["G", "G"], evidence, config)
    assert result.shadows == {}
    assert result.group_meta[1]["ambiguous_taxids"] == ["B"]



def test_multiple_signature_pairs_require_support_on_both_species() -> None:
    species = [
        SpeciesEvidence(total_weight=4.0, read_count=4),
        SpeciesEvidence(total_weight=4.0, read_count=4),
    ]
    # One signature in species A paired with three signatures in species B.
    # This is still only one independently supported locus on the A side.
    pair = PairEvidence(
        shared_weight=3.0,
        shared_segments=3,
        shared_reads=3,
        locus_pairs={(10, 20), (10, 21), (10, 22)},
    )
    evidence = EvidenceResult(
        species, {(0, 1): pair}, Counter(), Counter(), Counter()
    )
    groups = finalize_species_groups(
        ["A", "B"],
        ["G", "G"],
        evidence,
        permissive_group_config(min_link_loci=2),
    )
    assert groups.groups == {}
    row = groups.edge_rows[0]
    assert row["SHARED_SIGNATURE_PAIRS"] == 3
    assert row["SHARED_INDEPENDENT_LOCI"] == 1
    assert row["PASS_MIN_SIGNATURES"] == 0


def test_default_thresholds_reject_single_read_bridge_chain() -> None:
    species = [SpeciesEvidence(total_weight=1.0, read_count=1) for _ in range(5)]
    pairs = {}
    for idx in range(4):
        pairs[(idx, idx + 1)] = PairEvidence(
            shared_weight=1.0,
            shared_segments=1,
            shared_reads=1,
            locus_pairs={(idx, idx + 1)},
        )
    evidence = EvidenceResult(species, pairs, Counter(), Counter(), Counter())
    groups = finalize_species_groups(
        [str(idx) for idx in range(5)],
        ["G"] * 5,
        evidence,
        GroupConfig(),
    )
    assert groups.groups == {}
    assert all(row["QUALIFIED"] == 0 for row in groups.edge_rows)

def test_storage_round_trip_memory_and_sqlite(tmp_path: pathlib.Path) -> None:
    observations = [
        obs("read-b", 0, 100, 0, 10, 100),
        obs("read-a", 0, 100, 1, 20, 90),
        obs("read-a", 0, 100, 0, 10, 100),
    ]
    for mode in ("memory", "sqlite"):
        store = LinkageObservationStore(mode=mode, temp_dir=str(tmp_path))
        path = store.path
        store.add_many(observations)
        groups = {read: rows for read, rows in store.iter_reads()}
        assert set(groups) == {"read-a", "read-b"}
        assert len(groups["read-a"]) == 2
        store.close()
        if mode == "sqlite":
            assert path is not None
            assert not pathlib.Path(path).exists()


def test_make_read_key_separates_read_groups_and_mates() -> None:
    assert make_read_key("q", "rg1", 1) != make_read_key("q", "rg1", 2)
    assert make_read_key("q", "rg1", 1) != make_read_key("q", "rg2", 1)


def test_diagnostic_writers_emit_consistent_tables(tmp_path: pathlib.Path) -> None:
    reads = [
        (
            "r1",
            [
                obs("r1", 0, 100, 0, 10, 200),
                obs("r1", 0, 100, 1, 20, 200, alignment_class=1),
            ],
        )
    ]
    evidence = build_linkage_evidence(
        reads, 2, permissive_evidence_config()
    )
    groups = finalize_species_groups(
        ["A", "B"],
        ["G", "G"],
        evidence,
        permissive_group_config(),
    )
    groups_path = tmp_path / "groups.tsv"
    links_path = tmp_path / "links.tsv"
    evidence_path = tmp_path / "evidence.tsv"
    write_species_groups(groups, str(groups_path))
    write_species_links(groups.edge_rows, str(links_path))
    write_species_evidence(
        ["A", "B"], ["G", "G"], evidence, groups, str(evidence_path)
    )
    assert "ANCHOR_TAXID" in groups_path.read_text().splitlines()[0]
    links_header = links_path.read_text().splitlines()[0].split("\t")
    assert "SHARED_INDEPENDENT_LOCI" in links_header
    assert "SHADOW_DIRECTION" in links_header
    assert "INDEPENDENT_FRACTION" in evidence_path.read_text().splitlines()[0]


class FakeAlignment:
    def __init__(
        self,
        *,
        name: str,
        reference_start: int,
        cigartuples,
        secondary: bool = False,
        supplementary: bool = False,
        reverse: bool = False,
        tags=None,
    ) -> None:
        self.query_name = name
        self.reference_start = reference_start
        self.cigartuples = cigartuples
        self.is_secondary = secondary
        self.is_supplementary = supplementary
        self.is_reverse = reverse
        self.is_unmapped = False
        self.is_duplicate = False
        self.is_qcfail = False
        self.mapping_quality = 60
        self.query_length = sum(
            length for operation, length in cigartuples
            if operation in (0, 1, 4, 7, 8)
        )
        self.query_alignment_length = sum(
            length for operation, length in cigartuples
            if operation in (0, 1, 7, 8)
        )
        self.alen = self.query_alignment_length
        self.is_paired = False
        self.is_read1 = False
        self.is_read2 = False
        self._tags = dict(tags or {})

    def has_tag(self, tag: str) -> bool:
        return tag in self._tags

    def get_tag(self, tag: str):
        return self._tags[tag]


class FakeBam:
    lengths = [1000]

    def __init__(self, alignments) -> None:
        self.alignments = alignments

    def fetch(self, _rname: str, _start: int, _end: int):
        return iter(self.alignments)


def import_process_bam_with_stubs(monkeypatch):
    pysam_stub = types.ModuleType("pysam")
    pysam_stub.AlignedSegment = object
    pysam_stub.AlignmentFile = object
    monkeypatch.setitem(sys.modules, "pysam", pysam_stub)

    taxonomy_stub = types.ModuleType("taxonomy")
    taxonomy_stub.taxid2taxidOnRank = lambda taxid, target_rank: str(taxid)
    monkeypatch.setitem(sys.modules, "taxonomy", taxonomy_stub)

    sys.modules.pop("process_bam", None)
    return importlib.import_module("process_bam")


def test_worker_can_exclude_secondary_from_coverage_but_keep_for_linkage(
    monkeypatch,
) -> None:
    process_bam = import_process_bam_with_stubs(monkeypatch)
    primary = FakeAlignment(
        name="read-1",
        reference_start=10,
        cigartuples=[(7, 100)],
        tags={"AS": 200, "NM": 0},
    )
    secondary = FakeAlignment(
        name="read-2",
        reference_start=200,
        cigartuples=[(7, 100)],
        secondary=True,
        tags={"AS": 198, "NM": 1},
    )
    process_bam._BAM = FakeBam([primary, secondary])
    process_bam._CFG = {
        "coverage": {
            "min_mapq": 0,
            "min_frac": 0.0,
            "min_idt": 0.0,
            "min_alen": 0,
            "include_secondary": False,
            "include_supplementary": False,
            "include_duplicates": False,
            "include_qcfail": False,
        },
        "linkage": {
            "min_mapq": 0,
            "min_frac": 0.0,
            "min_idt": 0.0,
            "min_alen": 0,
            "include_secondary": True,
            "include_supplementary": False,
            "include_duplicates": False,
            "include_qcfail": False,
        },
        "split_read_flag": False,
    }
    metrics, observations = process_bam._process_chunk(
        ("ref|1|x", 0, 500, 0, 0)
    )
    assert metrics[3] == 1  # primary only in NUMREADS
    assert len(observations) == 2
    assert {item.alignment_class for item in observations} == {0, 1}
