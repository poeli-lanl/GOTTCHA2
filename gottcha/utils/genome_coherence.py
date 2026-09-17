#!/usr/bin/env python3
"""Assembly-resolved, molecule-aware evidence diagnostics for GOTTCHA2.

All internal and output intervals are 0-based, half-open. This module never
changes abundance estimates, calls a species present, or reconstructs sequence.
SAM (including headerless SAM and SAM.gz) needs only the Python standard library.
BAM/CRAM additionally needs pysam. Run this file with --help for the CLI.
"""
from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import hashlib
import itertools
import json
import logging
import math
import os
from pathlib import Path
import re
import sqlite3
import sys
import tempfile
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass, field
from typing import Iterable, Optional

VERSION = "0.1.0"
LOG = logging.getLogger("gottcha.coherence")
CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
Interval = tuple[int, int]
MANIFEST_COLUMNS = ["signature_id", "assembly_id", "contig_id", "start0", "end0",
                    "species_taxid", "genus_taxid", "species_name", "contig_length", "topology"]


def text_open(path, mode="rt"):
    """Open UTF-8 text, optionally gzip compressed (including BGZF FASTA)."""
    return (gzip.open if str(path).endswith(".gz") else open)(path, mode, encoding="utf-8")


def read_tsv(path, required=()):
    with text_open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not set(required).issubset(reader.fieldnames):
            raise ValueError(f"{path}: required columns: {', '.join(required)}")
        for number, row in enumerate(reader, 2):
            if None in row or any(v is None for v in row.values()):
                raise ValueError(f"{path}:{number}: inconsistent TSV column count")
            yield row


def write_tsv(path, rows, columns):
    with text_open(path, "wt") as handle:
        writer = csv.DictWriter(handle, columns, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: "" if v is None else v for k, v in row.items()})


def merge_intervals(intervals: Iterable[Interval]) -> list[Interval]:
    merged = []
    for start, end in sorted(intervals):
        if end <= start:
            continue
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(end, merged[-1][1]))
        else:
            merged.append((start, end))
    return merged


def subtract_intervals(intervals, masks):
    """Subtract a union of masks from a union of intervals."""
    masks = merge_intervals(masks)
    result, j = [], 0
    for start, end in merge_intervals(intervals):
        while j < len(masks) and masks[j][1] <= start:
            j += 1
        cursor, k = start, j
        while k < len(masks) and masks[k][0] < end:
            a, b = masks[k]
            if a > cursor:
                result.append((cursor, min(a, end)))
            cursor = max(cursor, b)
            if cursor >= end:
                break
            k += 1
        if cursor < end:
            result.append((cursor, end))
    return result


def interval_length(intervals):
    return sum(b - a for a, b in intervals)


def portions(start, end, width):
    """Yield (tile_start, overlapping bases) in fixed physical-coordinate tiles."""
    p = (start // width) * width
    while p < end:
        yield p, min(end, p + width) - max(start, p)
        p += width


@dataclass(frozen=True)
class Signature:
    signature_id: str
    assembly_id: str
    contig_id: str
    start0: int
    end0: int
    species_taxid: str
    genus_taxid: str = ""
    species_name: str = ""
    contig_length: Optional[int] = None
    topology: str = "unknown"


class Manifest:
    """Complete eligible signatures; never infer opportunities from observed hits."""
    def __init__(self, signatures):
        self.signatures = {}
        self.contigs = {}
        self.assemblies = {}
        self.opportunities = defaultdict(list)
        species_genus = {}
        for s in signatures:
            if not all((s.signature_id, s.assembly_id, s.contig_id, s.species_taxid)):
                raise ValueError("Empty signature, assembly, contig, or species identifier")
            if s.start0 < 0 or s.end0 <= s.start0:
                raise ValueError(f"Invalid interval: {s.signature_id}")
            if s.topology not in {"linear", "circular", "unknown"}:
                raise ValueError(f"Invalid topology: {s.topology}")
            if s.contig_length is not None and s.contig_length < s.end0:
                raise ValueError(f"Signature exceeds contig length: {s.signature_id}")
            if s.topology == "circular" and s.contig_length is None:
                raise ValueError("Circular topology requires an explicit contig_length")
            old = self.signatures.get(s.signature_id)
            if old is not None:
                raise ValueError(f"Duplicate signature identifier: {s.signature_id}")
            tax = (s.species_taxid, s.genus_taxid, s.species_name)
            if s.assembly_id in self.assemblies and self.assemblies[s.assembly_id] != tax:
                raise ValueError(f"Inconsistent taxonomy for assembly {s.assembly_id}")
            if s.genus_taxid:
                if s.species_taxid in species_genus and species_genus[s.species_taxid] != s.genus_taxid:
                    raise ValueError(f"Inconsistent genus for species {s.species_taxid}")
                species_genus[s.species_taxid] = s.genus_taxid
            self.assemblies[s.assembly_id] = tax
            key = (s.assembly_id, s.contig_id)
            info = (s.contig_length, s.topology)
            if key in self.contigs and self.contigs[key] != info:
                raise ValueError(f"Inconsistent contig metadata for {key}")
            self.contigs[key] = info
            self.signatures[s.signature_id] = s
            self.opportunities[key].append((s.start0, s.end0))
        if not self.signatures:
            raise ValueError("The signature manifest is empty")
        self.opportunities = {k: merge_intervals(v) for k, v in self.opportunities.items()}

    @classmethod
    def load(cls, path):
        signatures = []
        for r in read_tsv(path, MANIFEST_COLUMNS[:6]):
            signatures.append(Signature(
                signature_id=r["signature_id"], assembly_id=r["assembly_id"], contig_id=r["contig_id"],
                start0=int(r["start0"]), end0=int(r["end0"]), species_taxid=r["species_taxid"],
                genus_taxid=r.get("genus_taxid", ""), species_name=r.get("species_name", ""),
                contig_length=int(r["contig_length"]) if r.get("contig_length") else None,
                topology=r.get("topology") or "unknown"))
        return cls(signatures)


def parse_signature_name(name, coordinates="1-based-inclusive", assembly_field=None):
    """Parse name|start|end|assembly, optionally with extra taxon fields.

    assembly_field is 1-based. Without it, prefer a unique GCF_/GCA_ field;
    otherwise accept exactly four nonempty fields. No taxid is guessed.
    """
    fields = name.rstrip("|").split("|")
    if len(fields) < 4:
        raise ValueError(f"Cannot parse signature name {name!r}; supply an explicit manifest")
    try:
        start, end = int(fields[1]), int(fields[2])
    except ValueError as exc:
        raise ValueError(f"Non-integer coordinates in {name!r}") from exc
    if coordinates == "1-based-inclusive":
        start -= 1
    elif coordinates != "0-based-half-open":
        raise ValueError(f"Unknown coordinate convention: {coordinates}")
    if start < 0 or end <= start:
        raise ValueError(f"Invalid signature interval in {name!r}")
    if assembly_field is not None:
        if assembly_field < 1 or assembly_field > len(fields):
            raise ValueError(f"assembly_field out of range for {name!r}")
        assembly = fields[assembly_field - 1]
    else:
        assemblies = [p for p in fields[3:] if re.fullmatch(r"GC[AF]_\d+(?:\.\d+)?", p)]
        if len(assemblies) == 1:
            assembly = assemblies[0]
        elif len(fields) == 4:
            assembly = fields[3]
        else:
            raise ValueError(f"Ambiguous assembly field in {name!r}; specify --assembly-field")
    return fields[0], start, end, assembly


def fasta_entries(path):
    """Stream (identifier, sequence_length) from FASTA or a FASTA .fai index."""
    if str(path).endswith(".fai"):
        with text_open(path) as handle:
            for line in handle:
                cols = line.rstrip().split("\t")
                if len(cols) < 2:
                    raise ValueError("Malformed FASTA index")
                yield cols[0], int(cols[1])
        return
    name, length = None, 0
    with text_open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    yield name, length
                name, length = line[1:].split()[0], 0
            elif line.strip():
                if name is None:
                    raise ValueError("Expected a FASTA header")
                length += len("".join(line.split()))
    if name is not None:
        yield name, length


def build_manifest(signatures_path, output, assembly_table=None, contig_table=None,
                   coordinates="1-based-inclusive", assembly_field=None, taxonomy_file=None):
    taxa, contigs = {}, {}
    if assembly_table:
        for row in read_tsv(assembly_table, ("assembly_id", "species_taxid")):
            key = row["assembly_id"]
            if key in taxa:
                raise ValueError(f"Duplicate assembly in taxonomy table: {key}")
            taxa[key] = row
    tax_module = None
    if taxonomy_file:
        try:
            from gottcha.utils import taxonomy as tax_module
        except ImportError as exc:
            raise ValueError("--taxonomy requires the installed GOTTCHA2 taxonomy module; use --assembly-table instead") from exc
        tax_module.loadTaxonomy(cus_taxonomy_file=str(taxonomy_file), auto_download=False)
    if not assembly_table and not taxonomy_file:
        raise ValueError("Provide --assembly-table or --taxonomy; assembly IDs are not species IDs")
    if contig_table:
        for row in read_tsv(contig_table, ("assembly_id", "contig_id", "contig_length")):
            key = (row["assembly_id"], row["contig_id"])
            if key in contigs:
                raise ValueError(f"Duplicate contig metadata: {key}")
            contigs[key] = (int(row["contig_length"]), row.get("topology") or "unknown")
    result = []
    for name, length in fasta_entries(signatures_path):
        contig, start, end, assembly = parse_signature_name(name, coordinates, assembly_field)
        if length != end - start:
            raise ValueError(f"FASTA length {length} disagrees with coordinates for {name}")
        if assembly not in taxa and tax_module is not None:
            species = tax_module.taxid2taxidOnRank(assembly, target_rank="species")
            genus = tax_module.taxid2taxidOnRank(assembly, target_rank="genus")
            if species and str(species) not in {"unknown", "0", "None"}:
                taxa[assembly] = dict(species_taxid=str(species), genus_taxid=str(genus or ""),
                                      species_name=tax_module.taxid2name(species))
        if assembly not in taxa:
            raise ValueError(f"No species mapping for assembly {assembly}; provide an assembly table")
        tax = taxa[assembly]
        contig_length, topology = contigs.get((assembly, contig), (None, "unknown"))
        result.append(Signature(name, assembly, contig, start, end, tax["species_taxid"],
                                tax.get("genus_taxid", ""), tax.get("species_name", ""),
                                contig_length, topology))
    manifest = Manifest(result)
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    write_tsv(output, (asdict(x) for x in manifest.signatures.values()), MANIFEST_COLUMNS)
    metadata = {"module_version": VERSION, "coordinates": "0-based-half-open",
                "source_coordinate_convention": coordinates,
                "signature_source": str(signatures_path), "signature_count": len(result),
                "opportunity_scope": "all entries in supplied signature source; verify it matches the mapped database"}
    Path(str(output) + ".json").write_text(json.dumps(metadata, indent=2) + "\n")
    return manifest


@dataclass
class Config:
    min_aligned_bases: int = 50
    min_identity: float = 0.0
    min_mapq: int = 0
    include_secondary: bool = True
    include_supplementary: bool = True
    include_duplicates: bool = False
    include_qcfail: bool = False
    window_sizes: tuple[int, ...] = (5000, 20000, 100000)
    min_link_molecules: int = 3
    min_link_fraction: float = 0.0
    same_genus_only: bool = True
    alternative_overlap: float = 0.5
    max_query_overlap: int = 20
    gap_tolerance: int = 100
    gap_relative_tolerance: float = 0.15
    pair_orientation: str = "FR"
    min_insert: Optional[int] = None
    max_insert: Optional[int] = None
    chunk_step: Optional[int] = None
    explainers: tuple[str, ...] = ()
    max_alignments_per_molecule: int = 2000
    html_link_limit: int = 2000
    min_evidence_molecules: int = 10

    def validate(self):
        for name in ("min_aligned_bases", "min_mapq", "gap_tolerance", "max_query_overlap", "html_link_limit"):
            if getattr(self, name) < 0:
                raise ValueError(f"{name} must be nonnegative")
        if self.min_mapq > 254:
            raise ValueError("min_mapq must be <=254; MAPQ 255 denotes unavailable")
        for name in ("min_identity", "min_link_fraction", "alternative_overlap", "gap_relative_tolerance"):
            if not 0 <= getattr(self, name) <= 1:
                raise ValueError(f"{name} must be between 0 and 1")
        if self.alternative_overlap <= 0:
            raise ValueError("alternative_overlap must be positive")
        if not self.window_sizes or any(w <= 0 for w in self.window_sizes):
            raise ValueError("window_sizes must contain positive integers")
        if len(set(self.window_sizes)) != len(self.window_sizes):
            raise ValueError("window_sizes cannot contain duplicates")
        if self.min_link_molecules < 1 or self.max_alignments_per_molecule < 2 or self.min_evidence_molecules < 1:
            raise ValueError("Invalid molecule count configuration")
        if self.chunk_step is not None and self.chunk_step < 1:
            raise ValueError("chunk_step must be positive")
        if self.pair_orientation not in {"FR", "RF", "FF", "any"}:
            raise ValueError("pair_orientation must be FR, RF, FF, or any")
        if self.min_insert is not None and self.min_insert < 0:
            raise ValueError("min_insert must be nonnegative")
        if self.max_insert is not None and self.max_insert < (self.min_insert or 0):
            raise ValueError("max_insert must be >= min_insert")


@dataclass
class Alignment:
    template: str
    qname: str
    read_end: int
    flag: int
    signature_id: str
    assembly: str
    contig: str
    species: str
    genus: str
    start: int
    end: int
    qstart: int
    qend: int
    query_length: int
    reverse: bool
    mapq: int
    score: Optional[int]
    identity: Optional[float]
    # Original-read query interval and reference interval for each M/= /X block.
    # The two intervals have equal length; reverse says they run oppositely.
    segments: list

    @property
    def secondary(self):
        return bool(self.flag & 0x100)

    @property
    def supplementary(self):
        return bool(self.flag & 0x800)

    @property
    def primary(self):
        return not self.secondary and not self.supplementary

    @property
    def aligned_bases(self):
        return sum(x[3] - x[2] for x in self.segments)

    def key(self):
        return (self.read_end, self.qstart, self.qend, self.assembly, self.contig,
                self.start, self.end, self.reverse, self.flag, self.signature_id,
                tuple(tuple(s) for s in self.segments))


def parse_cigar(cigar, reference_start, reverse=False, query_offset=0):
    """Return exact aligned blocks and original-orientation query coordinates.

    H/S count toward the original query length. D/N advance only reference;
    I advances only query. Coverage uses M, = and X, not bounding spans.
    """
    ops = [(int(n), op) for n, op in CIGAR_RE.findall(cigar)]
    if not ops or "".join(f"{n}{op}" for n, op in ops) != cigar or any(n <= 0 for n, _ in ops):
        raise ValueError(f"Invalid CIGAR: {cigar}")
    # Reject internal clipping, while permitting the usual H,S ... S,H arrangement.
    aligned_indices = [i for i, (_, op) in enumerate(ops) if op not in "HS"]
    if not aligned_indices:
        raise ValueError(f"CIGAR has no alignment operations: {cigar}")
    lo, hi = min(aligned_indices), max(aligned_indices)
    if any(op in "HS" for _, op in ops[lo:hi + 1]):
        raise ValueError(f"Internal clipping in CIGAR: {cigar}")
    qlength = sum(n for n, op in ops if op in "MIS H=X".replace(" ", ""))
    qpos, rpos, segments = 0, reference_start, []
    for n, op in ops:
        if op in "M=X":
            qa, qb = (qlength - qpos - n, qlength - qpos) if reverse else (qpos, qpos + n)
            segments.append((qa + query_offset, qb + query_offset, rpos, rpos + n))
        if op in "MISH=X":
            qpos += n
        if op in "MDN=X":
            rpos += n
    if not segments:
        raise ValueError(f"CIGAR has no aligned query bases: {cigar}")
    return dict(segments=segments, qstart=min(s[0] for s in segments), qend=max(s[1] for s in segments),
                query_length=qlength, reference_end=rpos, ops=ops)


def read_sam_records(path, cram_reference=None):
    """Yield SAM fields. BAM/CRAM are streamed without requiring an index."""
    suffix = str(path).lower()
    if suffix.endswith((".bam", ".cram")):
        try:
            import pysam
        except ImportError as exc:
            raise RuntimeError("BAM/CRAM input requires pysam. Install it, or provide SAM/SAM.gz.") from exc
        kwargs = {"reference_filename": str(cram_reference)} if cram_reference else {}
        with pysam.AlignmentFile(str(path), "rc" if suffix.endswith(".cram") else "rb", **kwargs) as handle:
            for record in handle.fetch(until_eof=True):
                yield record.to_string().split("\t")
    else:
        with text_open(path) as handle:
            for line_no, line in enumerate(handle, 1):
                if not line.strip() or line.startswith("@"):
                    continue
                fields = line.rstrip("\r\n").split("\t")
                if len(fields) < 11:
                    raise ValueError(f"{path}:{line_no}: SAM record needs at least 11 fields")
                yield fields


def optional_tags(fields):
    result = {}
    for value in fields[11:]:
        parts = value.split(":", 2)
        if len(parts) != 3:
            raise ValueError(f"Malformed optional SAM tag: {value}")
        tag, kind, text = parts
        if tag in result:
            raise ValueError(f"Duplicate SAM tag {tag}")
        if tag in {"NM", "AS"} and kind not in {"c", "C", "s", "S", "i", "I"}:
            raise ValueError(f"{tag} must use an integer SAM tag type")
        result[tag] = int(text) if kind in "cCsSiI" else float(text) if kind == "f" else text
    return result


def load_read_map(path):
    result = {}
    if path:
        for r in read_tsv(path, ("qname", "template_id", "query_offset")):
            if r["qname"] in result:
                raise ValueError(f"Duplicate qname in read map: {r['qname']}")
            if int(r["query_offset"]) < 0:
                raise ValueError("query_offset must be nonnegative")
            if r.get("read_end") and int(r["read_end"]) not in (0, 1, 2):
                raise ValueError("read_end must be 0, 1 or 2")
            result[r["qname"]] = r
    return result


def parse_alignment(fields, namespace, manifest, config, read_map, counters):
    flag = int(fields[1])
    if flag < 0 or flag > 65535:
        raise ValueError("SAM FLAG outside uint16 range")
    if flag & 0x4 or fields[2] == "*":
        counters["filtered_unmapped"] += 1
        return None
    for mask, allowed, label in ((0x100, config.include_secondary, "secondary"),
                                  (0x800, config.include_supplementary, "supplementary"),
                                  (0x400, config.include_duplicates, "duplicate"),
                                  (0x200, config.include_qcfail, "qcfail")):
        if flag & mask and not allowed:
            counters["filtered_" + label] += 1
            return None
    rname = fields[2]
    if rname not in manifest.signatures:
        raise ValueError(f"Mapped reference {rname!r} is absent from the complete manifest")
    sig = manifest.signatures[rname]
    tags = optional_tags(fields)
    mapq = int(fields[4])
    if not 0 <= mapq <= 255:
        raise ValueError("MAPQ outside [0,255]")
    if config.min_mapq and (mapq == 255 or mapq < config.min_mapq):
        counters["filtered_mapq"] += 1
        return None
    qname, offset, read_end = fields[0], 0, (1 if flag & 0x40 else 2 if flag & 0x80 else 0)
    if (flag & 0x40) and (flag & 0x80):
        raise ValueError(f"Both read1 and read2 flags set for {qname}")
    template_id = qname
    if qname in read_map:
        info = read_map[qname]
        template_id, offset = info["template_id"], int(info["query_offset"])
        if info.get("read_end"):
            read_end = int(info["read_end"])
    elif config.chunk_step is not None:
        match = re.fullmatch(r"(.+)\|chunk=(\d+)", qname)
        if match:
            template_id, offset = match[1], int(match[2]) * config.chunk_step
    elif re.search(r"\|chunk=\d+$", qname):
        raise ValueError("Chunked read names detected: supply --chunk-step or an explicit --read-map")
    template = json.dumps([namespace, str(tags.get("RG", "")), template_id], separators=(",", ":"))
    pos = int(fields[3])
    if pos < 1:
        raise ValueError(f"Mapped SAM POS must be positive: {qname}")
    reverse = bool(flag & 0x10)
    parsed = parse_cigar(fields[5], sig.start0 + pos - 1, reverse, offset)
    if parsed["reference_end"] > sig.end0:
        raise ValueError(f"Alignment extends beyond signature: {qname} / {rname}")
    ops = parsed["ops"]
    aligned = sum(n for n, op in ops if op in "M=X")
    if aligned < config.min_aligned_bases:
        counters["filtered_aligned_bases"] += 1
        return None
    denominator = sum(n for n, op in ops if op in "MID=X")
    if "NM" in tags:
        nm = int(tags["NM"])
        known_edits = sum(n for n, op in ops if op in "IDX")
        if nm < known_edits or nm > denominator:
            raise ValueError(f"NM inconsistent with CIGAR: {qname}")
        identity = (denominator - nm) / denominator
    elif not any(op == "M" for _, op in ops):
        identity = sum(n for n, op in ops if op == "=") / denominator
    else:
        identity = None
    if config.min_identity and (identity is None or identity < config.min_identity):
        counters["filtered_identity_missing" if identity is None else "filtered_identity"] += 1
        return None
    if identity is None:
        counters["retained_identity_unavailable"] += 1
    if "SA" in tags:
        counters["retained_records_with_SA_tag"] += 1
    # Explicit SA records are analyzed; SA tags are not synthesized into alignments.
    counters["retained_secondary" if flag & 0x100 else "retained_supplementary" if flag & 0x800 else "retained_primary"] += 1
    return Alignment(template, qname, read_end, flag, rname, sig.assembly_id, sig.contig_id,
                     sig.species_taxid, sig.genus_taxid, sig.start0 + pos - 1, parsed["reference_end"],
                     parsed["qstart"], parsed["qend"], parsed["query_length"], reverse, mapq,
                     int(tags["AS"]) if "AS" in tags else None, identity, parsed["segments"])


def deduplicate_alignments(alignments, max_placements=None, counters=None):
    """One copy per placement, preferring primary/nonsecondary records."""
    chosen, query_lengths = {}, {}
    for a in alignments:
        query_key = (a.qname, a.read_end)
        if query_key in query_lengths and query_lengths[query_key] != a.query_length:
            raise ValueError(f"Inconsistent original query length for {a.qname}, read end {a.read_end}; "
                             "check supplementary hard/soft clipping or provide distinct chunk IDs and offsets")
        query_lengths[query_key] = a.query_length
        key = (a.read_end, a.assembly, a.contig, a.start, a.end, a.reverse,
               tuple(tuple(x) for x in a.segments))
        priority = (a.secondary, a.supplementary, -a.mapq,
                    -(a.score if a.score is not None else -10**12), a.signature_id, a.qname)
        if key in chosen and counters is not None:
            counters["duplicate_placements_removed"] += 1
        if key not in chosen or priority < chosen[key][0]:
            chosen[key] = (priority, a)
        if max_placements is not None and len(chosen) > max_placements:
            raise ValueError(f"Molecule {a.template} exceeds max_alignments_per_molecule. "
                             "Increase the limit or inspect candidate-reporting settings. "
                             "No silent truncation is performed.")
    return sorted((v[1] for v in chosen.values()), key=lambda a: a.key())


def best_query_chain(alignments, max_overlap=20):
    """Deterministic maximum-covered-query chain of nonsecondary placements.

    This is an observed nonsecondary chain, NOT source inference. Tie breaks are
    lexical; alternatives remain in ambiguity tracks. An overlap must satisfy
    both the absolute limit and 20% of the shorter query interval.
    """
    items = sorted((a for a in alignments if not a.secondary), key=lambda a: a.key())
    if not items:
        return []
    scores, predecessor = [], []
    for i, current in enumerate(items):
        best, parent = current.aligned_bases, -1
        for j in range(i):
            previous = items[j]
            if previous.qstart >= current.qstart or previous.qend >= current.qend:
                continue
            overlap = max(0, previous.qend - current.qstart)
            allowed = min(max_overlap, 0.2 * min(previous.qend - previous.qstart, current.qend - current.qstart))
            if overlap > allowed:
                continue
            score = scores[j] + current.aligned_bases - overlap
            if score > best:
                best, parent = score, j
        scores.append(best)
        predecessor.append(parent)
    index = max(range(len(items)), key=lambda i: (scores[i], -i))
    chain = []
    while index >= 0:
        chain.append(items[index])
        index = predecessor[index]
    return chain[::-1]


def overlap_bases(a, b):
    aa = merge_intervals((x[0], x[1]) for x in a.segments)
    bb = merge_intervals((x[0], x[1]) for x in b.segments)
    return interval_length(aa) - interval_length(subtract_intervals(aa, bb))


def overlapping_pairs(items):
    """Sweep potential query overlaps; avoid all-pairs scans of long chunks."""
    active = []
    for a in sorted(items, key=lambda x: (x.qstart, x.qend, x.key())):
        active = [b for b in active if b.qend > a.qstart]
        for b in active:
            yield b, a
        active.append(a)


def comparable_contig(a, b, manifest):
    if a.assembly != b.assembly:
        return None, "different_assembly"
    if a.contig != b.contig:
        return None, "different_contig"
    return manifest.contigs[(a.assembly, a.contig)], ""


def split_geometry(a, b, manifest, config):
    info, reason = comparable_contig(a, b, manifest)
    qgap = b.qstart - a.qend
    result = dict(status="not_evaluable", reason=reason, query_gap=qgap,
                  genome_gap=None, gap_error=None, insert_size=None, wraps_origin=False)
    if info is None:
        return result
    length, topology = info
    if a.reverse != b.reverse:
        return dict(result, status="discordant", reason="strand_switch")
    gap = a.start - b.end if a.reverse else b.start - a.end
    candidates = [(gap, False)]
    if topology == "circular":
        candidates += [(gap + length, True), (gap - length, True)]
    gap, wrap = min(candidates, key=lambda x: (abs(x[0] - qgap), x[1]))
    error = abs(gap - qgap)
    tolerance = config.gap_tolerance + config.gap_relative_tolerance * max(abs(qgap), 1)
    ok = error <= tolerance
    status = "concordant" if ok else "not_evaluable" if topology == "unknown" else "discordant"
    reason = "gap_consistent" if ok else "gap_mismatch_topology_unknown" if topology == "unknown" else "gap_mismatch"
    return dict(result, status=status, reason=reason, genome_gap=gap, gap_error=error, wraps_origin=wrap)


def pair_geometry(a, b, manifest, config):
    info, reason = comparable_contig(a, b, manifest)
    result = dict(status="not_evaluable", reason=reason, query_gap=None,
                  genome_gap=None, gap_error=None, insert_size=None, wraps_origin=False)
    if info is None:
        return result
    length, topology = info
    shifts = [0, -length, length] if topology == "circular" else [0]
    candidates = []
    for shift in shifts:
        left, right = sorted([(a.start, a.end, a.reverse), (b.start + shift, b.end + shift, b.reverse)])
        span = max(left[1], right[1]) - min(left[0], right[0])
        orient = config.pair_orientation
        orient_ok = (orient == "any" or
                     (orient == "FR" and not left[2] and right[2]) or
                     (orient == "RF" and left[2] and not right[2]) or
                     (orient == "FF" and left[2] == right[2]))
        size_ok = span >= (config.min_insert or 0) and (config.max_insert is None or span <= config.max_insert)
        candidates.append((not orient_ok, not size_ok, span, shift != 0))
    bad_orientation, bad_size, span, wrap = min(candidates)
    if not bad_orientation and not bad_size:
        status = "concordant" if config.max_insert is not None else "orientation_only"
        reason = "fragment_consistent" if config.max_insert is not None else "insert_bounds_not_supplied"
    else:
        status = "not_evaluable" if topology == "unknown" else "discordant"
        reason = "orientation_mismatch" if bad_orientation else "insert_size_mismatch"
        if topology == "unknown":
            reason += "_topology_unknown"
    return dict(result, status=status, reason=reason, insert_size=span, wraps_origin=wrap)


LINK_COLUMNS = ["molecule_id", "kind", "species_a", "species_b", "assembly_a", "contig_a",
                "start_a", "end_a", "query_start_a", "query_end_a", "read_end_a", "strand_a",
                "assembly_b", "contig_b", "start_b", "end_b", "query_start_b", "query_end_b",
                "read_end_b", "strand_b", "status", "reason", "query_gap", "genome_gap",
                "gap_error", "insert_size", "wraps_origin", "overlap_bases"]


def link_record(a, b, kind, manifest, config, molecule_id):
    r = {"molecule_id": molecule_id, "kind": kind, "overlap_bases": None}
    for suffix, aln in (("a", a), ("b", b)):
        r.update({f"species_{suffix}": aln.species, f"assembly_{suffix}": aln.assembly,
                  f"contig_{suffix}": aln.contig, f"start_{suffix}": aln.start, f"end_{suffix}": aln.end,
                  f"query_start_{suffix}": aln.qstart, f"query_end_{suffix}": aln.qend,
                  f"read_end_{suffix}": aln.read_end, f"strand_{suffix}": "-" if aln.reverse else "+"})
    if kind == "split":
        r.update(split_geometry(a, b, manifest, config))
    elif kind == "pair":
        r.update(pair_geometry(a, b, manifest, config))
    else:
        r.update(status="not_applicable", reason="overlapping_query_placements", query_gap=None,
                 genome_gap=None, gap_error=None, insert_size=None, wraps_origin=False,
                 overlap_bases=overlap_bases(a, b))
    return r


def molecule_links(alignments, manifest, config, molecule_id):
    by_end = defaultdict(list)
    for a in alignments:
        by_end[a.read_end].append(a)
    for read_end in sorted(by_end):
        for a, b in overlapping_pairs(by_end[read_end]):
            if a.species == b.species:
                continue
            overlap = overlap_bases(a, b)
            if overlap >= config.alternative_overlap * min(a.aligned_bases, b.aligned_bases):
                yield link_record(a, b, "alternative", manifest, config, molecule_id), a, b
    chains = {end: best_query_chain(items, config.max_query_overlap) for end, items in by_end.items()}
    for end in sorted(chains):
        for a, b in zip(chains[end], chains[end][1:]):
            yield link_record(a, b, "split", manifest, config, molecule_id), a, b
    # Different read ends are physical links, not competing placements.
    for a in chains.get(1, []):
        for b in chains.get(2, []):
            yield link_record(a, b, "pair", manifest, config, molecule_id), a, b


def project_unmasked(aln, masks):
    """Project query bases left after masking, preserving indels and strand."""
    out = []
    masks = merge_intervals(masks)
    for qa, qb, ra, rb in aln.segments:
        for x, y in subtract_intervals([(qa, qb)], masks):
            if aln.reverse:
                out.append((rb - (y - qa), rb - (x - qa)))
            else:
                out.append((ra + x - qa, ra + y - qa))
    return out


def molecule_layers(alignments, explainers=()):
    """All-compatible, query-exclusive, and optional conditional residual bases.

    Query-exclusive means not aligned to another retained species on the same
    read end. It is NOT reference-independent species uniqueness. Conditional
    masks use only user-specified explainer species; no explainer is inferred.
    """
    query = defaultdict(list)
    for a in alignments:
        query[(a.read_end, a.species)].extend((s[0], s[1]) for s in a.segments)
    query = {k: merge_intervals(v) for k, v in query.items()}
    masks, residual_masks = {}, {}
    for read_end, species in query:
        masks[(read_end, species)] = merge_intervals(
            p for (end, sp), spans in query.items() if end == read_end and sp != species for p in spans)
        if explainers:
            residual_masks[(read_end, species)] = merge_intervals(
                p for (end, sp), spans in query.items()
                if end == read_end and sp != species and sp in explainers for p in spans)
    layers = defaultdict(list)
    for a in alignments:
        layers[(a.assembly, a.contig, "all")].extend((s[2], s[3]) for s in a.segments)
        layers[(a.assembly, a.contig, "query_exclusive")].extend(project_unmasked(a, masks[(a.read_end, a.species)]))
        if explainers:
            layers[(a.assembly, a.contig, "conditional_residual")].extend(project_unmasked(a, residual_masks[(a.read_end, a.species)]))
    return {key: merge_intervals(spans) for key, spans in layers.items() if spans}


def wilson_lower(successes, trials):
    if not trials:
        return 0.0
    z, p = 1.959963984540054, successes / trials
    z2 = z * z
    return (p + z2 / (2 * trials) - z * math.sqrt((p * (1 - p) + z2 / (4 * trials)) / trials)) / (1 + z2 / trials)


def build_groups(species_counts, genera, edges, config):
    parent = {s: s for s in species_counts}
    def root(s):
        while parent[s] != s:
            parent[s] = parent[parent[s]]
            s = parent[s]
        return s
    rows = []
    for (a, b), counts in sorted(edges.items()):
        n, na, nb = counts["any"], species_counts[a], species_counts[b]
        fa, fb = n / na, n / nb
        same_genus = bool(genera.get(a)) and genera.get(a) == genera.get(b)
        kept = (n >= config.min_link_molecules and max(fa, fb) >= config.min_link_fraction and
                (same_genus or not config.same_genus_only))
        if kept:
            ra, rb = root(a), root(b)
            parent[max(ra, rb)] = min(ra, rb)
        rows.append(dict(species_a=a, species_b=b, any_molecules=n,
                         alternative_molecules=counts["alternative"], split_molecules=counts["split"],
                         pair_molecules=counts["pair"], a_molecules=na, b_molecules=nb,
                         a_link_fraction=fa, b_link_fraction=fb,
                         a_link_wilson_lower=wilson_lower(n, na), b_link_wilson_lower=wilson_lower(n, nb),
                         a_primary_to_b_alternative=counts["a_primary_to_b"],
                         b_primary_to_a_alternative=counts["b_primary_to_a"],
                         same_genus=same_genus, kept=kept))
    groups = defaultdict(list)
    for sp in sorted(parent):
        groups[root(sp)].append(sp)
    groups = {f"group{i + 1}": values for i, (_, values) in enumerate(sorted(groups.items()))}
    return groups, rows


def init_spool(path):
    db = sqlite3.connect(path)
    db.executescript("""
        PRAGMA journal_mode=OFF;
        PRAGMA synchronous=OFF;
        PRAGMA temp_store=FILE;
        PRAGMA cache_size=-65536;
        CREATE TABLE alignments (template TEXT NOT NULL, data TEXT NOT NULL);
        CREATE TABLE blocks (assembly TEXT, contig TEXT, layer TEXT, start INTEGER, end INTEGER);
        CREATE TABLE hits (assembly TEXT, contig TEXT, layer TEXT, width INTEGER, start INTEGER,
                           molecule TEXT, bases INTEGER);
        CREATE TABLE memberships (assembly TEXT, layer TEXT, molecule TEXT);
    """)
    return db


def stream_merged(rows):
    current = None
    for start, end in rows:
        if current is None:
            current = (start, end)
        elif start <= current[1]:
            current = (current[0], max(end, current[1]))
        else:
            yield current
            current = (start, end)
    if current is not None:
        yield current


def file_sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


@dataclass
class AnalysisResult:
    output_dir: Path
    groups: dict
    assembly_summary: list
    species_summary: list
    counters: dict
    files: dict = field(default_factory=dict)


def load_occupancy_model(path, manifest, widths):
    """Optional empirically calibrated P(molecule touches window | assembly hit).

    Probabilities need not sum to 1: one physical molecule can touch many windows.
    A supplied model must cover every eligible window for any modeled assembly /
    width combination. It is used for the all-compatible layer only.
    """
    result = {}
    if not path:
        return result
    for r in read_tsv(path, ("assembly_id", "contig_id", "window_size", "start0", "touch_probability")):
        key = (r["assembly_id"], r["contig_id"], int(r["window_size"]), int(r["start0"]))
        probability = float(r["touch_probability"])
        if key in result or not 0 <= probability <= 1 or key[2] not in widths or key[3] < 0 or key[3] % key[2]:
            raise ValueError(f"Invalid or duplicate occupancy model row: {key}")
        if key[:2] not in manifest.opportunities:
            raise ValueError(f"Occupancy model contains an unknown assembly/contig: {key}")
        result[key] = probability
    for assembly, width in sorted({(key[0], key[2]) for key in result}):
        expected = {(a, c, width, p) for (a, c), spans in manifest.opportunities.items() if a == assembly
                    for start, end in spans for p, _ in portions(start, end, width)}
        actual = {k for k in result if k[0] == assembly and k[2] == width}
        if actual != expected:
            raise ValueError(f"Occupancy model must cover every eligible window for {assembly}, width {width}")
    return result


def summarize(db, manifest, config, geometries, identity_stats, observed_assemblies, occupancy_model,
              coverage_writer):
    """Summarize per assembly, never pool coordinate systems across strains."""
    all_windows, assembly_rows, tracks = [], [], []
    layers = ["all", "query_exclusive"] + (["conditional_residual"] if config.explainers else [])
    memberships = {(a, layer): n for a, layer, n in db.execute(
        "SELECT assembly,layer,COUNT(*) FROM memberships GROUP BY assembly,layer")}
    hits = {(a, c, layer, w, p): (n, bases) for a, c, layer, w, p, n, bases in db.execute(
        "SELECT assembly,contig,layer,width,start,COUNT(*),SUM(bases) FROM hits GROUP BY assembly,contig,layer,width,start")}
    for assembly in sorted(observed_assemblies):
        species, genus, name = manifest.assemblies[assembly]
        contig_keys = sorted(k for k in manifest.contigs if k[0] == assembly)
        opportunities = {}
        covered_by_window = Counter()
        covered_total = Counter()
        signature_total = 0
        for key in contig_keys:
            _, contig = key
            signature_total += interval_length(manifest.opportunities[key])
            for width in config.window_sizes:
                for start, end in manifest.opportunities[key]:
                    for p, bases in portions(start, end, width):
                        opportunities[(contig, width, p)] = opportunities.get((contig, width, p), 0) + bases
            coverage_preview, preview_truncated = {}, set()
            for layer in layers:
                coverage_preview[layer] = []
                intervals = stream_merged(db.execute(
                    "SELECT start,end FROM blocks WHERE assembly=? AND contig=? AND layer=? ORDER BY start,end",
                    (assembly, contig, layer)))
                for start, end in intervals:
                    if len(coverage_preview[layer]) < 5000:
                        coverage_preview[layer].append((start, end))
                    else:
                        preview_truncated.add(layer)
                    covered_total[layer] += end - start
                    coverage_writer.writerow(dict(assembly_id=assembly, contig_id=contig, layer=layer, start0=start, end0=end))
                    for width in config.window_sizes:
                        for p, bases in portions(start, end, width):
                            covered_by_window[(contig, layer, width, p)] += bases
            length, topology = manifest.contigs[key]
            tracks.append(dict(assembly_id=assembly, species_taxid=species, species_name=name, contig_id=contig,
                               contig_length=length, display_extent=length or manifest.opportunities[key][-1][1],
                               extent_is_signature_bound=length is None, topology=topology,
                               signatures=manifest.opportunities[key], coverage=coverage_preview,
                               coverage_preview_truncated=sorted(preview_truncated)))
        for layer in layers:
            molecule_count = memberships.get((assembly, layer), 0)
            for width in config.window_sizes:
                selected = [(c, p, bp) for (c, w, p), bp in opportunities.items() if w == width]
                bins, total_support = [], 0
                for contig, p, opportunity in sorted(selected):
                    n, bases = hits.get((assembly, contig, layer, width, p), (0, 0))
                    covered = covered_by_window[(contig, layer, width, p)]
                    if covered > opportunity:
                        raise AssertionError("Coverage exceeds signature opportunity")
                    length = manifest.contigs[(assembly, contig)][0]
                    row = dict(assembly_id=assembly, species_taxid=species, contig_id=contig, layer=layer,
                               window_size=width, start0=p, end0=min(p + width, length) if length else p + width,
                               signature_bp=opportunity, covered_bp=covered, breadth=covered / opportunity,
                               molecules=n, aligned_bp=bases, mean_signature_depth=bases / opportunity)
                    bins.append(row)
                    all_windows.append(row)
                    total_support += bases
                mean_depth = total_support / signature_total
                cv = (math.sqrt(sum(r["signature_bp"] * (r["mean_signature_depth"] - mean_depth)**2 for r in bins) /
                                signature_total) / mean_depth) if mean_depth else None
                probs = [occupancy_model.get((assembly, r["contig_id"], width, r["start0"])) for r in bins]
                expected = None
                if layer == "all" and probs and all(p is not None for p in probs):
                    expected = sum((1.0 if p == 1 and molecule_count else
                                    -math.expm1(molecule_count * math.log1p(-p)) if p < 1 else 0.0) for p in probs)
                occupied = sum(r["molecules"] > 0 for r in bins)
                identity_sum, identity_weight = identity_stats.get(assembly, (0.0, 0))
                row = dict(assembly_id=assembly, species_taxid=species, genus_taxid=genus, species_name=name,
                           layer=layer, window_size=width, molecules=molecule_count, signature_bp=signature_total,
                           covered_bp=covered_total[layer], signature_breadth=covered_total[layer] / signature_total,
                           aligned_bp=total_support, mean_signature_depth=mean_depth, callable_windows=len(bins),
                           occupied_windows=occupied, occupied_fraction=occupied / len(bins),
                           max_tiled_support_fraction=max((r["aligned_bp"] for r in bins), default=0) / total_support if total_support else None,
                           callable_depth_cv=cv, expected_occupied_windows=expected,
                           occupancy_ratio=occupied / expected if expected else None,
                           mean_alignment_identity=identity_sum / identity_weight if identity_weight and layer == "all" else None,
                           evidence_level="low_molecule_count" if molecule_count < config.min_evidence_molecules else "inspectable_not_classified")
                for kind in ("split", "pair"):
                    for status in ("concordant", "discordant", "mixed", "orientation_only", "not_evaluable"):
                        row[f"{kind}_{status}_molecules"] = geometries[(assembly, kind, status)] if layer == "all" else None
                assembly_rows.append(row)
    return assembly_rows, all_windows, tracks


def analyze_alignments(inputs, manifest_path, output_dir, config=None, read_map_path=None,
                       occupancy_model_path=None, tmpdir=None, namespaces=None,
                       cram_reference=None, overwrite=False, make_html=True):
    """Run diagnostics from SAM/BAM/CRAM paths and an explicit signature manifest.

    Input order may be arbitrary: SQLite groups complete physical molecules on
    disk. Default namespaces separate files; supply identical explicit namespaces
    only when records in different files truly belong to the same library.
    Returns AnalysisResult. No existing GOTTCHA2 profiles are modified.
    """
    config = config or Config()
    config.validate()
    inputs = [Path(x) for x in ([inputs] if isinstance(inputs, (str, Path)) else inputs)]
    if not inputs:
        raise ValueError("At least one alignment input is required")
    if len(set(str(x.resolve()) for x in inputs)) != len(inputs):
        raise ValueError("The same alignment file was supplied more than once")
    namespaces = list(namespaces) if namespaces is not None else [str(x.resolve()) for x in inputs]
    if len(namespaces) != len(inputs) or not all(namespaces):
        raise ValueError("Provide one nonempty namespace per input")
    manifest = Manifest.load(manifest_path)
    read_map = load_read_map(read_map_path)
    occupancy_model = load_occupancy_model(occupancy_model_path, manifest, config.window_sizes)
    if config.explainers and not set(config.explainers).issubset({x[0] for x in manifest.assemblies.values()}):
        raise ValueError("An explainer species is absent from the manifest")
    output_dir = Path(output_dir)
    if output_dir.exists() and any(output_dir.iterdir()) and not overwrite:
        raise FileExistsError(f"Output directory is not empty: {output_dir}; use --overwrite explicitly")
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    if tmpdir:
        Path(tmpdir).mkdir(parents=True, exist_ok=True)
    counters, species_counts, species_exclusive, species_residual = Counter(), Counter(), Counter(), Counter()
    genera = {tax[0]: tax[1] for tax in manifest.assemblies.values() if tax[1]}
    species_names = {tax[0]: tax[2] for tax in manifest.assemblies.values()}
    edges, geometries, identity_stats, observed_assemblies = defaultdict(Counter), Counter(), {}, set()
    browser_links, input_stats = [], []
    with tempfile.TemporaryDirectory(prefix="gottcha-coherence-", dir=tmpdir) as work:
        stage = Path(work) / "reports"
        stage.mkdir()
        db = init_spool(str(Path(work) / "evidence.sqlite"))
        try:
            batch = []
            for path, namespace in zip(inputs, namespaces):
                stat = path.stat()
                input_stats.append(dict(path=str(path.resolve()), namespace=namespace, size=stat.st_size, mtime_ns=stat.st_mtime_ns))
                LOG.info("Reading %s", path)
                for fields in read_sam_records(path, cram_reference):
                    counters["input_records"] += 1
                    try:
                        a = parse_alignment(fields, namespace, manifest, config, read_map, counters)
                    except (ValueError, IndexError) as exc:
                        raise ValueError(f"{path}: record {counters['input_records']} ({fields[0]}): {exc}") from exc
                    if a is not None:
                        counters["retained_records_before_deduplication"] += 1
                        batch.append((a.template, json.dumps(asdict(a), separators=(",", ":"))))
                    if len(batch) >= 5000:
                        db.executemany("INSERT INTO alignments VALUES (?,?)", batch)
                        db.commit()
                        batch.clear()
            if batch:
                db.executemany("INSERT INTO alignments VALUES (?,?)", batch)
            db.execute("CREATE INDEX template_idx ON alignments(template)")
            db.commit()
            with text_open(stage / "links.tsv.gz", "wt") as link_handle, text_open(stage / "molecules.tsv.gz", "wt") as molecule_handle:
                link_writer = csv.DictWriter(link_handle, LINK_COLUMNS, delimiter="\t", lineterminator="\n")
                link_writer.writeheader()
                molecule_columns = ["molecule_id", "namespace", "read_group", "template_id", "qnames", "species_taxids", "alignments"]
                molecule_writer = csv.DictWriter(molecule_handle, molecule_columns, delimiter="\t", lineterminator="\n")
                molecule_writer.writeheader()
                cursor = db.execute("SELECT template,data FROM alignments ORDER BY template")
                for template, records in itertools.groupby(cursor, key=lambda r: r[0]):
                    raw = (Alignment(**json.loads(row[1])) for row in records)
                    alignments = deduplicate_alignments(raw, config.max_alignments_per_molecule, counters)
                    counters["retained_placements"] += len(alignments)
                    counters["eligible_molecules"] += 1
                    molecule_id = hashlib.sha256(template.encode()).hexdigest()[:24]
                    species = sorted({a.species for a in alignments})
                    species_counts.update(species)
                    observed_assemblies.update(a.assembly for a in alignments)
                    for a in alignments:
                        if a.identity is not None:
                            s, n = identity_stats.get(a.assembly, (0.0, 0))
                            identity_stats[a.assembly] = (s + a.identity * a.aligned_bases, n + a.aligned_bases)
                    ns, rg, template_id = json.loads(template)
                    molecule_writer.writerow(dict(molecule_id=molecule_id, namespace=ns, read_group=rg,
                                                  template_id=template_id, qnames=json.dumps(sorted({a.qname for a in alignments})),
                                                  species_taxids=json.dumps(species), alignments=len(alignments)))
                    # Counts are unique physical molecules, never alignment-segment totals.
                    local_edges, local_directions, local_geometry = defaultdict(set), set(), defaultdict(set)
                    for link, a, b in molecule_links(alignments, manifest, config, molecule_id):
                        link_writer.writerow(link)
                        counters["link_records"] += 1
                        if len(browser_links) < config.html_link_limit:
                            browser_links.append(link)
                        if a.species != b.species:
                            sa, sb = sorted((a.species, b.species))
                            local_edges[(sa, sb)].add(link["kind"])
                            if link["kind"] == "alternative":
                                if a.primary:
                                    local_directions.add((a.species, b.species))
                                if b.primary:
                                    local_directions.add((b.species, a.species))
                        if link["kind"] != "alternative":
                            for assembly in {a.assembly, b.assembly}:
                                local_geometry[(assembly, link["kind"])].add(link["status"])
                    for pair, kinds in local_edges.items():
                        edges[pair]["any"] += 1
                        edges[pair].update(kinds)
                    for sa, sb in local_directions:
                        pair = tuple(sorted((sa, sb)))
                        edges[pair]["a_primary_to_b" if sa == pair[0] else "b_primary_to_a"] += 1
                    for (assembly, kind), statuses in local_geometry.items():
                        informative = statuses & {"concordant", "discordant"}
                        status = "mixed" if len(informative) == 2 else next(iter(informative)) if informative else \
                                 "orientation_only" if "orientation_only" in statuses else "not_evaluable"
                        geometries[(assembly, kind, status)] += 1
                    layers = molecule_layers(alignments, set(config.explainers))
                    memberships, tile_bases = set(), Counter()
                    exclusive_species, residual_species = set(), set()
                    for (assembly, contig, layer), spans in layers.items():
                        memberships.add((assembly, layer, molecule_id))
                        taxid = manifest.assemblies[assembly][0]
                        if layer == "query_exclusive":
                            exclusive_species.add(taxid)
                        elif layer == "conditional_residual":
                            residual_species.add(taxid)
                        db.executemany("INSERT INTO blocks VALUES (?,?,?,?,?)",
                                       ((assembly, contig, layer, start, end) for start, end in spans))
                        for width in config.window_sizes:
                            for start, end in spans:
                                for p, bases in portions(start, end, width):
                                    tile_bases[(assembly, contig, layer, width, p, molecule_id)] += bases
                    species_exclusive.update(exclusive_species)
                    species_residual.update(residual_species)
                    db.executemany("INSERT INTO memberships VALUES (?,?,?)", sorted(memberships))
                    db.executemany("INSERT INTO hits VALUES (?,?,?,?,?,?,?)", (key + (bases,) for key, bases in tile_bases.items()))
                    if counters["eligible_molecules"] % 10000 == 0:
                        db.commit()
                        LOG.info("Analyzed %s molecules", format(counters["eligible_molecules"], ","))
            db.execute("CREATE INDEX block_idx ON blocks(assembly,contig,layer,start,end)")
            db.commit()
            with text_open(stage / "coverage_intervals.tsv.gz", "wt") as coverage_handle:
                coverage_writer = csv.DictWriter(coverage_handle, ["assembly_id", "contig_id", "layer", "start0", "end0"],
                                                 delimiter="\t", lineterminator="\n")
                coverage_writer.writeheader()
                assembly_rows, windows, tracks = summarize(db, manifest, config, geometries, identity_stats,
                                                           observed_assemblies, occupancy_model, coverage_writer)
        finally:
            db.close()
        groups, edge_rows = build_groups(species_counts, genera, edges, config)
        membership = {s: group for group, species in groups.items() for s in species}
        species_rows = []
        for species in sorted(species_counts):
            members = [a for a in observed_assemblies if manifest.assemblies[a][0] == species]
            species_rows.append(dict(species_taxid=species, species_name=species_names.get(species, ""),
                                     genus_taxid=genera.get(species, ""), component=membership[species],
                                     molecules=species_counts[species], query_exclusive_molecules=species_exclusive[species],
                                     conditional_residual_molecules=species_residual[species] if config.explainers else None,
                                     assemblies_with_support=len(members), assembly_ids=";".join(sorted(members)),
                                     source_status="not_inferred"))
        write_tsv(stage / "assemblies.tsv", assembly_rows, ASSEMBLY_COLUMNS)
        write_tsv(stage / "windows.tsv.gz", windows, WINDOW_COLUMNS)
        write_tsv(stage / "species.tsv", species_rows, SPECIES_COLUMNS)
        write_tsv(stage / "edges.tsv", edge_rows, EDGE_COLUMNS)
        (stage / "groups.json").write_text(json.dumps(groups, indent=2) + "\n")
        warnings = ["Diagnostic evidence only: no species presence, absence, novelty, or abundance is inferred.",
                    "Query-exclusive and residual support depend on retained alignments and the candidate database; neither proves species uniqueness.",
                    "All-compatible depth is not an abundance estimate and can support multiple alternative loci or assemblies.",
                    "Signature completeness and correspondence to the mapped search space are caller responsibilities.",
                    "A graph component may contain multiple true species; components are not collapsed."]
        if config.min_mapq or not config.include_secondary or not config.include_supplementary:
            warnings.append("Evidence-retention filters may hide alternative placements or physical links.")
        if not counters["retained_secondary"]:
            warnings.append("No secondary records were retained; missing alternatives are not evidence of uniqueness.")
        if config.max_insert is None:
            warnings.append("No maximum insert size supplied: paired-end geometry can report orientation only, not full fragment concordance.")
        if not occupancy_model:
            warnings.append("No calibrated occupancy model supplied: expected occupancy and its ratio are intentionally blank.")
        if not observed_assemblies:
            warnings.append("No eligible mapped alignments remained after filtering.")
        metadata = dict(module_version=VERSION, completed=True, coordinate_convention="0-based-half-open",
                        manifest_path=str(Path(manifest_path).resolve()), manifest_sha256=file_sha256(manifest_path),
                        manifest_signature_count=len(manifest.signatures), inputs=input_stats, config=asdict(config),
                        counters=dict(counters), warnings=warnings,
                        occupancy_model_sha256=file_sha256(occupancy_model_path) if occupancy_model_path else None,
                        read_map_sha256=file_sha256(read_map_path) if read_map_path else None,
                        browser_link_count=len(browser_links), browser_links_truncated=counters["link_records"] > len(browser_links),
                        browser_link_selection="first links in deterministic molecule/placement order; tables retain all links")
        (stage / "run.json").write_text(json.dumps(metadata, indent=2) + "\n")
        if make_html:
            try:
                from .coherence_browser import write_browser
            except ImportError:
                from coherence_browser import write_browser
            write_browser(stage / "evidence.html", dict(metadata=metadata, assemblies=assembly_rows, species=species_rows,
                                                        edges=edge_rows, groups=groups, windows=windows, tracks=tracks,
                                                        links=browser_links))
        # Analysis is complete before any user-facing output files are replaced.
        output_dir.mkdir(parents=True, exist_ok=True)
        import shutil
        if not make_html and (output_dir / "evidence.html").exists():
            (output_dir / "evidence.html").unlink()
        for file in stage.iterdir():
            destination = output_dir / file.name
            # Copy to a same-filesystem temporary file before atomic replacement.
            temporary = output_dir / ("." + file.name + ".tmp")
            shutil.copyfile(file, temporary)
            os.replace(temporary, destination)
        files = {p.name: str(p) for p in output_dir.iterdir() if p.is_file()}
    return AnalysisResult(output_dir, groups, assembly_rows, species_rows, dict(counters), files)


ASSEMBLY_COLUMNS = [
    "assembly_id", "species_taxid", "genus_taxid", "species_name", "layer", "window_size", "molecules",
    "signature_bp", "covered_bp", "signature_breadth", "aligned_bp", "mean_signature_depth",
    "callable_windows", "occupied_windows", "occupied_fraction", "max_tiled_support_fraction",
    "callable_depth_cv", "expected_occupied_windows", "occupancy_ratio", "mean_alignment_identity", "evidence_level",
] + [f"{kind}_{status}_molecules" for kind in ("split", "pair")
     for status in ("concordant", "discordant", "mixed", "orientation_only", "not_evaluable")]
WINDOW_COLUMNS = ["assembly_id", "species_taxid", "contig_id", "layer", "window_size", "start0", "end0",
                  "signature_bp", "covered_bp", "breadth", "molecules", "aligned_bp", "mean_signature_depth"]
SPECIES_COLUMNS = ["species_taxid", "species_name", "genus_taxid", "component", "molecules", "query_exclusive_molecules",
                   "conditional_residual_molecules", "assemblies_with_support", "assembly_ids", "source_status"]
EDGE_COLUMNS = ["species_a", "species_b", "any_molecules", "alternative_molecules", "split_molecules", "pair_molecules",
                "a_molecules", "b_molecules", "a_link_fraction", "b_link_fraction", "a_link_wilson_lower", "b_link_wilson_lower",
                "a_primary_to_b_alternative", "b_primary_to_a_alternative", "same_genus", "kept"]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--version", action="version", version=VERSION)
    sub = parser.add_subparsers(dest="command", required=True)
    build = sub.add_parser("build-manifest", help="Create the complete signature opportunity manifest")
    build.add_argument("--signatures", required=True, help="Actual signature FASTA/FASTA.gz or .fai used for mapping")
    taxonomy_args = build.add_mutually_exclusive_group(required=True)
    taxonomy_args.add_argument("--assembly-table", help="TSV: assembly_id, species_taxid, optional genus_taxid/species_name")
    taxonomy_args.add_argument("--taxonomy", help="GOTTCHA2 taxonomy file; assembly IDs must be resolvable as taxonomy nodes")
    build.add_argument("--contigs", help="TSV: assembly_id, contig_id, contig_length, optional topology")
    build.add_argument("--coordinates", choices=["1-based-inclusive", "0-based-half-open"], default="1-based-inclusive")
    build.add_argument("--assembly-field", type=int, help="1-based RNAME field containing assembly ID")
    build.add_argument("-o", "--output", required=True)
    run = sub.add_parser("analyze", help="Analyze molecule links and assembly-resolved signature support")
    run.add_argument("-i", "--input", nargs="+", required=True, help="SAM/SAM.gz/BAM/CRAM; arbitrary alignment order")
    run.add_argument("--manifest", required=True)
    run.add_argument("-o", "--outdir", required=True)
    run.add_argument("--namespace", nargs="+", help="One library namespace per input; default separates files")
    run.add_argument("--read-map", help="TSV: qname, template_id, query_offset, optional read_end")
    run.add_argument("--chunk-step", type=int, help="Original read offsets for GOTTCHA2 |chunk=N names: N * step")
    run.add_argument("--cram-reference", help="Signature FASTA required to decode reference-dependent CRAM")
    run.add_argument("--min-aligned-bases", type=int, default=50, help="M/= /X bases, not AS or CIGAR bounding span (default: 50)")
    run.add_argument("--min-identity", type=float, default=0.0, help="NM/CIGAR identity; unknown identity fails any positive threshold")
    run.add_argument("--min-mapq", type=int, default=0, help="Default 0 retains ambiguous evidence; 255 is not high confidence")
    run.add_argument("--exclude-secondary", action="store_true")
    run.add_argument("--exclude-supplementary", action="store_true")
    run.add_argument("--include-duplicates", action="store_true")
    run.add_argument("--include-qcfail", action="store_true")
    run.add_argument("--window-sizes", default="5000,20000,100000", help="Comma-separated physical window sizes")
    run.add_argument("--min-link-molecules", type=int, default=3)
    run.add_argument("--min-link-fraction", type=float, default=0.0,
                     help="Grouping requires this fraction at either endpoint; physical links need no reciprocal primary assignment")
    run.add_argument("--allow-cross-genus", action="store_true", help="Permit cross-genus grouping; all raw links are always reported")
    run.add_argument("--alternative-overlap", type=float, default=0.5, help="Required overlap / shorter aligned query length")
    run.add_argument("--max-query-overlap", type=int, default=20)
    run.add_argument("--gap-tolerance", type=int, default=100)
    run.add_argument("--gap-relative-tolerance", type=float, default=0.15)
    run.add_argument("--pair-orientation", choices=["FR", "RF", "FF", "any"], default="FR")
    run.add_argument("--min-insert", type=int)
    run.add_argument("--max-insert", type=int, help="Required for full paired-fragment concordance; use library-specific bounds")
    run.add_argument("--explainer-species", default="", help="Comma-separated species IDs for optional conditional query-overlap subtraction")
    run.add_argument("--occupancy-model", help="Empirical per-window molecule-touch probabilities; no model is guessed")
    run.add_argument("--max-alignments-per-molecule", type=int, default=2000)
    run.add_argument("--min-evidence-molecules", type=int, default=10, help="Evidence-size annotation, NOT a species detection cutoff")
    run.add_argument("--html-link-limit", type=int, default=2000, help="Display limit only; links.tsv.gz remains complete")
    run.add_argument("--no-html", action="store_true")
    run.add_argument("--tmpdir", help="Scratch disk for SQLite molecule grouping")
    run.add_argument("--overwrite", action="store_true")
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    try:
        if args.command == "build-manifest":
            manifest = build_manifest(args.signatures, args.output, args.assembly_table, args.contigs,
                                      args.coordinates, args.assembly_field, args.taxonomy)
            LOG.info("Wrote %d signatures to %s", len(manifest.signatures), args.output)
        else:
            config = Config(min_aligned_bases=args.min_aligned_bases, min_identity=args.min_identity,
                            min_mapq=args.min_mapq, include_secondary=not args.exclude_secondary,
                            include_supplementary=not args.exclude_supplementary, include_duplicates=args.include_duplicates,
                            include_qcfail=args.include_qcfail, window_sizes=tuple(int(w) for w in args.window_sizes.split(",")),
                            min_link_molecules=args.min_link_molecules, min_link_fraction=args.min_link_fraction,
                            same_genus_only=not args.allow_cross_genus, alternative_overlap=args.alternative_overlap,
                            max_query_overlap=args.max_query_overlap, gap_tolerance=args.gap_tolerance,
                            gap_relative_tolerance=args.gap_relative_tolerance, pair_orientation=args.pair_orientation,
                            min_insert=args.min_insert, max_insert=args.max_insert, chunk_step=args.chunk_step,
                            explainers=tuple(x.strip() for x in args.explainer_species.split(",") if x.strip()),
                            max_alignments_per_molecule=args.max_alignments_per_molecule,
                            min_evidence_molecules=args.min_evidence_molecules, html_link_limit=args.html_link_limit)
            result = analyze_alignments(args.input, args.manifest, args.outdir, config, args.read_map,
                                        args.occupancy_model, args.tmpdir, args.namespace, args.cram_reference,
                                        args.overwrite, not args.no_html)
            LOG.info("Analyzed %d molecules; reports: %s", result.counters.get("eligible_molecules", 0), result.output_dir)
    except (OSError, ValueError, RuntimeError) as exc:
        parser.exit(2, f"ERROR: {exc}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
