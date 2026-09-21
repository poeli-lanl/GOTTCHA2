"""Discover profiling outputs and prepare data for the coverage browser."""

import gzip
import re
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional

import pysam
from pysam import bcftools

@dataclass
class BrowserInputs:
    full: Path
    output: Path
    bam: Optional[Path]
    coverage: Optional[Path]
    reference: Optional[Path]
    vcfs: List[Path]


def _input_file(value, label):
    path = Path(value).expanduser().resolve()
    if not path.is_file():
        raise ValueError(f"{label} file not found: {path}")
    return path


def _sample_prefix(bam):
    return re.sub(r"\.gottcha_[^.]+$", "", bam.stem)


def resolve_inputs(*, results=None, prefix=None, bam=None, coverage=None,
                   full=None, reference=None, vcfs=(), output=None,
                   no_variants=False):
    """Match files by sample prefix; never choose arbitrarily among samples."""
    directory = Path(results).expanduser().resolve() if results else None
    if directory is not None and not directory.is_dir():
        raise ValueError(f"Result directory not found: {directory}")
    if prefix and directory is None:
        raise ValueError("--prefix requires --results")
    if prefix and Path(prefix).name != prefix:
        raise ValueError("--prefix must be a sample name, not a path")

    bam = _input_file(bam, "BAM") if bam else None
    coverage = _input_file(coverage, "Coverage") if coverage else None
    full = _input_file(full, "Full taxonomy") if full else None
    reference = _input_file(reference, "Reference") if reference else None
    vcfs = [_input_file(vcf, "VCF") for vcf in vcfs]

    if directory is not None:
        prefix = prefix or (
            full.name.removesuffix(".full.tsv") if full else
            _sample_prefix(bam) if bam else
            re.sub(r"\.gottcha_[^.]+\.coverage\.tsv$", "", coverage.name) if coverage else None
        )
        if prefix is None:
            candidates = sorted(directory.glob("*.gottcha_*.bam"))
            if len(candidates) != 1:
                raise ValueError(
                    f"Expected one GOTTCHA2 BAM in {directory}; found {len(candidates)}. "
                    "Use --prefix to select a sample or --bam to specify its BAM."
                )
            bam = candidates[0]
            prefix = _sample_prefix(bam)
        if bam is None and coverage is None:
            candidates = sorted(
                path for path in directory.glob("*.bam")
                if _sample_prefix(path) == prefix
            )
            if len(candidates) != 1:
                raise ValueError(
                    f"Expected one BAM for sample '{prefix}'; found {len(candidates)}. "
                    "Specify --bam explicitly."
                )
            bam = candidates[0]
        if full is None:
            full = _input_file(directory / f"{prefix}.full.tsv", "Full taxonomy")
        if not no_variants:
            if not vcfs and bam is not None:
                candidate = directory / f"{bam.stem}.vcf.gz"
                if candidate.is_file():
                    vcfs = [candidate]
            if reference is None and not vcfs:
                # Prefer an already prepared reference over the original gzip.
                for suffix in ("fa.bgz", "fa.gz", "fa", "fna.bgz", "fna.gz", "fna", "fasta"):
                    candidate = directory / f"{prefix}.sylph_extracted.{suffix}"
                    if candidate.is_file():
                        reference = candidate
                        break

    if bam is None and coverage is None:
        raise ValueError("Provide --results, --bam, or an existing --coverage file")
    if full is None:
        raise ValueError("--full is required when not using --results")
    if reference is not None and bam is None and not vcfs and not no_variants:
        raise ValueError("--reference requires --bam to call variants")

    if output:
        output = Path(output).expanduser().resolve()
    elif directory is not None:
        output = directory / f"{prefix}.coverage.html"
    elif bam is not None:
        output = bam.parent / f"{_sample_prefix(bam)}.coverage.html"
    else:
        output = Path("coverage_visualization.html").resolve()
    if output.is_dir():
        raise ValueError(f"Output must be an HTML file, not a directory: {output}")
    if output in [full, bam, coverage, reference, *vcfs]:
        raise ValueError("Output HTML must not overwrite an input file")

    return BrowserInputs(full, output, bam, coverage, reference, [] if no_variants else vcfs)


def _prepare_reference(reference, destination):
    """Stream plain or gzip FASTA into an indexed BGZF copy, leaving input intact."""
    with reference.open("rb") as source:
        compressed = source.read(2) == b"\x1f\x8b"
    opener = gzip.open if compressed else open
    with opener(reference, "rb") as source, pysam.BGZFile(str(destination), "wb") as target:
        shutil.copyfileobj(source, target)
    pysam.faidx(str(destination))


def prepare_inputs(inputs, workdir, *, threads=1, no_variants=False):
    """Run bundled samtools/bcftools with the same filters as the shell workflow.

    Large intermediates stay on disk in the caller's temporary directory.
    pysam raises SamtoolsError for a failed step, so subsequent steps cannot
    silently consume partial coverage, pileup, or variant files.
    """
    workdir = Path(workdir)
    coverage = inputs.coverage
    if coverage is None:
        coverage = workdir / "coverage.tsv"
        pysam.coverage("-o", str(coverage), str(inputs.bam), catch_stdout=False)

    vcfs = list(inputs.vcfs)
    if no_variants:
        return coverage, []
    if vcfs or inputs.reference is None:
        return coverage, vcfs

    reference = workdir / "reference.fa.bgz"
    _prepare_reference(inputs.reference, reference)
    with pysam.AlignmentFile(str(inputs.bam), "rb") as bam, pysam.FastaFile(str(reference)) as fasta:
        lengths = dict(zip(fasta.references, fasta.lengths))
        mismatched = [name for name, length in zip(bam.references, bam.lengths)
                      if lengths.get(name) != length]
        if mismatched:
            raise ValueError(
                "Reference FASTA must contain the BAM's signature contigs with matching lengths; "
                f"missing or mismatched: {', '.join(mismatched[:5])}"
            )

    pileup = workdir / "pileup.bcf"
    vcf = workdir / f"{inputs.bam.stem}.vcf.gz"
    bcftools.mpileup(
        "-Ob", "-o", str(pileup), "-f", str(reference), "-q", "20", "-Q", "20",
        "-a", "FORMAT/DP,FORMAT/AD", "--threads", str(threads), str(inputs.bam),
        catch_stdout=False,
    )
    bcftools.call(
        "-mv", "--ploidy", "1", "-Oz", "--threads", str(threads),
        "-o", str(vcf), str(pileup), catch_stdout=False,
    )
    bcftools.index("-t", str(vcf), catch_stdout=False)
    return coverage, [vcf]
