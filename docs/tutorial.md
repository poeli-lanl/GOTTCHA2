[![logo](https://gottcha2.readthedocs.io/en/latest/_images/gottcha_icon.png)](https://gottcha2.readthedocs.io/en/latest/_images/gottcha_icon.png)

# Genomic Origin Through Taxonomic CHAllenge (GOTTCHA2)

GOTTCHA2 is a gene-independent, signature-based metagenomic taxonomic profiler designed to reduce false discoveries while remaining practical to run on a workstation or laptop. It maps reads to precomputed unique signature fragments and estimates abundance from signature coverage and depth rather than marker genes.

> GOTTCHA v1 databases are not compatible with GOTTCHA2.

---

## Table of contents

- [What's new](#whats-new)
- [Installation](#installation)
- [Dependencies](#dependencies)
- [Databases](#databases)
- [Quick start](#quick-start)
- [Command overview](#command-overview)
- [Profiling](#profiling)
- [Fast profile mode](#fast-profile-mode)
- [Read extraction](#read-extraction)
- [Output files](#output-files)
- [Thresholds and filtering](#thresholds-and-filtering)
- [Full report fields](#full-report-fields)
- [Troubleshooting](#troubleshooting)
- [License and citation](#license-and-citation)

---

## What's new

Recent GOTTCHA2 releases through v2.4.0 include several workflow changes that are worth knowing before you start:

- **Fast prefiltering mode**: `fast-profile` uses `sylph` to prefilter the reference set before read mapping, which can substantially reducing both runtime and memory usage.
- **Current CLI**: the supported entry points are `profile`, `fast-profile`, `extract`, `sam2bam`, `download`, and `version`.
- **Updated identity handling**: the reported `SNI_SCORE` is based on consensus identity rather than the legacy read-weighted identity metric.
- **BAM-based workflow**: runs use sorted and indexed BAM for downstream processing instead of keeping SAM as the main intermediate.
- **Legacy compatibility**: the older `gottcha2.py` workflow (SAM-based) is still available for compatibility, but it is frozen at v2.2.3.

---

## Installation

### Option A: Conda

```bash
conda install -c bioconda gottcha2
```

### Option B: Install from source

Install the external tools first, then install the Python package:

```bash
# required for profile
# add sylph as well if you plan to use fast-profile
git clone https://github.com/poeli/GOTTCHA2
cd GOTTCHA2
python -m pip install .

# development install
python -m pip install -e .
```

Confirm the installation:

```bash
gottcha2 version
gottcha2 profile --help
gottcha2 fast-profile --help
```

For containerized usage, see [DOCKER.md](../DOCKER.md).

---

## Dependencies

GOTTCHA2 requires Python 3.9+.

Runtime dependencies:

- `minimap2` for mapping
- `samtools` and `pysam` for BAM conversion and parsing
- `numpy` and `pandas`
- `requests`
- `tqdm`
- `biom-format` if you use `--format biom`
- `sylph` if you use `fast-profile`

A Conda environment file is provided as `environment.yml`.

---

## Databases

### Prebuilt databases

The default download target used by `gottcha2 download` is:

```text
https://ref-db.edgebioinformatics.org/gottcha2/latest/gottcha_db.species.tar
```

You can also download database bundles manually from the same host if you prefer.

### Database bundle contents

A standard profiling database should include these files with the same prefix:

- `gottcha_db.<level>.fna.mmi` for `profile`
- `gottcha_db.<level>.fna.tax.tsv` taxonomy mapping
- `gottcha_db.<level>.fna.stats` signature and genome statistics

Additional files used by fast mode:

- `gottcha_db.<level>.fna.syldb` `sylph` database for prefiltering
- `gottcha_db.<level>.fna.zip` archived signature sequences used to build the reduced reference

You should pass the shared prefix or database directory to `-d/--database`. The required files will be automatically located. For example:

```text
/path/to/db/gottcha_db.species.fna
```

or 

```text
/path/to/db
```

### Download helper

[Not available yet] Use the built-in downloader to fetch the default database tarball into a new `database/` directory:

```bash
gottcha2 download
```

See available options with:

```bash
gottcha2 download --help
```

---

## Quick start

### 1) Profile Illumina paired-end reads

```bash
gottcha2 profile \
  -d /path/to/db/gottcha_db.species.fna \
  -i sample_R1.fastq.gz sample_R2.fastq.gz \
  -t 8 \
  -o out \
  -p sample
```

### 2) Profile Illumina single-end reads

```bash
gottcha2 profile \
  -d /path/to/db/gottcha_db.species.fna \
  -i sample.fastq.gz \
  -t 8 \
  -o out
```

### 3) Profile Oxford Nanopore reads

Nanopore mode expects a single input file:

```bash
gottcha2 profile \
  -d /path/to/db/gottcha_db.species.fna \
  -i ont_reads.fastq.gz \
  --nanopore \
  -t 8 \
  -o out
```

### 4) Re-run profiling from an existing BAM

If you already have a sorted and indexed BAM from a previous GOTTCHA2 run, you can re-aggregate with different cutoffs without remapping:

```bash
gottcha2 profile \
  -b sample.gottcha_species.bam \
  -d /path/to/db/gottcha_db.species.fna \
  -Mc 0.01 \
  -Mr 10 \
  -mi 0.95 \
  -t 8 \
  -o out
```

### 5) Run the faster prefiltering workflow

```bash
gottcha2 fast-profile \
  -d /path/to/db/gottcha_db.species.fna \
  -i sample.fastq.gz \
  -t 8 \
  -o out
```

### 6) Extract reads for a taxon from an existing BAM

```bash
gottcha2 extract \
  -b sample.gottcha_species.bam \
  -e 562
```

---

## Command overview

GOTTCHA2 uses a subcommand-style CLI:

```text
gottcha2 <command> [options]
```

Commands:

- `profile` - map reads and generate taxonomic profiles
- `fast-profile` - prefilter the database with `sylph`, then run `profile`
- `extract` - extract reads for one or more taxa from an existing BAM
- `sam2bam` - convert legacy GOTTCHA2 SAM output to sorted and indexed BAM
- `download` - download the default database bundle
- `version` - print the installed version

---

## Profiling

### Key concepts

Metagenomic samples are profiled by mapping sequencing reads directly to taxon-specific signature fragments. GOTTCHA2 consolidates alignments across each genome's signature space to compute coverage and depth statistics, then derives an ANI-like metric called the signature nucleotide identity score (`SNI_SCORE`). Genome-level results are subsequently aggregated to higher taxonomic ranks.

### Oxford Nanopore mode

With `--nanopore`, GOTTCHA2 first converts the input read file into a temporary FASTA of fixed-length chunks to make long reads easier to map. The current implementation:

- requires exactly one input file
- splits reads into non-overlapping 150 bp chunks
- drops the trailing remainder if it is shorter than 150 bp
- removes inconsistent chunk assignments after mapping

If you do not override them explicitly, Nanopore mode uses these defaults:

- `--matchIdentity 0.85`
- `--matchFraction 0.85`
- `--matchLength 100`
- `--errorRate 0.03`

### Accessions of interest

Use `--accList` to provide a text file containing one accession or signature ID per line. This is useful for plasmids, spike-ins, or other targets you want to track during profiling.

Use `--accListAction` to control how those reads are handled:

- `report_only` keeps all reads and reports the count in `AOI_READ_COUNT`
- `filter_out` removes reads matching listed accessions
- `filter_in` keeps only reads matching listed accessions

### Reporting level and database level

The database level is usually auto-detected from the database prefix or BAM name. For example, `gottcha_db.species.fna` implies `species`, and `sample.gottcha_species.bam` implies `species`.

If auto-detection is not possible, set it explicitly with `-l/--dbLevel`.

---

## Fast profile mode

`fast-profile` is a convenience wrapper for `profile --fast`. It adds a prefiltering step utilized by `sylph` before read mapping:

1. query the `.syldb` database against the input sample
2. collect the subset of candidate signatures
3. extract those signatures from the `.zip` archive
4. map reads only against that reduced reference

This mode is ideal when you need faster execution with a minimal memory footprint. It still produces the standard GOTTCHA2 outputs, including the BAM and summary reports.

---

## Read extraction

GOTTCHA2 can extract reads for one or more taxa from an existing BAM file. Taxa can be provided as:

- comma-separated taxids, for example `-e "666,562"`
- comma-separated taxon names, for example `-e "Vibrio cholerae,Escherichia coli"`
- a file prefixed with `@`, for example `-e "@taxids.txt"`

The `extract` command is shorthand for running `profile` with `--extract` and `--extractOnly`.

### Example usages

Extract reads mapping to taxid `666`:

```bash
gottcha2 extract \
  -b sample.gottcha_species.bam \
  -e 666
```

Extract with explicit match thresholds:

```bash
gottcha2 extract \
  -b sample.gottcha_species.bam \
  -e 666 \
  -mi 0.9 \
  -mf 0.9
```

Extract multiple taxa:

```bash
gottcha2 extract -b sample.gottcha_species.bam -e "1234,5678"
gottcha2 extract -b sample.gottcha_species.bam -e "@taxids.txt"
```

Limit the number of reads per taxon and choose the output format with `:N:FORMAT`:

```bash
# up to 1000 reads per taxon, FASTQ output
gottcha2 extract -b sample.gottcha_species.bam -e "@taxids.txt:1000:fastq"
```

Extract up to 20 representative sequences per profiled reference:

```bash
gottcha2 extract -b sample.gottcha_species.bam -ef
```

### Extracted record format

Each extracted FASTA or FASTQ header encodes the matched reference, interval, taxon, and match statistics:

```text
>{READ_NAME}{MATE}|{REFERENCE}:{START}..{END} LEVEL={LEVEL} NAME={NAME} TAXID={TAXID} AOI={AOI} MG={MG} MI={MI} MF={MF}
```

Field definitions:

- `READ_NAME`: read identifier
- `MATE`: paired-end suffix (`.1`, `.2`, or empty)
- `REFERENCE`: matched reference sequence name
- `START..END`: mapped reference interval (1-based)
- `LEVEL`: extracted taxonomic rank
- `NAME`: extracted taxon name
- `TAXID`: extracted taxonomy ID
- `AOI`: accession-of-interest flag
- `MG`: alignment length
- `MI`: mapping identity
- `MF`: mapping fraction

Example:

```text
>read123.1|chrA|1|300|GCF10000:10..120 LEVEL=species NAME=Escherichia_coli TAXID=562 AOI=False MG=148 MI=98.65 MF=0.99
ACGT...
```

---

## Output files

By default, outputs go to `--outdir` and use `--prefix` derived from the first input filename or the BAM name.

Typical outputs:

- `*.tsv`, `*.csv`, or `*.biom` - summary report at the requested reporting level
- `*.full.tsv` - full report including filtered taxa and notes
- `*.lineage.tsv` - lineage table for qualified taxa
- `*.mpa.tsv` - MetaPhlAn-style output when `--mpa` is enabled
- `*.extract.fasta` or `*.extract.fastq` - extracted reads when `--extract` or `extract` is used
- `*.gottcha_<level>.bam` and `*.bai` - sorted BAM and index for reuse
- `*.gottcha_<level>.log` - run log including thresholds and processing steps

---

## Thresholds and filtering

Most taxonomic cutoffs default to `0` and are disabled unless you set them explicitly. Alignment thresholds are always applied unless you lower them yourself.

Use `--noCutoff` to disable taxonomic profiling cutoffs. This is equivalent to:

```text
-Mc 0 -Mr 0 -Ml 0 -Mz 0 -ss 0,0,0
```

### Alignment thresholds

- `-mi, --matchIdentity <FLOAT>`
  Minimum alignment identity for a valid match. Default: `0.95` for short reads, `0.85` for Nanopore mode.

- `-mf, --matchFraction <FLOAT>`
  Minimum aligned fraction for a valid match. Default: `0.95` for short reads, `0.85` for Nanopore mode.

- `-mg, --matchLength <INT>`
  Minimum alignment length in bp. Default: `100`.

- `-er, --errorRate <FLOAT>`
  Estimated sequencing error rate. Default: `0.005` for short reads, `0.03` for Nanopore mode.

### Taxonomic profiling cutoffs

- `-ss, --sniScore <FLOAT>[,<FLOAT>,<FLOAT>]`
  SNI-score thresholds for `other,species,strain`. Default: `0.9,0.95,0.99`.

- `-Mc, --minCov <FLOAT>`
  Minimum signature coverage required for abundance calculation. Default: `0`.

- `-Mr, --minReads <INT>`
  Minimum number of mapped reads. Default: `0`.

- `-Ml, --minLen <INT>`
  Minimum covered signature length. Default: `0`.

- `-Mz, --maxZscore <FLOAT>`
  Maximum z-score for mapped-region depth distribution. Default: `0` (disabled).

Filtered taxa remain visible in `*.full.tsv`, with the reason recorded in `NOTE`.

---

## Full report fields

The full report (`<prefix>.full.tsv`) contains all computed metrics. The summary report contains the qualified rows shown at the requested reporting level.

| Field Name             | Description |
| ---------------------- | ----------- |
| LEVEL                  | Taxonomic rank (`superkingdom` through `strain`) |
| NAME                   | Taxon name |
| TAXID                  | NCBI taxonomy ID |
| READ_COUNT             | Reads mapped to this taxon |
| TOTAL_BP_MAPPED        | Total mapped bases across this taxon's signatures |
| SNI_SCORE              | Signature nucleotide identity used during filtering and aggregation |
| COVERED_SIG_LEN        | Total covered signature length |
| BEST_SIG_COV           | Highest signature coverage among rolled-up members |
| DEPTH                  | Depth of coverage (`TOTAL_BP_MAPPED / TOTAL_SIG_LEN`) |
| REL_ABUNDANCE_GC       | Relative abundance from genomic-content estimate |
| REL_ABUNDANCE          | Relative abundance from the field selected by `--relAbu` |
| PARENT_NAME            | Parent taxon name |
| PARENT_TAXID           | Parent taxonomy ID |
| AOI_READ_COUNT         | Reads matched to `--accList` entries |
| TOTAL_READ_LEN         | Total aligned read length |
| TOTAL_BP_MISMATCH      | Total mismatched bases |
| TOTAL_BP_INDEL         | Total inserted and deleted bases |
| READ_WT_SNI            | Read-weighted identity estimate |
| CONSENSUS_SEQ_SNI      | Consensus-sequence identity estimate |
| SNI_CI95_LH            | Low and high 95% confidence bounds for identity |
| SIG_COV                | Signature coverage (`COVERED_SIG_LEN / TOTAL_SIG_LEN`) |
| MAPPED_SIG_LEN         | Signature length with at least one mapped read |
| TOTAL_SIG_LEN          | Total signature length for the taxon |
| COVERED_SIG_DEPTH      | Depth across covered signature only |
| COVERED_MAPPED_SIG_COV | Covered fraction of mapped signature |
| ZSCORE                 | Depth-distribution z-score |
| GENOMIC_CONTENT_EST    | Genomic-content estimate |
| ABUNDANCE              | Raw abundance value from `--relAbu` |
| REL_ABUNDANCE_DEPTH    | Relative abundance computed from depth |
| SIG_LEVEL              | Signature rank used for mapping |
| GENOME_COUNT           | Number of rolled-up genomes |
| GENOME_SIZE            | Combined genome size used for GC normalization |
| NOTE                   | Filtering or rollup note |

---

## Troubleshooting

### BAM input must be sorted and indexed

If you provide `-b/--bam`, the BAM must already be coordinate-sorted and indexed.

For legacy GOTTCHA2 SAM output, convert it with:

```bash
gottcha2 sam2bam -i sample.sam -o sample.bam -t 8
```

### Database sidecar files are required

For `profile`, keep the database sidecar files next to the database prefix. At minimum, GOTTCHA2 expects:

```text
<db>.mmi
<db>.tax.tsv
<db>.stats
```

For `fast-profile`, it additionally expects:

```text
<db>.syldb
<db>.zip
<db>.tax.tsv
<db>.stats
```

### `--nanopore` requires one input file

Nanopore mode only accepts a single FASTA or FASTQ input file. If you have multiple files, merge them first or process them separately.

### Python and external dependency checks happen at runtime

GOTTCHA2 checks for:

- Python 3.9+
- `minimap2`
- `samtools`
- `sylph` when `fast-profile` is used

If one of these tools is missing from `PATH`, the run will stop before mapping begins.

### Identity and SNI changed from older releases

Modern GOTTCHA2 releases report `SNI_SCORE` from consensus identity. If you compare output against older `gottcha2.py` runs, expect differences in SNI-related columns and filtering behavior.

---

## License and citation

- License: [TBD]
- If you use GOTTCHA2 in publications, cite the GOTTCHA or GOTTCHA2 project, the database source, and the exact software version reported by `gottcha2 version`.
