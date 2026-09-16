[![logo](https://gottcha2.readthedocs.io/en/latest/_images/gottcha_icon.png)](https://gottcha2.readthedocs.io/en/latest/_images/gottcha_icon.png)

# Genomic Origin Through Taxonomic CHAllenge (GOTTCHA2)

GOTTCHA2 is a gene-independent, signature-based metagenomic taxonomic profiler for sequencing reads. It is designed to reduce false discoveries while remaining practical to run on a workstation or laptop. Instead of relying on marker genes, it maps reads to precomputed unique signature fragments and estimates abundance from signature coverage and depth.

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
- [Running unit tests](#running-unit-tests)
- [Troubleshooting](#troubleshooting)
- [License and citation](#license-and-citation)

---

## What's new

The current development version is v2.5.0. It includes several workflow changes that are worth knowing before you start:

- **Direct Oxford Nanopore profiling**: `-np/--nanopore` maps intact ONT reads by default. The earlier 150 bp chunk workflow remains available with `--ont-chunk`.
- **Shared mapping controls**: `--secondary yes|no`, `--max-secondary`, `--secondary-ratio`, and `--m2-options` apply to short reads and both ONT workflows. Secondary candidates are disabled by default.
- **Updated identity and coverage reporting**: `SNI_SCORE` uses consensus differences, sequencing error, and covered signature space. Species with both species-level and strain-level signatures combine their evidence before SNI estimation. Reports distinguish `ALN_IDENTITY`, `CONSENSUS_SEQ_IDENTITY`, `SIG_COV`, and `SIG_COV_RAW`.
- **Experimental reciprocal groups**: `--reciprocal-groups yes` enables grouping based on reciprocal species mappings. It is disabled by default and does not enable secondary candidates automatically.

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

- `minimap2` 2.27 or newer for mapping
- `samtools` and `pysam` for BAM conversion and parsing
- `numpy`, `pandas`, and `scipy`
- `requests`
- `tqdm`
- `biom-format` if you use `--format biom`
- `sylph` if you use `fast-profile`

A Conda environment file is provided as `environment.yml`.

---

## Databases

### Prebuilt databases

The database bundles can be downloaded using `gottcha2 download`. The default download target is:

```text
https://ref-db.edgebioinformatics.org/gottcha2/RefSeq-GTDB-multi-domain+/gottcha_db_fast.tar
```

This is the 32 GB prebuilt GOTTCHA2 database for the `gottcha2 fast-profile` mode. If you prefer, you can also download the database bundles manually from the same host. For the standard profiling mode, the 198 GB prebuilt database [`gottcha_db_standard.tar`](https://ref-db.edgebioinformatics.org/gottcha2/RefSeq-GTDB-multi-domain%2B/gottcha_db_standard.tar), is also available for download.

### Database bundle contents

A standard profiling database should include these files with the same prefix:

- `gottcha_db.<level>.fna.mmi` for `profile`
- `gottcha_db.<level>.fna.tax.tsv` taxonomy mapping
- `gottcha_db.<level>.fna.stats` signature and genome statistics

Fast mode uses the shared `.tax.tsv` and `.stats` files above together with:

- `gottcha_db.<level>.fna.syldb` `sylph` database for prefiltering
- `gottcha_db.<level>.fna.zip` archived signature sequences used to build the reduced reference

The `.mmi` index is required by `profile` but is not used by `fast-profile`.

Pass either the shared database prefix or the database directory to `-d/--database`. GOTTCHA2 locates the required sidecar files from that path. For example:

```text
/path/to/db/gottcha_db.species.fna
```

or

```text
/path/to/db
```

### Download helper

Use the built-in downloader to fetch and verify a database tarball, then extract it into a new `database/` directory. Fast mode is the default download:

```bash
gottcha2 download
```

Download the larger standard database instead with:

```bash
gottcha2 download -d standard
```

The downloader stops if a `database/` directory already exists in the current working directory. It verifies the archive against the published SHA-256 checksum before extraction.

See available options with:

```bash
gottcha2 download --help
```

---

## Quick start

These examples use the most common workflows. Replace the example paths and filenames with your own sample and database locations.

### Example conventions

- `gottcha2 profile` maps reads to the selected signature database and produces taxonomic reports.
- `gottcha2 fast-profile` first narrows the reference set, then runs the profiling workflow on the reduced reference.
- `gottcha2 extract` pulls reads assigned to selected taxa from an existing GOTTCHA2 BAM file.
- `-d/--database` points to either a database prefix, such as `/path/to/db/gottcha_db.species.fna`, or to a directory that contains the matching database files.
- `-i/--input` supplies one or more read files. Use two files for paired-end Illumina reads and one file for single-end or Nanopore reads.
- `-b/--bam` reuses an existing sorted and indexed BAM instead of remapping reads.
- `-np/--nanopore` enables ONT mode and maps intact long reads by default.
- `--ont-chunk`, used together with `-np`, selects the earlier 150 bp chunk workflow.
- `-t/--threads` controls the number of CPU threads used by mapping and related processing steps.
- `-o/--outdir` chooses the output directory. GOTTCHA2 creates it if needed.
- `-p/--prefix` sets the output filename prefix. If you omit it, GOTTCHA2 derives a prefix from the input filename or BAM name.
- The backslash (`\`) at the end of a line lets long shell commands continue on the next line. You can also write each example as a single line.

### 1) Profile Illumina paired-end reads

Use this when your sample has forward and reverse FASTQ files.

```bash
gottcha2 profile \
  -d /path/to/db/gottcha_db.species.fna \
  -i sample_R1.fastq.gz sample_R2.fastq.gz \
  -t 8 \
  -o out \
  -p sample
```

What this command does:

- loads the species-level database specified with `-d`
- maps both paired-end read files supplied after `-i`
- uses 8 threads because of `-t 8`
- writes results into the `out/` directory
- names output files with the `sample` prefix, for example `sample.tsv`, `sample.full.tsv`, and `sample.gottcha_species.bam`

Use a sample-specific prefix whenever you process multiple samples into the same output directory.

### 2) Profile Illumina single-end reads

Use this when each sample has one FASTQ file.

```bash
gottcha2 profile \
  -d /path/to/db/gottcha_db.species.fna \
  -i sample.fastq.gz \
  -t 8 \
  -o out
```

Because `-p/--prefix` is omitted, GOTTCHA2 derives the output prefix from `sample.fastq.gz`. Add `-p sample_name` if you want a shorter or more explicit prefix.

### 3) Profile Oxford Nanopore reads

Nanopore mode expects exactly one input file. Add `-np` (short for `--nanopore`) so GOTTCHA2 maps the intact long reads with ONT-oriented settings.

```bash
gottcha2 profile \
  -d /path/to/db/gottcha_db.species.fna \
  -i ont_reads.fastq.gz \
  -np \
  -t 8 \
  -o out \
  -p ont_sample
```

By default, `-np` uses direct ONT mode with the `lr:hq` minimap2 preset and secondary candidates disabled. To use the earlier chunk-based workflow instead, add `--ont-chunk`:

```bash
gottcha2 profile \
  -d /path/to/db/gottcha_db.species.fna \
  -i ont_reads.fastq.gz \
  -np --ont-chunk \
  -t 8 \
  -o out \
  -p ont_sample.chunk
```

See [Oxford Nanopore mode](#oxford-nanopore-mode) for the different defaults and tuning options.

### 4) Re-run profiling from an existing BAM

Use this when you already have a sorted and indexed GOTTCHA2 BAM and want to re-aggregate results with different thresholds. This avoids the slower read-mapping step.

```bash
gottcha2 profile \
  -b sample.gottcha_species.bam \
  -d /path/to/db/gottcha_db.species.fna \
  -Mc 0.01 \
  -Mr 10 \
  -mi 0.95 \
  -t 8 \
  -o out \
  -p sample.refiltered
```

What the non-default options mean:

- `-b sample.gottcha_species.bam` reads alignments from an existing BAM instead of using `-i` input reads.
- `-Mc 0.01` requires at least 1% signature coverage for abundance calculation.
- `-Mr 10` requires at least 10 mapped reads.
- `-mi 0.95` keeps only matches with at least 95% alignment identity.
- `-p sample.refiltered` keeps the re-filtered output separate from the original run.

The BAM must be coordinate-sorted and indexed. Keep the database path consistent with the database used for the original mapping.

### 5) Run the faster prefiltering workflow

Use `fast-profile` when you want to reduce runtime and memory usage while producing results comparable to the standard `profile` workflow. It does this by preselecting likely reference sequences before read mapping.

```bash
gottcha2 fast-profile \
  -d /path/to/db/gottcha_db.species.fna \
  -i sample.fastq.gz \
  -t 8 \
  -o out \
  -p sample.fast
```

This command requires the fast-mode `.syldb` and `.zip` files plus the shared `.tax.tsv` and `.stats` files. It does not use the standard `.mmi` index. The outputs have the same general structure as `profile`, but mapping is performed against a reduced reference set selected by `sylph`.

### 6) Extract reads for a taxon from an existing BAM

Use `extract` after profiling when you want the reads assigned to one or more taxa. The example below extracts reads assigned to NCBI taxid `562`.

```bash
gottcha2 extract \
  -b sample.gottcha_species.bam \
  -e 562 \
  -o out \
  -p sample.ecoli
```

Key options:

- `-b` points to the GOTTCHA2 BAM created by `profile` or `fast-profile`.
- `-e` selects the taxon or taxa to extract. You can use taxids, names, or `@file` syntax.
- `-o` and `-p` control where the extracted FASTA or FASTQ output is written.

### Check the results

After any profiling run, start with these files:

```bash
ls out
less out/sample.gottcha_species.log
column -t -s $'\t' out/sample.tsv | less -S
column -t -s $'\t' out/sample.full.tsv | less -S
```

The TSV/CSV summary contains qualified rows across taxonomic ranks. The BIOM report contains the database rank selected with `--dbLevel`. TSV/CSV runs also write a full report (`*.full.tsv`) containing passing, filtered, and hidden rows, with reasons in the `NOTE` column.

---

## Command overview

GOTTCHA2 uses a subcommand-style CLI:

```text
gottcha2 <command> [options]
```

| Command | Use it when you need to... | Typical starting point |
| ------- | -------------------------- | ---------------------- |
| `profile` | Map reads or reuse a BAM and generate taxonomic profiles. | `gottcha2 profile -d DB -i reads.fastq.gz -o out` |
| `fast-profile` | Prefilter the database with `sylph`, then run profiling on a reduced reference set. | `gottcha2 fast-profile -d DB -i reads.fastq.gz -o out` |
| `extract` | Extract reads assigned to one or more taxa from an existing BAM. | `gottcha2 extract -b sample.bam -e 562` |
| `sam2bam` | Convert legacy GOTTCHA2 SAM output into sorted, indexed BAM. | `gottcha2 sam2bam -i sample.sam -o sample.bam` |
| `download` | Download the default database bundle, when supported by your build. | `gottcha2 download` |
| `version` | Print the installed GOTTCHA2 version. | `gottcha2 version` |

Use `--help` after any command to see command-specific options, defaults, and examples:

```bash
gottcha2 profile --help
gottcha2 extract --help
```

---

## Profiling

### Key concepts

GOTTCHA2 profiles metagenomic samples by mapping sequencing reads directly to taxon-specific signature fragments. GOTTCHA2 consolidates alignments across each genome's signature space to compute coverage and depth statistics, then derives an ANI-like metric called the signature nucleotide identity score (`SNI_SCORE`). Genome-level results are subsequently aggregated to higher taxonomic ranks.

`SNI_SCORE` is the center of a coverage-adjusted Wilson interval, using consensus differences after subtracting the estimated sequencing error rate. It is not identical to raw alignment identity. Sparse signature coverage widens `SNI_CI95_LH` and affects the score.

At the species rank, evidence from species-level and strain-level signatures is combined before estimating SNI. If only strain-level signatures have mapped reads, the species calculation also includes the representative species signature length from the database statistics. Other rollups retain the best constituent SNI score and its confidence interval. Filtering is evaluated at each rank, so a filtered strain can still contribute evidence to a qualifying species or genus.

### Oxford Nanopore mode

Use `-np/--nanopore` for a single ONT FASTA or FASTQ file. In v2.5.0 development, this selects direct mode by default. Direct mode maps each intact read to the signature database using `lr:hq` and the thresholds below.

The current profiling workflow does not run the earlier per-read species-support resolver. Mapping postprocessing selects the best primary alignment per read or mate before BAM conversion.

The earlier chunk workflow remains available with `-np --ont-chunk`. It splits reads into non-overlapping 150 bp pieces, drops a trailing piece shorter than 150 bp, and maps the pieces as short reads. The current profiling workflow does not run the earlier taxonomic consistency filter for chunks.

| Setting | Direct mode: `-np` | Chunk mode: `-np --ont-chunk` |
| ------- | ------------------ | ----------------------------- |
| Input used for mapping | Intact ONT reads | Non-overlapping 150 bp chunks |
| minimap2 preset | `lr:hq` | `sr` |
| `--matchIdentity` | `0.85` | `0.85` |
| `--matchFraction` | `0` | `0.85` |
| `--matchLength` | `100` bp | `100` bp |
| `--errorRate` | `0.01` | `0.03` |
| `--secondary` | `no` | `no` |
| Automatic `--m2-options` | `-n1 -m25 -s120 --no-long-join` | `-s120` |

The lower direct-mode match fraction is intentional: a short signature alignment can cover only a small fraction of an intact long read. An alignment passes this threshold when the aligned span covers the required fraction of either the read or the signature fragment.

### Mapping controls

These options apply to short reads and both Nanopore workflows:

- `--secondary yes|no`: request secondary candidates from minimap2 and allow secondary/supplementary alignments during BAM parsing. Default: `no`. Mapping postprocessing can still remove these candidates before BAM conversion.
- `--max-secondary <INT>`: maximum secondary alignments requested per primary alignment when `--secondary yes` is used. Default: `10`; must be non-negative.
- `--secondary-ratio <FLOAT>`: minimum secondary-to-primary minimap2 chaining-score ratio when secondary candidates are enabled. Default: `0.9`; accepted range: `0` to `1`.
- `-xm/--presetx <STR>`: override the minimap2 preset. Direct mode defaults to `lr:hq`; other accepted values are `sr`, `map-pb`, and `map-ont`.
- `--m2-options <STR>`: replace the automatically selected minimap2 tuning options. For values beginning with `-`, use the equals form, for example `--m2-options="-n2 -m25 -s150 --no-long-join"`.

For example, this command requests up to 25 secondary candidates with a minimum chaining-score ratio of 0.7:

```bash
gottcha2 profile \
  -d /path/to/db/gottcha_db.species.fna \
  -i ont_reads.fastq.gz \
  -np \
  --secondary yes \
  --max-secondary 25 \
  --secondary-ratio 0.7 \
  -t 8 \
  -o out
```

`--ont-chunk` itself is only valid together with `-np/--nanopore`. For direct ONT mapping in `fast-profile`, the automatic options use `-n2` instead of `-n1` because mapping uses the extracted reference.

### Experimental reciprocal groups

Enable reciprocal grouping with `--reciprocal-groups yes`. The grouping algorithm links species in the same genus only when mappings support both directions. By default, each direction needs at least three independent read names and a 95% Wilson lower bound of at least `0.005` for the directional association fraction. Connected species form a group, including species connected through intermediate members.

During aggregation, the observed species with the largest summed value of `--relAbu` represents its group. The default abundance field is `DEPTH`. Other members contribute to that species and its higher ranks; their strain identifiers are preserved, and `NOTE` records `Grouped with species ...`.

```bash
gottcha2 profile \
  -d /path/to/db/gottcha_db.species.fna \
  -i sample.fastq.gz \
  --secondary yes \
  --reciprocal-groups yes \
  -o out \
  -p sample.grouped
```

This option uses the intermediate SAM file available during mapping. Reusing a BAM alone does not reconstruct reciprocal groups. Mapping postprocessing removes competing alignments before the grouping step, so enabling this experimental option may leave species separate. Secondary candidate generation and reciprocal grouping are independent options, both disabled by default.

### Signature of interest

Use `--sigList` to provide a text file containing one accession or signature ID per line. This is useful for plasmids, spike-ins, or other targets you want to track during profiling.

Use `--sigListAction` to control how those reads are handled:

- `report_only` keeps all reads and reports the count in `AOI_READ_COUNT`
- `filter_out` removes reads matching listed accessions
- `filter_in` keeps only reads matching listed accessions

### Reporting level and database level

The database level is usually auto-detected from the database prefix or BAM name. For example, `gottcha_db.species.fna` implies `species`, and `sample.gottcha_species.bam` implies `species`.

If auto-detection is not possible, set it explicitly with `-l/--dbLevel`.

---

## Fast profile mode

`fast-profile` is a convenience wrapper for `profile --fast`. It adds a `sylph` prefiltering step before read mapping:

1. Query the `.syldb` database against the input sample.
2. Collect the subset of candidate signatures.
3. Extract those signatures from the `.zip` archive.
4. Map reads only against that reduced reference.

This mode is useful when you need faster execution with a smaller memory footprint. It still produces the standard GOTTCHA2 outputs, including the BAM and summary reports.

Nanopore selection works the same way in fast mode: `fast-profile -np` maps intact ONT reads by default, while `fast-profile -np --ont-chunk` uses the chunk workflow. Direct fast mode uses `-n2 -m25 -s120 --no-long-join` by default for the reduced reference extracted by the prefilter.

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

By default, outputs go to `--outdir` and use a prefix derived from `--prefix`, the first input filename, or the BAM name.

Typical outputs:

- `*.tsv` or `*.csv` - summary report with qualified rows across taxonomic ranks
- `*.biom` - BIOM report at the database rank selected with `--dbLevel`
- `*.full.tsv` - full report including filtered taxa and notes, written for TSV/CSV runs; with `--format csv`, this file contains comma-separated data despite its `.tsv` suffix
- `*.lineage.tsv` - lineage table for qualified taxa
- `*.mpa.tsv` - MetaPhlAn-style output when `--mpa` is enabled
- `*.extract.fasta` or `*.extract.fastq` - extracted reads when `--extract` or `extract` is used
- `*.gottcha_<level>.bam` and `*.bai` - sorted BAM and index for reuse
- `*.gottcha_<level>.log` - run log including thresholds and processing steps

---

## Thresholds and filtering

Coverage, read-count, covered-length, and z-score cutoffs default to `0` and are disabled unless you set them explicitly. SNI and alignment thresholds are enabled by default.

`--noCutoff` sets SNI thresholds to `0,0,0`. With the other cutoffs left at their defaults, this disables profiling-stage cutoffs. Explicit `-Mc`, `-Mr`, `-Ml`, or `-Mz` values remain active. To disable all profiling cutoffs explicitly, use:

```text
-Mc 0 -Mr 0 -Ml 0 -Mz 0 -ss 0,0,0
```

### Alignment thresholds

- `-mi, --matchIdentity <FLOAT>`
  Minimum alignment identity for a valid match. Default: `0.95` for short reads and `0.85` for both Nanopore workflows.

- `-mf, --matchFraction <FLOAT>`
  Minimum aligned fraction of the read or signature fragment for a valid match. Default: `0.95` for short reads, `0` for direct Nanopore mode, and `0.85` for Nanopore chunk mode.

- `-mg, --matchLength <INT>`
  Minimum alignment length in bp. Default: `100`.

### Taxonomic profiling cutoffs

- `-er, --errorRate <FLOAT>`
  Estimated sequencing error rate used for SNI inference. Default: `0.005` for short reads, `0.01` for direct Nanopore mode, and `0.03` for Nanopore chunk mode.

- `-ss, --sniScore <FLOAT>[,<FLOAT>,<FLOAT>]`
  SNI-score thresholds for `other,species,strain`. Default: `0.9,0.95,0.99`. One value applies to all ranks; two values supply `other,species` and retain the strain default of `0.99`.

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

The full report (`<prefix>.full.tsv`) contains the fields below. The TSV/CSV summary contains qualified rows across ranks and the first 11 columns, from `LEVEL` through `REL_ABUNDANCE`. Rows whose `NOTE` contains `Filtered out` or `Not shown` remain only in the full report.

| Field Name             | Description |
| ---------------------- | ----------- |
| LEVEL                  | Taxonomic rank (`superkingdom` through `strain`) |
| NAME                   | Taxon name |
| TAXID                  | Taxonomy identifier from the selected database |
| READ_COUNT             | Read count accumulated from accepted reference alignments |
| TOTAL_BP_MAPPED        | Total mapped bases across this taxon's signatures |
| SNI_SCORE              | Coverage-adjusted, error-corrected consensus identity score used during filtering and aggregation |
| COVERED_SIG_LEN        | Total covered signature length |
| SIG_COV                | Coverage used for filtering: covered/total signature length for strains and recomputed species; otherwise the highest constituent coverage |
| DEPTH                  | Mapped bases / total signature length for each reference entry; summed across entries during rollup |
| REL_ABUNDANCE_GC       | Relative abundance from genomic-content estimate |
| REL_ABUNDANCE          | Relative abundance from the field selected by `--relAbu` |
| PARENT_NAME            | Parent taxon name |
| PARENT_TAXID           | Parent taxonomy ID |
| AOI_READ_COUNT         | Reads matched to `--sigList` entries |
| TOTAL_READ_LEN         | Total aligned query length across accepted alignments |
| TOTAL_BP_MISMATCH      | Total mismatched bases |
| TOTAL_BP_INDEL         | Total inserted and deleted bases |
| ALN_IDENTITY           | Alignment identity (`1 - TOTAL_BP_MISMATCH / TOTAL_BP_MAPPED`) |
| CONSENSUS_SEQ_IDENTITY  | Consensus identity before error correction (`1 - CONSENSUS_DIFF / COVERED_SIG_LEN`); `CONSENSUS_DIFF` is an internal count of covered positions with majority mismatches |
| SNI_CI95_LH            | Low and high 95% confidence bounds for SNI, formatted as `[low-high]` |
| SIG_COV_RAW            | Aggregate signature coverage (`COVERED_SIG_LEN / TOTAL_SIG_LEN`) |
| MAPPED_SIG_LEN         | Signature length with at least one mapped read |
| TOTAL_SIG_LEN          | Total signature length for the taxon |
| COVERED_SIG_DEPTH      | Depth across covered signature only |
| COVERED_MAPPED_SIG_COV | Covered fraction of mapped signature |
| ZSCORE                 | Depth-distribution z-score |
| GENOMIC_CONTENT_EST    | Genomic-content estimate |
| ABUNDANCE              | Raw abundance value from `--relAbu` |
| REL_ABUNDANCE_DEPTH    | Relative abundance computed from depth |
| SIG_LEVEL              | Signature rank used for mapping |
| GENOME_COUNT           | Number of reference entries aggregated at this rank |
| GENOME_SIZE            | Combined genome size, with a representative-species adjustment for species supported only by strain signatures |
| NOTE                   | Filtering or rollup note |

---

## Running unit tests

After installing the package and its Python dependencies, run the same unit-test discovery command used in CI from the repository root:

```bash
python -m unittest discover test -v
```

The unit suite covers CLI dispatch and defaults, ONT preprocessing, read extraction, reciprocal grouping, SNI and taxonomic aggregation, and TSV/CSV report fields. It uses small in-memory fixtures and temporary files, so it does not need a downloaded signature database or external mapping commands. The CI workflow separately runs functional profiling tests with the bundled Ebola test database.

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

For `fast-profile`, it expects:

```text
<db>.syldb
<db>.zip
<db>.tax.tsv
<db>.stats
```

### `-np/--nanopore` requires one input file

Nanopore mode only accepts a single FASTA or FASTQ input file. If you have multiple files, merge them first or process them separately.

### `--ont-chunk` requires Nanopore mode

The chunk workflow is selected with `-np --ont-chunk`. Using `--ont-chunk` without `-np/--nanopore` is rejected because chunk preprocessing and postprocessing are specific to ONT reads.

### Python and external dependency checks happen at runtime

GOTTCHA2 checks for:

- Python 3.9+
- `minimap2` 2.27 or newer
- `samtools`
- `sylph` when `fast-profile` is used

If one of these tools is missing from `PATH`, the run stops before mapping begins. Install the missing tool or activate the environment that contains it, then rerun the command.

### No taxa are reported

If the summary report is empty, check the full report and the log before rerunning:

- `*.full.tsv` shows filtered taxa and the reason in `NOTE`.
- `*.gottcha_<level>.log` records the thresholds and database files used.
- Lowering `-Mc`, `-Mr`, `-mi`, or `-mf` can increase sensitivity, but may also increase false positives.

### Fast profile cannot find `.syldb` or `.zip`

`fast-profile` requires `.syldb`, `.zip`, `.tax.tsv`, and `.stats` files with a shared prefix. It does not require the standard `.mmi` index. If the fast-mode files are missing, either download a fast-mode-compatible database bundle or use `profile` with the standard `.mmi` database.

### Identity and SNI changed from older releases

Modern GOTTCHA2 releases report `SNI_SCORE` from consensus identity. If you compare output against older `gottcha2.py` runs, expect differences in SNI-related columns and filtering behavior.

---

## License and citation

- License: BSD 3-Clause
- If you use GOTTCHA2 in publications, cite the GOTTCHA or GOTTCHA2 project, the database source, and the exact software version reported by `gottcha2 version`.
