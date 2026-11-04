# LEAP

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A523.04.0-23aa8c.svg)](https://www.nextflow.io/)
[![Python](https://img.shields.io/badge/python-%E2%89%A53.7-3776AB?logo=python&logoColor=white)](https://www.python.org/)
[![Run with conda](https://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![Run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

LEAP is a Nextflow DSL2 pipeline for extending current gene annotations with long-read RNA-seq evidence. It identifies transcript-end extensions from long-read splice structures, validates 5' ends with CAGE data and 3' ends with polyA data, then applies safety checks to avoid creating unsupported overlaps with existing annotation.

The pipeline is organised around chromosome- and strand-aware processing so large BAM and annotation datasets can be run in parallel on HPC infrastructure.

## Contents

- [What LEAP Does](#what-leap-does)
- [Requirements](#requirements)
- [Inputs](#inputs)
- [Quick Start](#quick-start)
- [Parameters](#parameters)
- [Workflow](#workflow)
- [Outputs](#outputs)
- [Project Structure](#project-structure)
- [Configuration](#configuration)
- [Development Notes](#development-notes)
- [License](#license)

## What LEAP Does

- Converts long-read BAM alignments into BED, GenePred, and GTF-derived transcript structures.
- Extracts first-exon and last-exon unique splice sites from reads and current annotation.
- Uses CAGE and polyA BED files as independent evidence for transcript-end extensions.
- Retrieves current annotation from an Ensembl-style MySQL/MariaDB database.
- Filters partial transcripts and separates multi-exon from single-exon models.
- Detects natural overlaps in existing annotation so legitimate overlaps can be preserved.
- Selects the longest supported extension per terminal exon.
- Checks and fixes unsafe overlaps before generating extended transcript GFF files.

## Requirements

- Nextflow 23.04.0 or newer
- Linux or HPC execution environment
- Singularity/Apptainer or Conda
- Python 3 with `pandas`, `numpy`, and `matplotlib` available to helper scripts
- Perl and MySQL client tools for annotation database access
- Access to an annotation database compatible with `bin/db2gxf.pl`

Several modules declare Biocontainers images for tools such as `samtools`, `bedtools`, `bedToGenePred`, and `genePredToGtf`. The default config enables both Singularity and Conda:

```groovy
singularity.enabled = true
conda.enabled = true
```

## Inputs

LEAP expects:

- A coordinate-sorted long-read RNA-seq BAM file.
- A BAM index file (`.bai` or `.csi`).
- CAGE peaks in BED-like format for 5' validation.
- polyA sites in BED-like format for 3' validation.
- Database connection details for the current annotation source.
- An output directory for published results.

The CAGE and polyA files are split by chromosome and strand using the chromosomes found in the read BAM.

## Quick Start

Edit [params.config](/Users/rmerritt/my_app/LEAP/params.config) or provide parameters on the command line.

```groovy
params {
    output_directory = "/path/to/leap_results"

    reads = "/path/to/long_reads.bam"
    reads_index = "/path/to/long_reads.bam.bai"
    cage = "/path/to/cage.bed"
    polya = "/path/to/polya.bed"

    minimum_read_support = 2
    overlap_penalty = 10
    extend_on_current_annotation = true

    host = "database.host"
    port = 3306
    user = "username"
    pass = "password"
    dbname = "annotation_database"
}
```

Run locally with the default configuration:

```bash
nextflow run main.nf
```

Run on the EBI Slurm profile:

```bash
nextflow run main.nf -profile hpc_slurm
```

Resume a previous run:

```bash
nextflow run main.nf -profile hpc_slurm -resume
```

You can also override parameters at run time:

```bash
nextflow run main.nf \
  --reads /path/to/long_reads.bam \
  --reads_index /path/to/long_reads.bam.bai \
  --cage /path/to/cage.bed \
  --polya /path/to/polya.bed \
  --output_directory /path/to/leap_results \
  --minimum_read_support 2 \
  --overlap_penalty 10 \
  --extend_on_current_annotation true \
  --host database.host \
  --port 3306 \
  --user username \
  --pass password \
  --dbname annotation_database \
  -profile hpc_slurm
```

## Parameters

| Parameter | Required | Description |
| --- | --- | --- |
| `reads` | Yes | Long-read RNA-seq BAM file. |
| `reads_index` | Yes | Index for `reads`; used when extracting chromosomes and splitting BAMs. |
| `cage` | Yes | CAGE validation data. Used for first-exon/5' extension support. |
| `polya` | Yes | polyA validation data. Used for last-exon/3' extension support. |
| `output_directory` | Yes | Directory where selected pipeline outputs are published. |
| `minimum_read_support` | Yes | Minimum read support passed to read unique-splice-site extraction. |
| `overlap_penalty` | Yes | Alpha/penalty value passed to unsafe-overlap fixing. |
| `extend_on_current_annotation` | Yes | Controls whether candidate read extensions may be matched against current annotation. |
| `host` | Yes | Annotation database host. |
| `port` | Yes | Annotation database port. |
| `user` | Yes | Annotation database user. |
| `pass` | Yes | Annotation database password. |
| `dbname` | Yes | Annotation database name. |
| `sample_sheet` | No | Present in `params.config`, but not currently consumed by `main.nf`. |

`params.config` contains empty placeholders. Set values before running, either in that file or with `--parameter value` overrides.

## Workflow

The top-level workflow in [main.nf](/Users/rmerritt/my_app/LEAP/main.nf) runs four subworkflows.

### 1. `READ_PROCESSING`

Defined in [subworkflows/local/read_processing/main.nf](/Users/rmerritt/my_app/LEAP/subworkflows/local/read_processing/main.nf).

1. Extract chromosomes from the input BAM.
2. Split reads by chromosome.
3. Convert BAM to BED, BED to GenePred, and GenePred to GTF.
4. Split read-derived GTF files by strand.
5. Extract first and last exons.
6. Extract unique splice sites from read terminal exons.

### 2. `PROCESS_CAGE_AND_POLYA`

Defined in [subworkflows/local/process_cage_and_polya/main.nf](/Users/rmerritt/my_app/LEAP/subworkflows/local/process_cage_and_polya/main.nf).

1. Split CAGE data by chromosome and strand.
2. Split polyA data by chromosome and strand.
3. Emit keyed channels for downstream joins with read splice-site evidence.

### 3. `PROCESS_ANNOTATION`

Defined in [subworkflows/local/process_annotation/main.nf](/Users/rmerritt/my_app/LEAP/subworkflows/local/process_annotation/main.nf).

1. Retrieve current annotation per chromosome from the configured database.
2. Remove readthrough genes and OTT IDs during annotation retrieval.
3. Build transcript lists that should not be extended.
4. Classify full-length and partial coding transcripts, retaining both for overlap processing.
5. Split the combined coding annotation into multi-exon and single-exon transcripts.
6. Split current annotation by strand.
7. Extract first and last exons from current annotation.
8. Identify natural overlaps in current annotation.
9. Extract unique splice sites from current annotation.

### 4. `OVERLAP_READS`

Defined in [subworkflows/local/overlap_reads/main.nf](/Users/rmerritt/my_app/LEAP/subworkflows/local/overlap_reads/main.nf).

1. Overlap read-derived splice-site evidence with CAGE and polyA data.
2. Match validated candidate extensions against current annotation.
3. Select the longest possible extension.
4. Add gene IDs to selected extensions.
5. Filter internally overlapping features.
6. Intersect first-exon and last-exon extensions.
7. Check extension safety against natural overlaps.
8. Fix unsafe overlaps.
9. Build final extended transcript GFF files and per-chromosome statistics, excluding partial transcript models.

Partial coding models remain in the intermediate annotation so their CDS, unique
exons, terminal exons and single-exon features participate in the existing LEAP
overlap and natural-overlap decisions. They also participate in annotation-based
grouping at shared splice sites. The overlap predicates, shorter-candidate
selection and annotation fallback are unchanged. Only full-length multi-exon
models are emitted in the final extended GFF and extension stats, matching the
previous output eligibility.

`filter_out_partials.py` retains its `_FT.gff`, `_PT.gff` and `_NPC.gff` outputs
and additionally writes `_PC.gff`, containing both full-length and partial
protein-coding models. Nextflow uses `_PC.gff` for the normal processing stages
and passes `_PT.gff` to the final writer. Direct calls to
`build_chromosome_gff_and_stats.py` now require `--partial-annotation <PT.gff>`;
an empty file is valid when no partials exist. Partial IDs are matched exactly,
including version, chromosome and strand. All records belonging to those models
are omitted before applying extensions, so they cannot contribute stats or final
gene boundaries. See [partial-transcript processing](docs/partial_transcript_processing.md)
for testing and replay limitations.

A more detailed workflow diagram is available in [BIOINFORMATICS_FLOWCHART.md](/Users/rmerritt/my_app/LEAP/BIOINFORMATICS_FLOWCHART.md).

## Outputs

Selected outputs are published under `params.output_directory`.

| Directory | Contents |
| --- | --- |
| `raw_annotation/` | Per-chromosome current annotation GFF3 files from the configured database. |
| `natural_overlaps/` | Existing annotation overlaps that are treated as natural/known overlaps. |
| `extension_safety/` | Safety-check outputs, including safe and unsafe extension BED files. |
| `unsafe_overlaps/` | Rejected unsafe overlaps, named like `<chromosome>_<strand>_rejected.bed`. |
| `fixed_extensions/` | First-exon and last-exon BED files after unsafe-overlap fixing. |
| `extended_transcripts/` | Final per-chromosome extended transcript GFF files. |
| `extended_transcripts_stats/` | Per-chromosome TSV statistics for the extended transcript output. |

Intermediate files remain in the Nextflow work directory unless a module publishes them.

## Project Structure

```text
LEAP/
├── main.nf
├── nextflow.config
├── params.config
├── conf/
│   └── ebi_slurm.config
├── bin/
│   ├── get_current_annotation.sh
│   ├── db2gxf.pl
│   └── *.py
├── assets/
│   └── *.py
├── modules/
│   └── local/
├── subworkflows/
│   └── local/
├── BIOINFORMATICS_FLOWCHART.md
├── LICENSE
└── README.md
```

## Configuration

The base [nextflow.config](/Users/rmerritt/my_app/LEAP/nextflow.config) includes [params.config](/Users/rmerritt/my_app/LEAP/params.config), enables Singularity and Conda, and defines the `hpc_slurm` profile.

The [conf/ebi_slurm.config](/Users/rmerritt/my_app/LEAP/conf/ebi_slurm.config) profile is configured for the EMBL-EBI Codon Slurm cluster. It sets the Slurm executor, queue behaviour, retry labels, and resource labels used by local modules.

## Development Notes

- The pipeline uses Nextflow DSL2 modules in `modules/local`.
- Reusable workflow units live in `subworkflows/local`.
- Helper scripts in `bin` are on the Nextflow process path and implement most of the annotation, overlap, and extension logic.
- Scripts in `assets` provide additional checking, filtering, and output-combination utilities.
- Keep module outputs keyed by chromosome, strand, and exon type when adding new steps; downstream joins depend on those keys.

## Transcript boundary datacheck

Run the standalone check on a raw annotation and its corresponding extended annotation:

```bash
python3 assets/check_for_shorter_transcripts.py current_annotation_1.gff3 1_positive_extended.gff comparison.csv
```

Both inputs accept plain GFF3 or `.gz` files. The check uses `transcript` features and exact `transcript_id` values, including version suffixes. Each matched transcript must retain its chromosome and strand and satisfy `extended_start <= raw_start` and `extended_end >= raw_end`, on either strand. Growth at one end does not compensate for shortening at the other. Exon structure is not checked.

The CSV includes all IDs from both inputs, both sets of chromosome, strand, start, end and inclusive genomic length (`end - start + 1`), and boundary-shortening and location-mismatch flags. Status is `unchanged`, `extended`, `shortened`, `location_mismatch`, `missing_from_extended` or `missing_from_raw`. Location mismatches take precedence over boundary comparisons; shortening flags are blank for these rows and for unmatched IDs. Missing annotation values are blank.

The console reports a pass/fail result and counts for every status. Missing transcripts are reported separately and do not fail the check, since raw annotations may contain transcripts filtered before extension. A chromosome-wide raw file can therefore be compared with a strand-specific extended file.

Exit codes are `0` for pass, `1` for shortened transcripts or chromosome/strand mismatches, and `2` for input or execution errors. Invalid transcript records, missing or duplicate IDs, invalid coordinates and files without transcript records are rejected. This script is not yet integrated into Nextflow.

Run the focused tests with `python3 -m unittest discover -s tests`.

## Internal exon structure datacheck

```bash
python3 assets/check_internal_exon_positions.py current_annotation_1.gff3 1_positive_extended.gff exon_comparison.csv
```

This standalone, standard-library check accepts plain or gzipped GFF3 and matches exact `transcript_id` values. It sorts exons by genomic coordinates, ignoring file order, exon IDs and exon numbers. Matched transcripts must retain chromosome, strand and exon count. Every exon start except the lowest-coordinate exon's start, and every exon end except the highest-coordinate exon's end, must be unchanged. This preserves whole internal exons and all splice-facing boundaries on both strands. Single-exon transcripts have no protected internal boundaries. Use the transcript boundary datacheck above as well to detect shortening at the outer ends; CDS and UTR features are outside this check.

The CSV contains one row per transcript ID, both chromosome/strand values, exon counts and semicolon-separated `start-end` exon lists. Status is `unchanged_internal_structure`, `changed_boundaries`, `exon_count_mismatch`, `location_mismatch`, `missing_from_extended` or `missing_from_raw`. Status precedence is location mismatch, exon-count mismatch, then changed boundaries. Independent mismatch flags retain simultaneous failures. Boundary details use one-based genomic exon positions, for example `exon_2_start:200->201`; these positions run from low to high coordinates on either strand. Boundary flags and details are blank when exon counts differ, because pairing would be ambiguous. Missing annotation values and unmatched-transcript flags are blank.

The console reports status counts, independent failure-flag counts and the number of failing transcripts (counted once each). Missing transcripts are reported without failing, allowing chromosome-wide raw files to be compared with strand-specific outputs. Exit codes are `0` for pass, `1` for structural or location failures, and `2` for input/output errors. The check rejects malformed transcript/exon records, missing IDs, duplicate transcripts, invalid or overlapping exon intervals, exons outside their transcript interval, inconsistent exon chromosome/strand, orphan exons and transcripts without exons. Exons may precede transcript records. The report cannot overwrite either input. No Nextflow integration is included.

Run both datachecks' automated tests with `python3 -m unittest discover -s tests`.

## Added UTR overlap datacheck

```bash
python3 assets/check_cds_overlaps.py current_annotation_1.gff3 1_positive_extended.gff utr_overlap_comparison.csv \
  --leap-output-dir /path/to/original_leap_run \
  --extension-stats /path/to/leap_run/extended_transcripts_stats
```

This standalone standard-library check accepts plain or gzipped pipeline GFF3 and exact versioned transcript IDs. The third positional argument is a **CSV filename**. The required `--leap-output-dir` supplies the original run's `natural_overlaps/` and `fixed_extensions/` files. Missing files, multiple candidate files or conflicting candidate metadata are input errors. Raw `transcript_type` and UTR features reconstruct LEAP's full-length protein-coding eligibility (both 5′ and 3′ UTRs required). Other transcript types are non-protein-coding for these rules, even when their gene is protein-coding.

The overlap rules come from LEAP's clipping and terminal-safety stages, applied to the checked annotation with the original run context. Clipping protects eligible CDS and non-protein-coding exons from **all genes, including the source gene**. A full-length same-gene CDS contact is therefore assessed using LEAP's grouped annotation boundary and candidate support, rather than an extra per-isoform prohibition.

There is one deliberate stricter check: **within-gene CDS from partial raw protein-coding transcripts remain protected**, even when LEAP filtered those transcripts out. This extra check uses added UTR bases relative to the actual raw isoform. Its existing shared-internal-splice-site and purely `current_annotation` support allowances are retained. Both shared boundaries must be real splice sites; single-exon targets cannot qualify. Leftward extension uses the exon end, rightward extension the exon start, on either strand. Mixed annotation/read support does not grant this extra annotation allowance.

LEAP's two overlap stages are:

- **Clipping:** protect full-length coding models' CDS, single-exon CDS, and non-protein-coding exons. For single-CDS genes in different-gene natural-overlap pairs, apply LEAP's additional flank restrictions when the candidate's complete gene-ID field is that singleton gene. The flank check uses all such special genes' exons on the appropriate side of their CDS, falling back to the CDS-containing exon when that side has no exon.
- **Terminal safety:** test FE–LE and terminal–single-exon contacts. A contact is safe when the ordered candidate gene-ID groups share an ID, or any **zipped ordered** FE/LE gene pair occurs in the recorded natural-overlap file. Overlapping gene spans alone grant no allowance. FE–FE, LE–LE and ordinary internal-exon contacts are outside this stage.

Annotation-supported cross-gene extensions bypass clipping but still fail unsafe terminal contacts. The original fixed BED candidate's support field is authoritative for LEAP clipping. Per-transcript stats cannot override a read-derived candidate to bypass clipping; stats enable the separate partial-CDS allowance. Other raw exon contacts are retained as **non-failing audit findings**, including unprotected protein-coding UTRs and ineligible partial models. The partial-CDS check is the only intentional additional overlap restriction.

The clipping guard includes the grouped original boundary recorded in the fixed-extension BED, which can differ from a shorter isoform's raw boundary. It uses inclusive comparisons, as `stop_cds_overlap.py` does. Terminal safety checks whole candidate exons, including their already annotated portions. Contacts use LEAP's positive-length BED overlap predicate: touching endpoints alone do not qualify. Reports use inclusive genomic coordinates; clipping findings can therefore include the grouped original boundary. The terminal tests use the annotation being checked, with original run metadata, so internally processed files can be evaluated without restoring reverted candidates. Use **processed stats** with processed annotations.

The CSV has one row per source/target, target-feature type and rule. `check_rule` explains each result; `intentional_partial_cds_check=true` identifies the extra partial-CDS rule; `gene_relationship`, `target_is_partial` and `target_missing_utr` identify the relationship and partial-CDS evidence. Separate `annotation_allowed_intervals`, `splice_site_allowed_intervals`, `leap_allowed_intervals` and `failing_intervals` retain the evidence and mixed outcomes. Counts are merged bases within each row; different rules can overlap, so row/base totals are not unique transcript counts. Failures have status `new_cds_overlap`, `new_cross_gene_cds_overlap`, `new_cross_gene_exon_overlap`, `new_within_gene_exon_overlap` or `unsafe_terminal_overlap`. Allowances use `allowed_current_annotation`, `allowed_shared_splice_site`, `allowed_multiple_reasons` or `allowed_leap_overlap`. Unmatched extended UTR transcripts are `missing_from_raw`, counted separately without failing.

Stats accept a TSV file/directory, plain or gzipped. Gene, transcript location and extension lengths must match the checked annotation; stale or conflicting stats are rejected before the report is opened. Missing entries grant no within-gene annotation allowance. Only same-chromosome, same-strand contacts count. Self-comparisons are excluded. The parser checks IDs, locations, containing exons and intervals; duplicate intervals are deduplicated and overlapping distinct exons rejected. Unrelated transcript-only raw records are accepted. Reports cannot overwrite annotations, stats or required run-context files.

Exit codes are `0` for pass, `1` for biological failures and `2` for input/output errors. The console reports `leap_rule_failures` and `partial_cds_failures` separately, as well as evaluated UTR transcripts and report rows by outcome. Production LEAP, internal-extension processing and Nextflow integration are unchanged. Inspectable examples, including all seven retained within-gene partial-CDS failures, are in [tests/data/leap_cross_gene](tests/data/leap_cross_gene/README.md). Earlier within-gene examples remain in [tests/data/cds_overlaps](tests/data/cds_overlaps/README.md). See [the LEAP overlap-rule contract](docs/leap_overlap_rules.md) for the source mapping and intentional exception. Run `python3 -m unittest discover -s tests` for regression and production-predicate parity checks.

## Internal-extension processing asset

```bash
python3 assets/process_internal_extensions.py /path/to/leap_run /path/to/new_processed_run
```

This standalone asset brings LEAP2GO's existing internal/external classification into LEAP without changing the Nextflow workflow. It restores internal terminal exons from raw annotation, saves original candidates separately, and writes strand-correct extension stats, including zero-extension rows for fully restored transcripts. Transcripts without CDS are excluded as before. It requires only Python's standard library and local `sort`, and leaves input files untouched.

See [the processing guide](docs/internal_extension_processing.md) for classification rules, output layout and the subsequent datacheck/JBrowse workflow. Small positive/negative-strand examples are in `tests/data/internal_extensions/`.

## JBrowse review bundles

Rebuild the visual review with `assets/build_jbrowse_review.py`:

```bash
python3 assets/build_jbrowse_review.py \
  --raw-annotation /path/to/leap_run/raw_annotation \
  --extended-annotation /path/to/leap_run/extended_transcripts \
  --reports /path/to/cds_overlap_reports \
  --output /path/to/new_jbrowse_review
```

Open the generated `LEAP-overlap-review.jbrowse` in JBrowse Desktop. The bundle includes indexed raw/extended annotations, CDS/UTR context, red/orange failure markers, green shared-site, teal annotation-supported and blue LEAP audit markers, example sessions and a table of review locations. It uses public GRCh38 and requires `sort`, `bgzip` and `tabix`. Inputs may be files or directories; the output directory must be new.

For future datachecks, omit `--reports` and add marker GFF3 files with repeatable `--extra-track 'NAME=PATH'` arguments. The full workflow, marker format and reuse notes are in [the JBrowse review guide](docs/jbrowse_reviews.md).

## Author

Ryan Merritt

## License

See [LICENSE](/Users/rmerritt/my_app/LEAP/LICENSE).
