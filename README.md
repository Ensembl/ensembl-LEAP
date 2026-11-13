# LEAP

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A523.04.0-23aa8c.svg)](https://www.nextflow.io/)
[![Python](https://img.shields.io/badge/python-%E2%89%A53.10-3776AB?logo=python&logoColor=white)](https://www.python.org/)
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
- [Automatic datachecks](#automatic-datachecks)
- [Manual reviewer preparation](#manual-reviewer-preparation)
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
- Python 3.10+ with SQLAlchemy 2 and PyMySQL 1 for annotation database access
- Access to an Ensembl-style annotation database compatible with the bundled Python exporter

Several modules declare Biocontainers images for tools such as `samtools`, `bedtools`, `bedToGenePred`, and `genePredToGtf`. The default config enables both Singularity and Conda:

```groovy
singularity.enabled = true
conda.enabled = true
```

See [the exporter migration guide](docs/db2gxf_upgrade.md) for JSON setup, bundled exporter provenance, runtime requirements and chromosome 21 validation. The old `--host`, `--port`, `--user`, `--pass` and `--dbname` parameters have been replaced by `--db_secrets`.

## Inputs

LEAP expects:

- A coordinate-sorted long-read RNA-seq BAM file.
- A BAM index file (`.bai` or `.csi`).
- CAGE peaks in BED-like format for 5' validation.
- polyA sites in BED-like format for 3' validation.
- A database connection JSON file for the current annotation source (`--db_secrets`).
- An output directory for published results.
- A species-matched FASTA `.fai` or chromosome-length file for automatic datachecks, with exact annotation sequence names.

The CAGE and polyA files are split by chromosome and strand using the chromosomes found in the read BAM.

## Quick Start

Edit [params.config](/Users/rmerritt/my_app/LEAP/params.config) or provide parameters on the command line.

```groovy
params {
    output_directory = "/path/to/leap_results"
    datacheck_reference_lengths = "/path/to/matching_reference.fa.fai"

    reads = "/path/to/long_reads.bam"
    reads_index = "/path/to/long_reads.bam.bai"
    cage = "/path/to/cage.bed"
    polya = "/path/to/polya.bed"

    minimum_read_support = 2
    overlap_penalty = 10
    extend_on_current_annotation = true

    db_secrets = "/private/path/db_secrets.json"
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
  --datacheck_reference_lengths /path/to/matching_reference.fa.fai \
  --minimum_read_support 2 \
  --overlap_penalty 10 \
  --extend_on_current_annotation true \
  --db_secrets /private/path/db_secrets.json \
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
| `run_datachecks` | No | Run all ten validation modules automatically; default `true`. |
| `datacheck_reference_lengths` | When checks enabled | Matching FASTA `.fai` or two-column chromosome-length file; exact sequence names required. |
| `minimum_read_support` | Yes | Minimum read support passed to read unique-splice-site extraction. |
| `overlap_penalty` | Yes | Alpha/penalty value passed to unsafe-overlap fixing. |
| `extend_on_current_annotation` | Yes | Controls whether candidate read extensions may be matched against current annotation. |
| `db_secrets` | Yes | Connection JSON: host, port, user, database and password; see the example and migration guide. |
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
4. Classify full-length and partial coding transcripts, keeping partials separately for protection.
5. Split the original full-length annotation into multi-exon and single-exon transcripts.
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
9. Restrict already selected read extensions against partial reference models.
10. Build final extended transcript GFF files and per-chromosome statistics, excluding partial transcript models.

Partial models are protection context only. Endpoint generation, shared-splice
annotation boundaries, gene groups, natural-overlap pairs and the original
single-CDS-gene eligibility use the original full-length reference. After normal
overlap fixing, partial CDS and applicable exon protections can shorten a
read-supported endpoint or restore its original annotation fallback. They cannot
create endpoints, lengthen a selected extension or free space for another end
to grow. Existing annotation-only extensions are preserved. Only full-length
multi-exon models are emitted, matching the original output eligibility.

`filter_out_partials.py` retains its `_FT.gff`, `_PT.gff` and `_NPC.gff` outputs
and additionally writes `_PC.gff`, containing both full-length and partial
protein-coding models, for reference/compatibility only. Nextflow uses `_FT.gff`
for normal processing and passes `_PT.gff` separately to
`protect_partial_extensions.py` and the final writer. Direct calls to
`build_chromosome_gff_and_stats.py` now require `--partial-annotation <PT.gff>`;
an empty file is valid when no partials exist. Partial IDs are matched exactly,
including version, chromosome and strand. All records belonging to those models
are omitted before applying extensions, so they cannot contribute stats or final
gene boundaries. See [partial-transcript processing](docs/partial_transcript_processing.md)
for the protection contract, command, tests and chromosome 12 validation.

The final writer takes each generated gene row's GFF `ID` from its transcripts' gene `Parent`, preserving identifiers such as `_PAR_Y` copies. The stable `gene_id` attribute and all transcript/child IDs and Parent links are unchanged. Gene bounds are grouped by this GFF parent identity. Legacy direct inputs lacking `Parent` retain the previous `gene_id` fallback; empty or multiple gene Parents are rejected.

A more detailed workflow diagram is available in [BIOINFORMATICS_FLOWCHART.md](/Users/rmerritt/my_app/LEAP/BIOINFORMATICS_FLOWCHART.md).

## Outputs

Selected outputs are published under `params.output_directory`.

| Directory | Contents |
| --- | --- |
| `raw_annotation/` | Per-chromosome current annotation GFF3 files from the configured database. |
| `natural_overlaps/` | Existing annotation overlaps that are treated as natural/known overlaps. |
| `extension_safety/` | Safety-check outputs, including safe and unsafe extension BED files. |
| `unsafe_overlaps/` | Rejected unsafe overlaps, named like `<chromosome>_<strand>_rejected.bed`. |
| `datachecks/` | Per-check CSVs, logs and statuses, plus combined HTML/CSV/JSON summaries and final verdict. |
| `fixed_extensions_before_partial_protection/` | Original selected first/last-exon candidates before the added partial checks. |
| `fixed_extensions/` | Final first/last-exon candidates after protection; these correspond to the final GFF and stats. |
| `partial_protection/` | TSV audit of shortened candidates, original endpoints, annotation fallback and blocking partial features. |
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
│   ├── gxf_dump.py
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

## Automatic datachecks

All ten datachecks now run automatically on standard LEAP outputs. Supply the matching reference lengths when launching LEAP:

```bash
nextflow run main.nf -profile hpc_slurm \
  --db_secrets /private/path/db_secrets.json \
  --datacheck_reference_lengths /path/to/matching_reference.fa.fai
```

A two-column chromosome-length file is also accepted, including gzip. Names must match exactly: `1` and `chr1` are distinct. The reference parameter is required while checks are enabled and is validated before production tasks are submitted. Set `--run_datachecks false` only to explicitly disable automatic validation; the default is `true`.

Each chromosome/strand gets separate boundary, internal-exon, complete-preservation, stats/GFF, structural-integrity, completeness, partial-protection, reference-coordinate, CDS/UTR-overlap and duplicate-recommendation tasks. Inputs come directly from declared workflow channels, independently of asynchronous published copies. Expected partitions come from the chromosomes selected from the read BAM, with both strands required. Missing output partitions therefore reach completeness rather than disappearing through joins.

Results are saved under `datachecks/<chromosome>_<strand>_<check>/`: `report.csv`, `stdout.txt`, `stderr.txt` and a uniquely named `.status.json`. The combined `datachecks/summary/summary.html`, `summary.csv` and `summary.json` distinguish `PASS`, `FAIL`, `INPUT_ERROR`, `NOT_APPLICABLE` and non-failing `REVIEW`. All checks finish and their reports and summary are published before the final gate marks LEAP as failed. Input errors remain distinct from biological failures. Missing/duplicate/unexpected result records and unexpected workflow scopes also fail validation. Partial-CDS contacts within verified original grouped annotation boundaries are reported as allowed; new read-driven contacts beyond those boundaries remain failures unless an existing explicit allowance applies.

For valid empty partitions, completeness, partial protection, reference coordinates and duplicate recommendations still run. Comparison/structural checks are marked `NOT_APPLICABLE` only after a scoped completeness preflight confirms zero eligible transcripts, empty final GFF and no stats models. Missing eligible output or malformed inputs are never hidden as an empty-partition pass. The reference check covers both the published raw annotation and final GFF. The stats/GFF check always receives selected-candidate provenance. Internal-extension processing remains standalone; the automatic suite checks the standard LEAP output stage.

The new modules require only Python 3 on the execution environment. Annotation, evidence and extension-selection logic are unchanged. Existing standalone commands below remain available. See [the integration guide](docs/datachecks.md) for workflow wiring, output meanings and real Nextflow smoke tests.

## Transcript boundary datacheck

Run the standalone check on a raw annotation and its corresponding extended annotation:

```bash
python3 assets/check_for_shorter_transcripts.py current_annotation_1.gff3 1_positive_extended.gff comparison.csv
```

Both inputs accept plain GFF3 or `.gz` files. The check uses `transcript` features and exact `transcript_id` values, including version suffixes. Each matched transcript must retain its chromosome and strand and satisfy `extended_start <= raw_start` and `extended_end >= raw_end`, on either strand. Growth at one end does not compensate for shortening at the other. Exon structure is not checked.

The CSV includes all IDs from both inputs, both sets of chromosome, strand, start, end and inclusive genomic length (`end - start + 1`), and boundary-shortening and location-mismatch flags. Status is `unchanged`, `extended`, `shortened`, `location_mismatch`, `missing_from_extended` or `missing_from_raw`. Location mismatches take precedence over boundary comparisons; shortening flags are blank for these rows and for unmatched IDs. Missing annotation values are blank.

The console reports a pass/fail result and counts for every status. Missing transcripts are reported separately and do not fail the check, since raw annotations may contain transcripts filtered before extension. A chromosome-wide raw file can therefore be compared with a strand-specific extended file.

Exit codes are `0` for pass, `1` for shortened transcripts or chromosome/strand mismatches, and `2` for input or execution errors. Invalid transcript records, missing or duplicate IDs, invalid coordinates and files without transcript records are rejected. This script also runs through the automatic datacheck module above.

Run the focused tests with `python3 -m unittest discover -s tests`.

## Internal exon structure datacheck

```bash
python3 assets/check_internal_exon_positions.py current_annotation_1.gff3 1_positive_extended.gff exon_comparison.csv
```

This standalone, standard-library check accepts plain or gzipped GFF3 and matches exact `transcript_id` values. It sorts exons by genomic coordinates, ignoring file order, exon IDs and exon numbers. Matched transcripts must retain chromosome, strand and exon count. Every exon start except the lowest-coordinate exon's start, and every exon end except the highest-coordinate exon's end, must be unchanged. This preserves whole internal exons and all splice-facing boundaries on both strands. Single-exon transcripts have no protected internal boundaries. Use the transcript boundary datacheck above as well to detect shortening at the outer ends; CDS and UTR features are outside this check.

The CSV contains one row per transcript ID, both chromosome/strand values, exon counts and semicolon-separated `start-end` exon lists. Status is `unchanged_internal_structure`, `changed_boundaries`, `exon_count_mismatch`, `location_mismatch`, `missing_from_extended` or `missing_from_raw`. Status precedence is location mismatch, exon-count mismatch, then changed boundaries. Independent mismatch flags retain simultaneous failures. Boundary details use one-based genomic exon positions, for example `exon_2_start:200->201`; these positions run from low to high coordinates on either strand. Boundary flags and details are blank when exon counts differ, because pairing would be ambiguous. Missing annotation values and unmatched-transcript flags are blank.

The console reports status counts, independent failure-flag counts and the number of failing transcripts (counted once each). Missing transcripts are reported without failing, allowing chromosome-wide raw files to be compared with strand-specific outputs. Exit codes are `0` for pass, `1` for structural or location failures, and `2` for input/output errors. The check rejects malformed transcript/exon records, missing IDs, duplicate transcripts, invalid or overlapping exon intervals, exons outside their transcript interval, inconsistent exon chromosome/strand, orphan exons and transcripts without exons. Exons may precede transcript records. The report cannot overwrite either input. It also runs through the automatic datacheck module above.

Run both datachecks' automated tests with `python3 -m unittest discover -s tests`.

## Manual reviewer preparation

LEAP now prepares `manual_review/` automatically, with raw, original LEAP, processed external and original internal-candidate annotations. Each group has positive, negative and unsplit GFF3 views across the selected chromosomes. Internal tracks show what the models would look like if their internal extensions were retained; external tracks restore those ends using the existing processor. The bundle includes HAVANA JBrowser samplesheets, metadata, corrected external stats/classifications and editable `run_jbrowser.sbatch` (public FTP) and `run_jbrowser_offline.sbatch` launchers. The offline launcher refreshes a separate samplesheet with every row marked offline, using the shared reference/path settings, and needs no FTP settings. The final datacheck gate waits for this bundle.

Supply `--review_reference /server/reference/genome.fa` (or a gzipped FASTA/2bit), alongside the matching `--datacheck_reference_lengths`. Set `--review_jbrowser_pipeline`, `--review_public_base_url`, `--review_ftp_target_dir` and `--review_hub_email` to prefill server settings; unknown settings remain editable. Use `--review_genome_id` and `--review_species` for explicit assembly/species names. Disable preparation with `--prepare_manual_review false`. LEAP generates the launch files; submission and FTP publishing remain manual.

See [manual review preparation](docs/manual_review.md) for the 12 views, launch commands and preparing an existing run without rerunning LEAP.

## Duplicate transcript recommendations

```bash
python3 assets/check_duplicate_transcripts.py /path/to/run/extended_transcripts duplicate_transcripts.csv --recommend-removals
# A single plain or gzipped GFF/GFF3 file is also accepted.
```

Different exact `transcript_id` values are grouped when their chromosome, strand, complete exon coordinates and complete CDS coordinates **and phases** are identical. Terminal exon boundaries must match too; sharing only splice junctions or transcript bounds is insufficient. Exons and CDS segments are sorted by genomic position, so record order, exon IDs, exon numbers and other metadata do not affect matching. Gene membership does not restrict grouping. Noncoding models can match other models with the same exon chain and no CDS; coding and noncoding models cannot match each other. X/Y copies remain separate by chromosome. This tests genomic structure, not reference sequence identity, and does not recommend which transcript to retain.

The CSV contains only duplicate groups, with one row per member. Filter `group_id` to find all members of a group. Each row includes group size, exact versioned transcript ID, gene identifiers/names, chromosome/strand, transcript bounds, ordered exon and CDS lists, CDS phases, inclusive spliced/CDS lengths, source file and transcript line, and a `jbrowse_label` in `chromosome:start-end` format. Group labels depend on structure and stay stable across input ordering and subset reports. The console counts groups, transcripts belonging to groups, extra members beyond one per group (`redundant_transcripts`) and transcripts outside groups (`unique_transcripts`).

Directory mode checks all immediate `.gff`, `.gff3` and gzip equivalents together. Ambiguous plain/gzip copies, repeated transcript records (including across files on the same chromosome), malformed records, missing IDs, orphan exon/CDS records, exonless transcripts, overlapping/duplicate exon or CDS intervals, inconsistent locations/gene IDs, transcript bounds inconsistent with exon spans and CDS outside exons are input errors. Children may precede transcripts. Empty/header-only files are accepted. Inputs are validated before opening the report, and overwrite protection covers symbolic and hard links. Exit codes are `0` for no duplicate groups, `1` for duplicate groups found, and `2` for input/execution errors. The tenth automatic module enables recommendation mode and records duplicate findings as non-failing `REVIEW`; malformed inputs still fail validation. Without `--recommend-removals`, the standalone command retains detection-only behaviour.

Recommendation mode keeps every MANE Select/Plus Clinical member. Otherwise it prefers non-canonical models when available, then models carrying `LEAP_blacklisted` or the original upstream-ATG/CAGE-supported-TSS/RP-supported-TIS tags, then the lowest numeric stable ID. Canonical-only groups retain one representative. The same priorities apply across genes; chromosome and strand still separate groups. Dates and transcript versions never break ties. Unreadable or tied numeric IDs produce manual-review candidates. CSV fields include observed tags, classification, `keep`/`remove_from_gff`/`review`, reasons, retained member IDs and a cross-gene flag. Recommendations apply to LEAP GFF entries only and do not edit GFFs or databases. Historical files lacking MANE metadata cannot support reliable MANE-priority decisions.

Run focused tests with `python3 -m unittest discover -s tests -p 'test_*duplicate*.py'`.

The standalone [annotation uploader](assets/update_LEAPed_annotation.pl) now reuses matching database models tagged `Ensembl_canonical_extended` when processing MANE extensions. An existing model already covering the requested ends is left intact; a shorter matching model is extended in place, keeping the farther endpoint at each end. The uploader creates a copy only when no matching canonical extended model exists. See [uploading LEAP annotations](docs/uploading_leap_annotations.md) for matching rules, dry-run usage and test limitations. Uploading remains separate from the pipeline and datachecks.

## Complete transcript preservation datacheck

```bash
python3 assets/check_transcript_preservation.py current_annotation_1.gff3 1_positive_extended.gff preservation.csv
```

This standalone, standard-library check accepts plain or gzipped pipeline GFF3 and matches exact versioned `transcript_id` values. It checks every transcript-associated feature, including internal and terminal exons, CDS, start/stop codons, UTRs and unfamiliar feature types. All original feature types, counts, coordinates, chromosome, strand, source, score, phase and attributes must be preserved. This includes feature IDs, `Parent`, gene membership, exon IDs and exon numbers. Only the LEAP-generated attributes `leap_tag` and `transcript_class` are excluded from comparison. Record order, attribute order, headers and gene rows are not compared; ordinary attribute values are compared exactly.

Only outward changes to the transcript span, the lowest-coordinate exon's start and the highest-coordinate exon's end are permitted. This also applies to single-exon transcripts. Corresponding terminal UTR features must follow those outer boundaries, using the correct 5′/3′ interpretation on each strand. Internal UTR segments and UTR boundaries facing CDS cannot move. For noncoding models without UTR features, exon/transcript extension alone is valid. Conceptually removing only these permitted terminal additions must recover the complete original feature collection. New features, including features wholly within added sequence, are failures; no records are silently deduplicated. Transcript bounds must agree with the outer exon bounds.

The CSV contains one row per ID with status `unchanged`, `extended`, `changed_features`, `missing_from_extended` or `missing_from_raw`, a `failed` flag, both annotations' locations and JSON feature-type counts. Independent flags identify location mismatches, shortening, feature-count changes, inconsistent extensions and differences in complete feature records. `extension_details` is a JSON list explaining transcript/exon/UTR inconsistencies. `feature_differences` is a JSON object containing original-only records, extended-only records after removing permitted additions, and actual extended records for affected feature types; records include coordinates, attributes, source, score and phase. These details avoid ambiguous pairing when IDs or feature counts change.

Raw transcripts omitted from output are reported without failing because LEAP filters annotations. An output transcript absent from raw **does fail**, including a changed version suffix. A chromosome-wide raw file can be compared with a strand-specific output; raw transcripts on the other strand are non-failing unmatched records. Summary counts and independent flag counts are printed. Exit codes are `0` for pass, `1` for preservation failures or unexpected output IDs, and `2` for input/output errors.

Input validation requires `transcript_id` and `gene_id` on transcript-associated records and rejects malformed records, duplicate transcript/feature records, invalid coordinates, scores or phases, inconsistent locations/gene IDs, overlapping exons, children without transcript records, unresolved child `Parent` references, features outside their transcript, CDS/codon/UTR features without a containing exon, and transcripts without exons. Children may precede transcripts. Both inputs are validated before opening the report, and reports cannot overwrite either input, including through symbolic or hard links. Invalid unmatched raw models remain input errors rather than being silently filtered out.

This check includes the preservation rules covered by the boundary and internal-exon checks and adds complete feature/metadata comparison; those existing scripts remain available unchanged. It does not replace the overlap datacheck, validate extension evidence, or compare gene-row bounds. It also runs through the automatic datacheck module above. Run its focused tests with `python3 -m unittest discover -s tests -p test_check_transcript_preservation.py`.

## Stats and GFF consistency datacheck

```bash
python3 assets/check_stats_gff_consistency.py current_annotation_1.gff3 1_positive_extended.gff stats_comparison.csv \
  --extension-stats 1_positive_stats.tsv \
  --leap-output-dir /path/to/leap_run
```

This standalone, standard-library check compares final transcript/exon coordinates, the corresponding stats and transcript `leap_tag` values. It accepts plain/gzipped GFF3 and TSV, plus a stats directory containing only the files corresponding to the checked final annotation. Exact versioned IDs are keyed by chromosome, so X/Y PAR copies remain distinct. Every final transcript must have a raw counterpart and exactly one stats row. Extra stats rows and unexpected output IDs fail; raw models intentionally absent from output are non-failing. The final transcript span must equal its outer exon boundaries; stats locations and gene IDs must agree with the GFF, and chromosome/strand/gene changes or inward boundaries fail. Raw coordinates come from transcript records; validating unrelated raw child features is left to the preservation check.

The expected 5′ extension is `raw_start - final_start` on `+` and `final_end - raw_end` on `-`; 3′ uses the opposite end. These are added-base counts, without an extra `+1`. A nonzero extension requires `reads` and/or `current_annotation` support. Read-supported entries require a positive count. Annotation-supported entries can have **zero or nonzero** counts: LEAP reports reads supporting the splice site, not necessarily reads supplying that endpoint. Aligned comma-separated support/read-count lists are supported and validated in order. A zero-extension end requires support `none` and count `0`, except blacklisted 5′ ends, which require zero extension, support `blacklisted` and literal read-count value `blacklisted`. Blacklisting is invalid at 3′. `LEAP_5` and `LEAP_3` must reflect actual extensions; `LEAP_blacklisted` must match the stats blacklist declaration; `LEAP_none` is used only when no other tag applies. Tag order is ignored; missing, unknown, duplicate and contradictory tags fail.

`--extension-stats` is required. Optional `--leap-output-dir` verifies nonzero ends against the original run's `fixed_extensions/`: endpoint, gene group, support category and splice-site read count must match the selected FE/LE record. This replays the final writer's ordered **last-record-per-splice-site** lookup, including coincident annotation/read candidates; the overlap check now uses the same ordered selection. Candidate fields are used as written by LEAP without a BED/GFF offset conversion. Missing/ambiguous candidate files or malformed records are input errors. With processed GFFs, supply **processed stats** and the original run directory: zero/restored ends are checked for cleared support/counts/tags, while surviving extensions retain their original candidate provenance. Without run context, syntax and internal consistency are checked, but support origins and exact numeric read counts cannot be established from a GFF alone. This check does not recount BAM reads or independently establish blacklist membership or annotation-donor eligibility.

The CSV includes all raw/final/stats keys, status (`consistent`, `inconsistent`, `missing_from_extended`), both GFF locations, reported and expected lengths, support/counts, observed/expected tags, independent mismatch flags and JSON `details` explaining each discrepancy. `candidate_provenance_checked` identifies rows where a nonzero end was checked against candidates. Summary counts and flag counts are printed; exit codes are `0` for pass, `1` for inconsistencies and `2` for malformed inputs or execution errors. Duplicated stats keys/columns, missing required fields, malformed integers/support lists, malformed GFF records, duplicate transcript records, orphan/missing/overlapping final exons and missing context files are rejected. Inputs are loaded before opening the report, which cannot overwrite annotations, stats or candidates, including through links.

The standalone rules and production LEAP logic are unchanged; automatic execution is provided by the datacheck subworkflow above. Run focused tests with `python3 -m unittest discover -s tests -p test_check_stats_gff_consistency.py`.

## GFF structural integrity datacheck

```bash
python3 assets/check_gff_structural_integrity.py 1_positive_extended.gff structural_issues.csv
```

This standalone standard-library check inspects one plain/gzipped GFF3 without requiring raw annotation or stats. It resolves the actual `ID` / `Parent` graph, independently of file order, and checks the gene → transcript/mRNA → child hierarchy. Gene and transcript IDs must identify single records. Leaf IDs are optional; compatible multi-part features such as CDS, UTRs and split codons may share an ID, including varying CDS phases and exon annotations. Exact duplicate records, conflicting ID reuse, duplicate Parent references, missing/unresolved/ambiguous parents and parent cycles fail. CDS segments are not rejected solely for overlap: programmed frameshifts can use overlapping CDS records. GFF IDs have file-wide scope; merged annotations must disambiguate colliding IDs, including PAR copies. Ordinary `exon_id` attributes may be shared and are not treated as global GFF IDs. Percent-encoded IDs are decoded when resolving links, with Parent lists split before decoding.

Each transcript/mRNA must have one gene Parent and at least one exon. Exons must have transcript/mRNA parents and must not overlap within a transcript. CDS, codon and UTR segments must have a transcript ancestor and lie within an exon of **each** associated transcript. Shared exons and other multiple-parent features are supported; every declared Parent must resolve. Nested unfamiliar features are checked through their ancestry. `gene_id` / `transcript_id` attributes, when supplied, must agree with the resolved gene/transcript metadata; these attributes are not substituted for missing GFF IDs or Parent links. Genes without transcript descendants also fail in this LEAP annotation profile.

Every child must have the same chromosome and strand as its parent and fit inside the parent's interval(s). The check also verifies containment in all ancestor transcript and gene bounds, so a misplaced descendant is reported even through an intermediate feature. For discontinuous parents, a child cannot span a gap between parent segments. Parent bounds may include margins: this check requires containment, not exact equality to the child union. Containment, root genes and single-gene transcript parents are requirements of this LEAP annotation profile, rather than rules for every possible GFF3 ontology feature. The stats/GFF consistency check separately enforces LEAP's exact transcript/exon spans. Numeric coordinates, scores, strand and phase syntax are validated; CDS records without numeric phase are reported as structural failures. Headers and embedded FASTA are not validated, and this is not a full Sequence Ontology, reference-length or translation validator. ID/Parent handling follows the [Sequence Ontology GFF3 specification](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md).

Internally processed outputs deliberately omit gene rows. Check those with the explicit option:

```bash
python3 assets/check_gff_structural_integrity.py processed/1_positive_extended.gff structural_issues.csv --allow-missing-genes
```

This option allows unresolved single gene Parents on transcript/mRNA records **only when the entire file has no gene rows**. All transcript/child links, exon structure and transcript containment remain required. It does not excuse a partially missing set of genes. Conflicting `gene_id` declarations for a shared omitted gene Parent still fail. The console reports `gene_boundaries_checked: False` and the number of omitted gene-parent checks; gene containment cannot be verified without gene records. No gene boundaries are inferred or repaired.

The CSV contains one row per structural issue: issue code, source line, feature type/ID, transcript/gene attributes, chromosome, coordinates, strand, relevant parent ID and JSON diagnostic details. A passing report contains only its header. The console reports parsed record/model counts, issue counts and gene-check coverage. Structural violations return `1`; a valid annotation returns `0`; malformed/unreadable/empty input or execution errors return `2`. Syntax errors are diagnosed before opening the report. Reports cannot overwrite the input, including through symbolic/hard links. Structural problems such as an exonless transcript are reported while the remaining readable records are still evaluated. Files without transcript/mRNA records are input errors.

Existing checks and production code remain unchanged; automatic execution is provided by the datacheck subworkflow above. Run focused tests with `python3 -m unittest discover -s tests -p test_check_gff_structural_integrity.py`.

## Output completeness datacheck

Check a complete run with an explicit chromosome selection:

```bash
python3 assets/check_output_completeness.py run/raw_annotation run/extended_transcripts completeness.csv \
  --extension-stats run/extended_transcripts_stats --chromosomes 1 2 3 X Y
```

Directory comparisons default to both strands. Use `--strand positive` or `negative` to select one strand. Single-file comparisons require exactly one chromosome and an explicit strand:

```bash
python3 assets/check_output_completeness.py current_annotation_Y.gff3 Y_positive_extended.gff Y_completeness.csv \
  --extension-stats Y_positive_stats.tsv --chromosomes Y --strand positive
```

This standalone standard-library check derives eligibility independently from the supplied published raw annotation, after upstream OTT/readthrough filtering. Every `protein_coding` transcript with both `five_prime_UTR` and `three_prime_UTR` features and more than one exon must occur in **both** final GFF and stats. Noncoding, partial and single-exon models are expected exclusions, with explicit reasons. Blacklisting, lack of reads, rejected candidates and zero extension do not excuse missing eligible models. Internal classification alone never excludes a model. Use `--stage processed` for internally processed output: this additionally excludes otherwise eligible models without CDS, matching the processor. The default `--stage leap` does not require CDS.

Identity is the exact versioned `transcript_id`, chromosome and strand. GFF primary IDs, including `_PAR_Y`, are retained as evidence and are not substituted for stable transcript IDs. Raw records outside the selected scope are ignored; final GFF/stats records outside the selected scope or placed in the wrong chromosome/strand file fail. Stats have no strand column, so their file slot supplies the strand. This checks downstream completeness against current eligibility rules, not database-export completeness or historical-release inclusion policy.

Use three corresponding files or three directories. Inputs accept plain/gzipped files. Directory slots require exactly one `current_annotation_<chromosome>.gff3` (also `.gff`), `<chromosome>_<positive|negative>_extended.gff` (also `.gff3`) and `<chromosome>_<positive|negative>_stats.tsv`, with optional `.gz` suffixes. Missing selected files fail even for zero-eligible partitions; multiple files matching a slot are input errors. Existing output directories may be empty. Empty final GFFs and header-only stats pass only when the scope has no eligible models. Empty/unclassifiable raw annotation is an input error; a valid chromosome raw file containing only opposite-strand models can establish a zero-eligible selected strand.

The CSV combines one row per chromosome/strand/transcript identity with file findings. It records eligibility, exclusion reasons, independent GFF/stats presence, raw feature counts, GFF primary IDs, paths and explanations. Statuses include `retained`, `expected_exclusion`, `missing_expected`, `unexpected_retained`, `outside_selected_scope`, `raw_unavailable`, missing-file statuses and misplaced-file findings. Missing raw files leave membership explicitly unverifiable. Console counts include expected retention, actual membership, statuses and exclusion reasons; exclusion reasons can overlap. File findings and transcript failures are counted separately.

Exit codes are `0` for agreement, `1` for missing/unexpected models or required files, and `2` for malformed, unreadable or ambiguous inputs. Malformed records, missing/conflicting classification metadata, missing IDs, duplicate transcripts/features/stats, orphan children, exonless transcripts, overlapping exons, invalid intervals and misplaced child intervals are input errors. Raw gene rows are outside transcript eligibility. This does not replace the structural, preservation or stats/GFF datachecks: support values and feature preservation are their responsibility. Inputs are validated before the report is opened, and overwrite protection includes symbolic and hard links. Production logic and existing checks remain unchanged; automatic execution is provided by the datacheck subworkflow above.

Run focused tests with `python3 -m unittest discover -s tests -p test_check_output_completeness.py`.

## Partial-protection audit datacheck

```bash
python3 assets/check_partial_protection.py /path/to/leap_run protection_comparison.csv \
  --chromosomes 1 2 X Y --strand both
```

This standalone standard-library audit compares original selected candidates with protected candidates, rather than only inspecting the final GFF. Chromosomes are required; `--strand` is `positive`, `negative` or `both` (default). Each selected chromosome/strand requires its published `raw_annotation/`, original FE/LE BEDs in `fixed_extensions_before_partial_protection/`, final FE/LE BEDs in `fixed_extensions/`, `natural_overlaps/<chromosome>_<strand>_NO.txt` and `partial_protection/<chromosome>_<strand>_partial_protection.tsv`. Existing pipeline filename prefixes are accepted. Each slot must have exactly one file, including plain/gzip alternatives. Missing or ambiguous context is an input error, so runs made before protection-stage context was published cannot be marked as passing.

The check independently classifies published raw models and derives the expected transition without importing or invoking production protection code. Full-length multi-exon protein-coding models determine original single-CDS-gene eligibility; partial protein-coding models supply additional CDS, applicable single-CDS-gene exon flanks, opposite-terminal and single-exon obstacles. Terminal checks preserve LEAP's shared-gene and ordered zipped natural-pair allowances. CDS/flank contacts include the annotation boundary; terminal contacts require positive-length overlap. Exon encounter order preserves the production FE/LE identity, including on the negative strand. This does not introduce blanket protection of every partial exon.

Read candidates may stay unchanged, use an existing shorter safe alternative in its recorded selection order, or return to the original grouped annotation fallback. Stale longer alternatives cannot be revived; no new endpoint may be invented. The fallback remains permitted even if that original annotation boundary contacts a partial feature. Annotation-only and single-exon context rows stay unchanged. Splice sites, ordered gene groups, splice-site read counts and other protected fields must remain identical; only the endpoint and permitted alternative-list/count updates may change. Candidate fallback retains its splice-site count; final GFF/stats support clearing remains the separate stats/GFF check's responsibility.

Candidate order and multiplicity are significant, because the writer uses the last record for a splice site. Legitimate duplicate rows are retained and compared by one-based record position. Unchanged records must also retain their original text, including line endings. Recorded audit rows must agree with independently derived outcomes and blocker details, including multiplicity; audit row order and blocker-list order may differ. Missing, extra or incorrect audit claims fail.

The CSV includes every original candidate, extra output candidates and audit/file findings, with original/expected/actual endpoints and BED fields, annotation fallback, outcome, blocker evidence, source paths and independent flags. Statuses include `unchanged`, `shorter_read_candidate`, `annotation_fallback`, `invalid_transition`, missing/extra candidates, `candidate_order_mismatch`, `missing_audit` and `unexpected_audit`. Console counts distinguish candidate failures, audit findings and file findings; several flags can describe one candidate. Passing requires agreement at both candidate and audit levels.

Exit codes are `0` for pass, `1` for transition/evidence violations and `2` for missing, malformed, unreadable or ambiguous inputs or execution errors. Invalid raw classification/location metadata, duplicate transcripts, orphan children, exonless models, invalid/overlapping exons and malformed BED/audit/natural-pair files are rejected. Raw records outside the selected chromosome/strand scope are ignored. A temporary report is assembled before opening the destination, so input errors preserve any existing report. Overwrite protection includes symbolic and hard links. The standalone rules and production code remain unchanged; automatic execution is provided by the datacheck subworkflow above.

Run focused tests with `python3 -m unittest discover -s tests -p test_check_partial_protection.py`.

## Reference-coordinate validity datacheck

```bash
python3 assets/check_reference_coordinates.py /path/to/run/extended_transcripts GRCh38.fa.fai coordinate_issues.csv
```

The first argument is a GFF/GFF3 file or a directory of immediate GFF/GFF3 files; gzip inputs are supported. Directory mode rejects ambiguous plain/gzip or GFF/GFF3 copies with the same basename. The second argument is an explicit five-column FASTA index (`.fai`) or two-column sequence-name/length file, optionally gzipped. Reference names must be unique and lengths positive. Chromosome names must match **exactly**: `1` and `chr1` are different names, with no automatic aliasing. The reference must contain every sequence used by the checked annotation.

Every feature type, including genes and unfamiliar features, must satisfy `1 <= start <= end <= reference_length`. Unknown chromosome names, coordinates before position 1, reversed intervals and coordinates beyond sequence ends fail. The CSV contains one row per violating feature with all applicable issue codes, path, line, feature/transcript/gene identifiers when present, chromosome, strand, coordinates, reference length and explanation. A passing report has only its header. Console counts include files, total features, violating features and independent issue counts.

Exit codes are `0` for pass, `1` for unknown sequence names or coordinate-range violations, and `2` for malformed/unreadable/ambiguous input or execution errors. Noninteger coordinates, malformed GFF columns/attributes/score/strand/phase, invalid index entries, duplicate reference names and empty reference files are input errors. Empty/header-only GFFs are accepted; missing eligible content is the completeness check's responsibility. Headers are skipped, and parsing stops at embedded FASTA. This check does not validate hierarchy, sequence content, reference assembly identity or index byte offsets against an actual FASTA. Matching length/name tables must be supplied for the intended assembly. Existing reports remain untouched on input errors, and input overwrite protection includes symbolic and hard links.

Run focused tests with `python3 -m unittest discover -s tests -p test_check_reference_coordinates.py`.

## Added UTR overlap datacheck

```bash
python3 assets/check_cds_overlaps.py current_annotation_1.gff3 1_positive_extended.gff utr_overlap_comparison.csv \
  --leap-output-dir /path/to/original_leap_run \
  --extension-stats /path/to/leap_run/extended_transcripts_stats
```

This standalone standard-library check accepts plain or gzipped pipeline GFF3 and exact versioned transcript IDs. The third positional argument is a **CSV filename**. The required `--leap-output-dir` supplies the original run's `natural_overlaps/` and `fixed_extensions/` files. Missing files, multiple candidate files or incompatible selected candidate context are input errors. Repeated rows at the same splice site follow the writer's last-row rule. Raw `transcript_type` and UTR features reconstruct LEAP's full-length protein-coding eligibility (both 5′ and 3′ UTRs required). Other transcript types are non-protein-coding for these rules, even when their gene is protein-coding.

The overlap rules come from LEAP's clipping and terminal-safety stages, applied to the checked annotation with the original run context. Clipping protects eligible CDS and non-protein-coding exons from **all genes, including the source gene**. A full-length same-gene CDS contact is therefore assessed using LEAP's grouped annotation boundary and candidate support, rather than an extra per-isoform prohibition.

There is one deliberate stricter check: **within-gene CDS from partial raw protein-coding transcripts remain protected**, even when LEAP filtered those transcripts out. This extra check uses added UTR bases relative to the actual raw isoform, but separates bases already present within its original grouped annotation boundary from new read-driven growth. Existing grouped sequence is reported as `allowed_grouped_annotation`, even if the overall end has `reads` or mixed support. The grouped extent is verified using full-length multi-exon raw terminal exons with the same chromosome, strand, terminal kind, splice site and candidate gene-group membership, and is capped by the recorded original boundary when available. Partial, noncoding, single-exon and internal-exon records cannot enlarge this allowance. Its existing shared-internal-splice-site and purely `current_annotation` support allowances are retained. Both shared boundaries must be real splice sites; single-exon targets cannot qualify. Leftward extension uses the exon end, rightward extension the exon start, on either strand. Mixed annotation/read support does not grant this extra annotation allowance.

LEAP's two overlap stages are:

- **Clipping:** protect full-length coding models' CDS, single-exon CDS, and non-protein-coding exons. For single-CDS genes in different-gene natural-overlap pairs, apply LEAP's additional flank restrictions when the candidate's complete gene-ID field is that singleton gene. The flank check uses all such special genes' exons on the appropriate side of their CDS, falling back to the CDS-containing exon when that side has no exon.
- **Terminal safety:** test FE–LE and terminal–single-exon contacts. A contact is safe when the ordered candidate gene-ID groups share an ID, or any **zipped ordered** FE/LE gene pair occurs in the recorded natural-overlap file. Overlapping gene spans alone grant no allowance. FE–FE, LE–LE and ordinary internal-exon contacts are outside this stage.

Annotation-supported cross-gene extensions bypass clipping but still fail unsafe terminal contacts. The final writer's last BED row per chromosome/strand/terminal-kind/splice-site key supplies the candidate support and ordered gene group for LEAP clipping. Differing earlier rows are not merged or treated as ambiguous. Per-transcript stats cannot override a read-derived candidate to bypass clipping; stats enable the separate partial-CDS allowance. Other raw exon contacts are retained as **non-failing audit findings**, including unprotected protein-coding UTRs and ineligible partial models. The partial-CDS check is the only intentional additional overlap restriction.

The clipping guard includes the grouped original boundary recorded in the fixed-extension BED, which can differ from a shorter isoform's raw boundary. It uses inclusive comparisons, as `stop_cds_overlap.py` does. Terminal safety checks whole candidate exons, including their already annotated portions. Contacts use LEAP's positive-length BED overlap predicate: touching endpoints alone do not qualify. Reports use inclusive genomic coordinates; clipping findings can therefore include the grouped original boundary. The terminal tests use the annotation being checked, with original run metadata, so internally processed files can be evaluated without restoring reverted candidates. Use **processed stats** with processed annotations.

The CSV has one row per source/target, target-feature type and rule. `check_rule` explains each result; `intentional_partial_cds_check=true` identifies the extra partial-CDS rule; `gene_relationship`, `target_is_partial` and `target_missing_utr` identify the relationship and partial-CDS evidence. Separate `grouped_annotation_allowed_intervals`, `annotation_allowed_intervals`, `splice_site_allowed_intervals`, `leap_allowed_intervals` and `failing_intervals` retain the evidence and mixed outcomes. `grouped_annotation_boundaries` records strand-aware 5′/3′ limits, and `grouped_annotation_donors` names supporting eligible raw transcripts. A contact crossing the original grouped boundary retains its allowed portion and fails only for the remaining unexempted read-driven bases. Counts are merged bases within each row; different rules can overlap, so row/base totals are not unique transcript counts. Failures have status `new_cds_overlap`, `new_cross_gene_cds_overlap`, `new_cross_gene_exon_overlap`, `new_within_gene_exon_overlap` or `unsafe_terminal_overlap`. Allowances use `allowed_grouped_annotation`, `allowed_current_annotation`, `allowed_shared_splice_site`, `allowed_multiple_reasons` or `allowed_leap_overlap`. Unmatched extended UTR transcripts are `missing_from_raw`, counted separately without failing.

Stats accept a TSV file/directory, plain or gzipped. Gene, transcript location and extension lengths must match the checked annotation; stale or conflicting stats are rejected before the report is opened. Missing entries grant no within-gene annotation allowance. Only same-chromosome, same-strand contacts count. Self-comparisons are excluded. The parser checks IDs, locations, containing exons and intervals; duplicate intervals are deduplicated and overlapping distinct exons rejected. Unrelated transcript-only raw records are accepted. Reports cannot overwrite annotations, stats or required run-context files.

Exit codes are `0` for pass, `1` for biological failures and `2` for input/output errors. The console reports `leap_rule_failures` and `partial_cds_failures` separately, as well as evaluated UTR transcripts and report rows by outcome. Production extension rules and standalone internal-extension processing are unchanged; this check now also runs automatically in Nextflow. Inspectable examples, including all seven retained within-gene partial-CDS failures, are in [tests/data/leap_cross_gene](tests/data/leap_cross_gene/README.md). Earlier within-gene examples remain in [tests/data/cds_overlaps](tests/data/cds_overlaps/README.md). See [the LEAP overlap-rule contract](docs/leap_overlap_rules.md) for the source mapping and intentional exception. Run `python3 -m unittest discover -s tests` for regression and production-predicate parity checks.

## Internal-extension processing asset

```bash
python3 assets/process_internal_extensions.py /path/to/leap_run /path/to/new_processed_run
```

This standalone asset brings LEAP2GO's existing internal/external classification into LEAP without changing the Nextflow workflow. It restores internal terminal exons from raw annotation, saves original candidates separately, and writes strand-correct extension stats, including zero-extension rows for fully restored transcripts. Already blacklisted 5′ ends retain support/read-count values `blacklisted` and their `LEAP_blacklisted` tag; other restored ends use `none / 0`. Internal classification does not create new blacklist entries. Transcripts without CDS are excluded as before. It requires only Python's standard library and local `sort`, and leaves input files untouched.

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

Open the generated `LEAP-overlap-review.jbrowse` in JBrowse Desktop. The bundle includes indexed raw/extended annotations, CDS/UTR context, red/orange failure markers, green shared-site, teal annotation-supported, separate grouped-annotation allowances and blue LEAP audit markers, example sessions and a table of review locations. It uses public GRCh38 and requires `sort`, `bgzip` and `tabix`. Inputs may be files or directories; the output directory must be new.

For future datachecks, omit `--reports` and add marker GFF3 files with repeatable `--extra-track 'NAME=PATH'` arguments. The full workflow, marker format and reuse notes are in [the JBrowse review guide](docs/jbrowse_reviews.md).

## Author

Ryan Merritt

## License

See [LICENSE](/Users/rmerritt/my_app/LEAP/LICENSE).
