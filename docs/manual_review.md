# Prepare annotations, read evidence and transcript tables for manual review

LEAP automatically prepares `manual_review/` after the final annotations and stats
are available. Set `--prepare_manual_review false` to disable this preparation.
It creates review copies; it does not replace LEAP's original outputs.

## Tracks

Each group has `positive.gff3`, `negative.gff3` and `unsplit.gff3`, combined across
all selected chromosomes. Unsplit includes both strands; raw unstranded records,
if any, appear only in its unsplit view. Records are coordinate-sorted for browser
conversion; coordinates, original IDs, parent relationships and metadata remain
intact.

| Directory under `annotations/` | Contents |
| --- | --- |
| `raw/` | Published raw annotation, before LEAP extension |
| `leaped/` | Original LEAP output, including internal extension candidates |
| `external/` | Existing internal-extension processor's output, with internal terminal ends restored to raw and retained external ends preserved |
| `internal/` | Original LEAP models with actual growth at an internal terminal end: what those models would look like if internal extensions were retained |

A transcript can appear in both external and internal views, with different
terminal boundaries. A model may have an internal extension at one end and an
external extension at the other. The internal view keeps the original whole
model, rather than only the newly added bases, and adds its classification as
`transcript_class` on the transcript row.

Classification uses the existing CDS-extreme policy described in
[internal-extension processing](internal_extension_processing.md). It identifies
internal terminal exons. The internal track requires a positive original
extension length at at least one such end, using strand-aware 5′/3′ stats.
Classification tables additionally retain internally classified models with no
growth; those models do not appear in the internal-extension track.
Models without CDS remain in raw/leaped views but are omitted from external and
internal views, following the processor's existing policy. Dropped models are
listed in the classification tables. Raw/leaped copies retain gene rows;
external/internal copies are transcript-oriented, as in the existing processor.

All 12 GFF3 files are produced, including header-only files for empty views.
Empty views remain in `track_manifest.tsv` but are omitted from the JBrowser
samplesheet because they cannot produce an annotation track. No empty view is
replaced with another strand's annotation.

## Read evidence

Automatic review preparation also uses LEAP's existing chromosome BAMs to build
`reads/positive.bam`, `reads/negative.bam` and `reads/unsplit.bam`, each with a
`.bam.csi` index. These views contain **all mapped input alignment records within
the run's selected chromosome/strand scope**. They are not restricted to reads
supporting a selected extension, or divided into external/internal support sets.

Positive/negative orientation follows BAM flag `0x10`, matching LEAP's existing
read-to-BED conversion. Unmapped records are excluded. Secondary, supplementary,
duplicate-marked, QC-failed and low-MAPQ mapped records remain; no new quality or
support filtering is applied. Repeated read names and identical alignment records
are preserved. Counts are **alignment records**, not unique molecules or LEAP's
per-extension support counts. SAM optional-tag order and processing-history
header lines may change during merging; alignment content, tag values and
read-group identities are retained.

The BAMs are coordinate-sorted and CSI-indexed, including empty views. Input
reference dictionaries/read groups must be compatible, and there must be one
chromosome BAM/index pair per expected chromosome. Missing, ambiguous, invalid or
incorrectly scoped inputs stop preparation. Source paths, selected-record counts
and output counts are retained in `reads/read-summary.json` and the bundle
summary. Empty views remain on disk and in the manifest but are omitted from
the samplesheets.

Nonempty read views use JBrowser's `type=bam` alignment route, with the matching
reference in `reference_file` and `data_category=alignment`. JBrowser validates
and converts them to indexed CRAM tracks. Both public and offline launchers
include the read rows; the offline launcher forces them offline too. If editing
the reference after bundle generation, update its reference row **and** the BAM
rows' `reference_file` entries in the shared samplesheet. The launchers check
these reference paths before starting JBrowser.

## CAGE and poly(A) evidence

Automatic review preparation includes the original `--cage` and `--polya` inputs,
with positive, negative and unsplit views in `sites/cage/` and `sites/polya/`.
These contain all supplied sites in the processed chromosome/strand scope, rather
than only sites selected to support an extension. Unstranded sites appear only
in the unsplit view. Original records and signal values are retained in each
source's `source_records.tsv`; source paths, counts, display format and coordinate
policy are recorded in `sites/site-summary.json` and `review-summary.json`.

BED inputs retain their zero-based starts; GFF/GFF3 starts are converted to BED
once. A leading `chr` prefix is removed, matching LEAP's evidence preparation.
No liftover or other chromosome alias inference is performed. Sites outside the
selected scope are counted and omitted. When reference lengths are provided,
selected sites must fit the reference. Malformed inputs or sorting errors stop
preparation without publishing a partial bundle.

Native BED9 CAGE features retain their thick bounds and colours. Display scores
are limited to integers from 0 to 1000. Each CAGE feature also carries an
`originalCageScore` field with the exact original score text, including values
above 1000, fractional scores and `.` for missing scores. This field is retained
in bigBed for JBrowse feature details; it is separate from the display score.
Original values also remain in the source record table. Eleven-column PolyA Atlas inputs use JBrowser's `polyAtlas` preset,
preserving fractional TPM, sample support, protocol counts, annotations and
motifs. Other supported BED/GFF evidence is normalized to BED6, retaining its
original fields in the source table.

All six BED files are produced, including empty views. Nonempty views become
`type=regions`, `data_category=regions` samplesheet rows using the appropriate
`cageBed9` (BED9+1), `cageBed6` (BED6+1), `polyAtlas`, `bed9` or `bed6` preset.
The bundle supplies the CAGE format definitions in `region_presets/` and an
`install_review_region_presets.py` helper. Both generated launchers install missing
definitions into the chosen JBrowser checkout before starting it. Existing
identical definitions are accepted; conflicting definitions are never overwritten
and stop the launcher with an explicit diagnostic. The checkout must be writable
for the first installation, or an administrator must copy the supplied definitions
into its `assets/region_presets/` directory beforehand.
Empty views remain in the manifest but are
omitted from the samplesheets. Public and offline launchers include these tracks;
the offline sheet sets every evidence row to `publication=offline` as usual.
Preparing these views does not alter extension selection or the original files.

## Transcript review spreadsheet

Every bundle includes **`spreadsheets/transcripts.csv`**, one combined spreadsheet
across the selected chromosomes and strands. It has one row per original LEAP
output transcript, including zero-extension models and models subsequently
excluded from external output because they have no CDS. Columns include exact
versioned transcript/gene IDs, chromosome, strand, 5′/3′
extension lengths, evidence and read counts, internal/external classification,
and the `jbrowse_label` (`chromosome:start-end`) for finding the original LEAP model.

This uses the same row builder as `impact/transcripts.csv` but omits the audit
columns `raw_start`, `raw_end`, `leap_start`, `leap_end`, `external_start` and
`external_end` from the shared spreadsheet. Those coordinates remain in the
detailed impact CSV; every other field and transcript row is retained, including
the JBrowse search label. It is assembled from declared annotation
and stats partitions without waiting for published impact files. Review and
impact reporting therefore stay independent before the final datacheck gate.

The table is automatically registered in **both** samplesheets as:

- `sample=LEAP_transcripts`, `type=spreadsheet`;
- `file=<published manual_review>/spreadsheets/transcripts.csv`;
- `data_category=leap_review`, with the bundle's FTP project and genome-ID namespace;
- the selected publication mode; the offline sheet always uses `publication=offline`.

JBrowser copies the local CSV unchanged. It is a downloadable spreadsheet, not
an annotation track, and needs no conversion, index or reference of its own.
After JBrowser preparation it is stored at
`<JBrowser output>/spreadsheet/LEAP_transcripts/transcripts.csv`. Its public
publishing script places it at
`spreadsheet/hg38/leap_review/YYYYMMDD_N/transcripts.csv` for human (`mm39` for
mouse), independently versioned from annotation/alignment data. The offline
bundle contains the same CSV alongside the prepared review data. Tags and the
description come from `sample_metadata.csv` as usual.

Header-only spreadsheets remain available when no LEAP models are present.
Missing, malformed or inconsistent per-partition spreadsheets stop bundle
publication. Generating the bundle and launch scripts does not upload anything;
JBrowser submission and its FTP publishing script remain manual steps.

## Reference and launch settings

Provide your server's matching FASTA (including gzip) or 2bit file to LEAP:

For an editable Slurm launcher with the paths and settings grouped at the top,
copy [`examples/run_leap.sbatch`](../examples/run_leap.sbatch), edit it, then submit
with `sbatch /path/to/run_leap.sbatch`. Set `RESUME=true` there to resume a run.
Its relative input/reference paths are resolved from `LEAP_DIR`.

```bash
nextflow run main.nf -profile hpc_slurm \
  --db_secrets /private/path/db_secrets.json \
  --review_reference /server/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
  --datacheck_reference_lengths /server/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz.fai \
  --review_jbrowser_pipeline /server/software/havana-jbrowser \
  --review_ftp_project LEAP_review \
  --review_public_base_url https://ftp.example.org/root \
  --review_ftp_target_dir /server/ftp/root \
  --review_hub_email reviewer@example.org \
  --reads /server/reads.bam --reads_index /server/reads.bam.csi \
  --cage /server/cage.bed --polya /server/polya.bed \
  --output_directory /server/leap_results
```

`--review_genome_id` sets JBrowser's assembly ID. `--review_species` sets the FTP
directory namespace: JBrowser calls this column `species`, but the online FTP
layout uses genome IDs such as `hg38`, not biological species names. Its default
is the resolved genome ID, so both columns contain `hg38` for human or `mm39` for
mouse. Explicit `--review_species` / standalone `--species` overrides are retained.
When the assembly ID is omitted, a reference filename containing `GRCh38`/`hg38`
supplies `hg38`; `GRCm39`/`mm39` supplies `mm39`. These are filename-based defaults,
not assembly verification.
For other references, or absent reference paths, unknown values are marked
`EDIT_ME` in the samplesheet. Exact reference sequence names must match annotation
names; no chromosome renaming or coordinate conversion is performed.

For an existing bundle using `homo_sapiens`, change the **`species` column on every
row** in `samplesheet.csv` to `hg38`, keeping `genome_id=hg38`. This includes the
reference, annotations and read-evidence rows. The offline launcher refreshes its
sheet from the shared sheet automatically. If JBrowser preparation already ran,
rerun that preparation with the corrected samplesheet before using the generated
FTP copy scripts; changing the sheet does not update an existing browser manifest.
No rerun of LEAP's annotation or extension processing is required.

The reference remains at its original server path; LEAP does not copy the genome
into the review bundle. If available, the datacheck length input is converted to
the two-column `reference.chrom.sizes` format JBrowser expects. A `.fa.gz.fai` is
an ordinary index input here. It is distinct from the FASTA/2bit passed through
`--review_reference`.

The bundle contains:

```text
manual_review/
  annotations/{raw,leaped,external,internal}/{positive,negative,unsplit}.gff3
  reads/{positive,negative,unsplit}.bam
  reads/{positive,negative,unsplit}.bam.csi
  reads/read-summary.json
  sites/{cage,polya}/{positive,negative,unsplit}.bed
  sites/{cage,polya}/source_records.tsv
  sites/site-summary.json
  region_presets/cageBed*.tsv              # Browser-visible original CAGE signal
  install_review_region_presets.py
  spreadsheets/transcripts.csv             # Combined transcript review table
  processing/<chromosome>_<strand>/stats/   # Corrected external-output stats
  processing/<chromosome>_<strand>/ids/     # Classification and dropped-model lists
  reference.chrom.sizes                   # When reference lengths are supplied
  samplesheet.csv
  samplesheet_offline.csv
  sample_metadata.csv
  track_manifest.tsv
  review-summary.json
  run_jbrowser.sbatch
  run_jbrowser_offline.sbatch
  README.txt
```

The two sheets follow the current HAVANA JBrowser schema: annotation, reference,
alignment, interval-evidence and spreadsheet rows in `samplesheet.csv`, and tags/descriptions separately in
`sample_metadata.csv`. Publication defaults to **public**, as requested. The FTP
project is per sample; the public URL and filesystem directory are FTP **roots**,
not project-specific subdirectories.

Before submitting, review `run_jbrowser.sbatch`: confirm the JBrowser checkout,
FTP roots, hub email, profile, and any site-specific account/partition/module
settings. It defaults to `slurm,singularity`, one CPU, 4 GB and 48 hours for the
Nextflow launcher. Individual JBrowser jobs use that pipeline's resource settings.
The script refuses to launch with unresolved samplesheet placeholders or missing
input files. It records the published review directory explicitly because Slurm
executes a spooled copy of the script.

```bash
sbatch /server/leap_results/manual_review/run_jbrowser.sbatch
# To resume JBrowser preparation:
sbatch /server/leap_results/manual_review/run_jbrowser.sbatch -resume
# Or choose private/offline preparation:
sbatch /server/leap_results/manual_review/run_jbrowser_offline.sbatch
```

The offline script requires no FTP roots or hub email. Edit reference paths and
assembly settings in the shared `samplesheet.csv`; on each launch, the offline
script regenerates `samplesheet_offline.csv`, preserving those edits and forcing
**every** row, including the reference, to `publication=offline`. The public
samplesheet remains unchanged. Both scripts share `sample_metadata.csv`.
Configure the JBrowser checkout/profile and Slurm settings in the script you
choose. Offline output, work, report and trace paths have an `_offline` suffix,
so public and offline preparation have separate locations. Add `-resume` to
either submission command when resuming that mode.

LEAP only writes these files: it does not submit the job or publish anything to
FTP. The generated job prepares JBrowser output and its FTP copy scripts. Review
those outputs before executing the copy scripts to publish. Annotation indexing
and bigBed/browser-session generation are handled by JBrowser, not this stage.
Paths in the sheets point to the published bundle. Update them and the script's
`REVIEW_DIR` in both scripts if moving the bundle.

## Workflow and existing runs

Five modules split raw annotation, process each LEAP chromosome/strand, prepare
read alignments and CAGE/poly(A) intervals, then merge the views and write the launch bundle. They consume declared workflow
channels, so they do not race asynchronous published copies. Missing, duplicated
or unexpected partitions fail clearly instead of silently producing an incomplete
review. Sorting errors stop the build. Bundle publication occurs only after a
complete build, including all BAMs and indexes. The final datacheck gate waits for review preparation, preserving
review files even when a datacheck reports a biological failure or input error.

To prepare the same bundle from an existing LEAP run, without a database or BAM
rerun:

```bash
python3 bin/prepare_manual_review.py run /server/leap_results /server/new_review \
  --chromosomes 1 2 21 X Y \
  --reference-genome /server/reference/genome.fa \
  --reference-lengths /server/reference/genome.fa.fai \
  --reads /server/reads.bam --reads-index /server/reads.bam.csi \
  --cage /server/CAGE.bed --polya /server/polyA.bed \
  --genome-id hg38 --species hg38 \
  --jbrowser-pipeline /server/software/havana-jbrowser
```

The standalone command also accepts `--public-base-url`, `--ftp-target-dir`,
`--hub-email`, `--ftp-project`, `--publication public|offline`, and
`--strand positive|negative|both`. It defaults to both strands and public
publication. Plain/gzipped pipeline inputs are supported; ambiguous copies are
errors. Choose a new output directory. Input files are never overwritten.

For standalone builds, `--reads` and `--reads-index` must be supplied together;
input indexes may be BAI or CSI, and are used explicitly. `--read-threads` controls
Samtools threads (default 1). The original BAM can contain other chromosomes;
only the selected chromosomes/strands enter the review views. Omit both read
options to retain the existing annotation-only command; the summary explicitly
records `read_evidence.status=not_supplied`. Automatic pipeline builds always
require read inputs when review preparation is enabled.

Standalone builds also accept optional `--cage` and `--polya` inputs, independently.
BED, GFF and GFF3 inputs may be gzipped. Omit both to keep the previous standalone
behaviour; the summary explicitly records `site_evidence.status=not_supplied`.
Automatic builds always stage both existing LEAP evidence inputs as declared
workflow inputs; no additional pipeline parameter is needed.

Requirements are Python 3.10+, local `sort`, and Samtools when including reads;
no database or browser installation is needed to prepare the bundle. The read
module declares Conda Python 3.11 and Samtools 1.18; deployments disabling Conda
must supply those tools. Building all views needs temporary disk space for
partition copies, sorting, BAM merging and the three delivered BAM/index pairs,
plus publication copies and sorted site-evidence views/source-record tables. Large input BAMs can substantially increase disk usage
and preparation time. The generated job requires your normal JBrowser/Nextflow
environment. No additional read option is needed for automatic LEAP runs:
`--reads` and `--reads_index` already provide the evidence.

## Impact statistics

The sibling `impact/` output is generated automatically even when manual-review preparation is disabled. Open `impact/index.html` for original/external/internal counts, read-supported sequence additions, length distributions, gene/transcript CSV downloads and datacheck outcomes. Its external/internal views use the same processor as this review bundle. See [impact statistics](impact_dashboard.md) for counting rules and limitations.
