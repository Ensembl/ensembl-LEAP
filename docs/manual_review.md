# Prepare annotation for manual review

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

## Reference and launch settings

Provide your server's matching FASTA (including gzip) or 2bit file to LEAP:

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

`--review_genome_id` and `--review_species` explicitly set JBrowser's assembly ID
and publication species directory. When omitted, a reference filename containing
`GRCh38`/`hg38` supplies `hg38` and `homo_sapiens`; `GRCm39`/`mm39` supplies `mm39`
and `mus_musculus`. These are filename-based defaults, not assembly verification.
For other references, or absent reference paths, unknown values are marked
`EDIT_ME` in the samplesheet. Exact reference sequence names must match annotation
names; no chromosome renaming or coordinate conversion is performed.

The reference remains at its original server path; LEAP does not copy the genome
into the review bundle. If available, the datacheck length input is converted to
the two-column `reference.chrom.sizes` format JBrowser expects. A `.fa.gz.fai` is
an ordinary index input here. It is distinct from the FASTA/2bit passed through
`--review_reference`.

The bundle contains:

```text
manual_review/
  annotations/{raw,leaped,external,internal}/{positive,negative,unsplit}.gff3
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

The two sheets follow the current HAVANA JBrowser schema: annotation/reference
rows in `samplesheet.csv`, and tags/descriptions separately in
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

Three new modules split raw annotation, process each LEAP chromosome/strand, then
merge the views and write the launch bundle. They consume declared workflow
channels, so they do not race asynchronous published copies. Missing, duplicated
or unexpected partitions fail clearly instead of silently producing an incomplete
review. Sorting errors stop the build. Bundle publication occurs only after a
complete build. The final datacheck gate waits for review preparation, preserving
review files even when a datacheck reports a biological failure or input error.

To prepare the same bundle from an existing LEAP run, without a database or BAM
rerun:

```bash
python3 bin/prepare_manual_review.py run /server/leap_results /server/new_review \
  --chromosomes 1 2 21 X Y \
  --reference-genome /server/reference/genome.fa \
  --reference-lengths /server/reference/genome.fa.fai \
  --genome-id hg38 --species homo_sapiens \
  --jbrowser-pipeline /server/software/havana-jbrowser
```

The standalone command also accepts `--public-base-url`, `--ftp-target-dir`,
`--hub-email`, `--ftp-project`, `--publication public|offline`, and
`--strand positive|negative|both`. It defaults to both strands and public
publication. Plain/gzipped pipeline inputs are supported; ambiguous copies are
errors. Choose a new output directory. Input files are never overwritten.

Requirements are Python 3.10+ and local `sort`; no database or browser installation
is needed to prepare the bundle. Building all views needs temporary disk space
for partition copies and whole-annotation sorting, in addition to the 12 final
files. The generated job requires your normal JBrowser/Nextflow environment.
