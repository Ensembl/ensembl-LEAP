# Build a JBrowse datacheck review

`assets/build_jbrowse_review.py` reproduces the LEAP review session: raw and extended transcripts, raw CDS and UTR segments, overlap failures, shared-splice-site allowances, and current-annotation allowances. The generated session format was checked with JBrowse Desktop 4.3.0.

## Build the overlap review

Requirements: Python 3 and `sort`, `bgzip`, and `tabix` on your command path. The latter two are supplied by HTSlib. No additional Python packages are required.

From the repository root:

```bash
python3 assets/build_jbrowse_review.py \
  --raw-annotation /path/to/leap_run/raw_annotation \
  --extended-annotation /path/to/leap_run/extended_transcripts \
  --reports /path/to/cds_overlap_reports \
  --output /path/to/new_jbrowse_review
```

Each annotation argument accepts one GFF3 file or a directory of `.gff`, `.gff3`, `.gff.gz`, or `.gff3.gz` files. A directory is read non-recursively. `--reports` accepts a single CSV or a directory of CSVs from `assets/check_cds_overlaps.py`; use reports generated from the same annotations. To apply the current-annotation exception, generate reports with the datacheck's `--extension-stats` option:

```bash
python3 assets/check_cds_overlaps.py raw.gff3 extended.gff3 report.csv \
  --leap-output-dir /path/to/original_leap_run \
  --extension-stats /path/to/leap_run/extended_transcripts_stats
```

The builder reads the resulting reason-specific intervals. It does not need the stats TSV files itself and cannot retroactively apply the exception to old reports. `--threads` sets compression threads (default 2).

The output directory must not already exist. The builder uses a temporary staging directory, publishes the bundle only after a successful build, and removes staging files on failure. Inputs are untouched. Large builds temporarily require enough disk space for uncompressed display copies and sorting, in addition to the compressed output.

The default reference is public human **GRCh38**, with chromosomes **1–22, X, Y**. Use `--assembly GRCm39` for mouse **mm39**, with chromosomes **1–19, X, Y**. The mouse config uses the UCSC two-bit reference and an explicit `chr`/unprefixed alias file; annotation coordinates and chromosome names are preserved. This script does not rename or lift over coordinates.

## Open and inspect

In JBrowse Desktop, choose **Open saved session (.jbrowse) file** and select `LEAP-overlap-review.jbrowse` from the output directory.

| Track | Purpose |
| --- | --- |
| Raw annotation | Original transcript, exon, CDS and UTR structure |
| Extended annotation | Extended structure, with transcript labels and LEAP tags |
| Raw CDS | Standalone coding segments for comparison |
| Raw UTR | Original UTR segments |
| Red: same-gene failure | LEAP clipping violation within the source gene, or the intentional partial-CDS failure |
| Orange: cross-gene failure | LEAP clipping or terminal-safety violation, including CDS hits |
| Green: shared splice site | An overlap portion allowed by the shared-site exception |
| Teal: current annotation | An overlap portion allowed because that extending end is annotation-supported |
| Blue: LEAP audit | Non-failing exon contacts outside protected regions, or terminal pairs allowed by LEAP |

Click markers for source/target transcript and gene IDs, rule, gene relationship, partial-target evidence, the intentional partial-CDS flag, splice sites and annotation-supported ends. Failure routing uses the source/target gene relationship, including cross-gene CDS hits. Annotation support never overrides an unsafe cross-gene terminal contact. Mixed comparisons have both allowed and failing markers. A green or teal portion does not make the entire comparison pass. Annotation provenance takes precedence for bases qualifying for both allowances, so they appear once in teal. Reports from before the annotation exception remain supported; their allowed intervals appear in green. Clipping intervals include the grouped original boundary; terminal contacts require positive-length BED overlap. Marker tracks reproduce the input CSV intervals; the builder does not run the datacheck or reinterpret its rules.

**File → Pre-configured sessions** contains one example for each nonempty marker track. If no markers exist, the session opens on the first extended annotation feature. The `review-loci.tsv` file lists overlap comparison locations you can paste into JBrowse's location box.

The generated `config.json` contains the same initial configuration as the `.jbrowse` file. `build-summary.json` records source file paths and counts; marker counts are individual intervals, not unique transcripts. `examples.json` records the selected example rows. JBrowse may normalize and save the session document after opening it.

These are **desktop** bundles with absolute local paths. Keep the bundle in place, or rebuild at its new location. The local data are BGZF-compressed and tabix-indexed. Internet access is required only for the public reference. To use a web-hosted JBrowse instance, the data must instead be served over HTTP with byte-range access and the local paths replaced with URLs.

## Reuse for another datacheck

The reusable pattern is:

1. Preserve raw and extended annotations as the comparison context.
2. Convert each issue or allowed region into a GFF3 feature with a descriptive `Name` and evidence attributes.
3. Sort by chromosome/start, compress with BGZF and index with tabix.
4. Add a separate feature track and an example session at the first issue.

The builder already implements steps 3–4. Use repeatable `--extra-track 'NAME=PATH'` arguments to add any datacheck's GFF3 markers. `PATH` may be a file or directory, including gzipped GFF3. Omit `--reports` if the review is unrelated to CDS/UTR overlaps:

```bash
python3 assets/build_jbrowse_review.py \
  --raw-annotation /path/to/raw_annotation \
  --extended-annotation /path/to/extended_transcripts \
  --extra-track 'Internal exon changes=/path/to/exon_changes.gff3' \
  --extra-track 'Shortened transcripts=/path/to/shortened_transcripts.gff3' \
  --output /path/to/new_structure_review
```

A minimal extra marker record has nine tab-separated fields:

```text
1	datacheck	sequence_feature	200	205	.	+	.	ID=issue_1;Name=Changed_exon;transcript_id=ENST_example;reason=internal_boundary_changed
```

GFF3 coordinates are **one-based and inclusive**. Convert a zero-based half-open interval `[start, end)` to GFF3 with `start + 1` and `end`. Give every independent marker a unique ID; percent-encode reserved characters in attributes. Extra marker records are retained as supplied, including evidence attributes and any valid parent/child hierarchy. They are not required to have `transcript_id` attributes.

The first two datachecks currently emit comparison CSVs, not these extra marker GFF3 files. Their CSV-to-marker conversion is a separate future addition. Do not pass those CSVs as overlap reports: their report format and meaning differ.

For annotation display copies, the builder removes gene rows and makes transcript records top-level. Transcript/exon/CDS/UTR coordinates are preserved. Child IDs are made unique, and transcript and LEAP metadata remain available. This display preparation does not change the original annotations.

## Validation

```bash
python3 -m unittest discover -s tests
```

Builder integration tests use small synthetic annotations and check real compression/indexing, coordinate preservation, empty marker tracks, mixed allowed/failing outcomes, custom tracks, gzip and directory input, output protection, and failed-build cleanup. These tests are skipped if `sort`, `bgzip`, or `tabix` is unavailable.

JBrowse references: [GFF3 tabix adapter](https://jbrowse.org/jb2/docs/config/gff3tabixadapter/), [reference adapter](https://jbrowse.org/jb2/docs/config/bgzipfastaadapter/), [default sessions](https://jbrowse.org/jb2/docs/config_guides/default_session/).
