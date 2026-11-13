# Process internal extensions within LEAP

`assets/process_internal_extensions.py` is a standalone LEAP asset based on LEAP2GO's `external_extensions.py`. It prepares for future Nextflow integration; it is not currently called by the LEAP workflow. LEAP2GO itself is unchanged.

```bash
python3 assets/process_internal_extensions.py /path/to/leap_run /path/to/new_processed_run
```

Only Python's standard library and a local `sort` command are required. The output directory must be new. Input GFF/GFF3 and stats files may be plain text or gzipped. Input names follow LEAP's conventions: `current_annotation_<chromosome>.gff3`, `<chromosome>_<positive|negative>_extended.gff`, and `<chromosome>_<positive|negative>_stats.tsv`.

## Classification and restoration

The classification preserves LEAP2GO's existing policy. Within each chromosome/strand input, collect the earliest-starting and latest-ending CDS segments per gene, keeping the first encountered segment on ties. A transcript matches an extreme if any of its CDS segments shares **either** the start or end of that extreme segment. This is a boundary-OR test, not a requirement for identical CDS intervals.

Transcripts matching both extremes are external. Otherwise their left and/or right terminal exon is internal. Terminal exons are selected with transcript-direction `exon_number` values, including reversed numbering on the negative strand. The processor substitutes the corresponding raw terminal-exon records and their associated features, and recalculates the transcript bounds from the resulting exons. Other exon coordinates are preserved. Features are sorted locally with checked return codes; no cluster-specific sorter path is required.

Transcripts without CDS are intentionally omitted, as in LEAP2GO. Gene rows are omitted from these transcript-oriented outputs. Full-version transcript IDs remain in GFF and stats. Classification text lists retain LEAP2GO's unversioned-ID convention; the additional classification TSV preserves full IDs and records dropped transcripts explicitly.

## Output layout

```text
new_processed_run/
  processing.log
  processing-summary.json
  external_extensions/
    gff/                 # Processed transcript annotation
    stats/               # Corrected per-end extension statistics
    ids/                 # Classification lists and full-ID classification TSVs
    potential_internal/  # Original extended candidates before restoration
```

The build is staged and published only after all input files succeed. Missing raw models, mismatched exon numbers, stale input stats or sorting failures stop the build and leave the original files untouched.

## Corrected statistics

Every retained transcript has a stats row, including transcripts restored at both ends. Transcript locations and extension lengths are calculated from the actual processed exon boundaries relative to raw annotation.

A restored end has zero extension, support `none`, and read count `0`, except an already blacklisted 5′ end: its support and read-count fields remain `blacklisted`, consistently with its retained `LEAP_blacklisted` tag. Blacklisting and internal/external classification are independent; processing does not blacklist additional transcripts. The left end maps to 5′ on `+` and 3′ on `-`; the right end maps the other way. Original support/read counts remain intact for surviving extensions. Extension tags (`LEAP_5` / `LEAP_3`) are removed for ends that no longer extend. `transcript_class` identifies the internal/external classification in the processed transcript records.

For a single-exon transcript, substituting its raw terminal exon restores that entire exon. Both transcript bounds and stats then follow the resulting exon, even if only one side triggered classification; there is no independent second terminal exon to retain.

## Datacheck and JBrowse follow-up

Run the CDS/UTR datacheck using processed GFFs and **processed** stats, while keeping the original raw annotation as the reference:

```bash
python3 assets/check_cds_overlaps.py \
  /path/to/leap_run/raw_annotation/current_annotation_1.gff3 \
  /path/to/new_processed_run/external_extensions/gff/1_positive_extended.gff \
  /path/to/new_overlap_report.csv \
  --leap-output-dir /path/to/leap_run \
  --extension-stats /path/to/new_processed_run/external_extensions/stats/1_positive_stats.tsv
```

Build a JBrowse review with the processed GFF directory and new overlap-report directory using `assets/build_jbrowse_review.py`. The processed transcript's `transcript_class` remains visible in feature details. See [the JBrowse review guide](jbrowse_reviews.md) for adding other tracks.

For before/after comparisons, compare **failing base intervals** for each directed source/target pair and target-feature type. Keep resolved, remaining and newly failing portions separate. A source transcript disappearing from processed annotation must be reported as dropped, not as a resolved overlap. Original overlaps allowed by the datacheck, shared-splice-site allowances and annotation-supported allowances remain governed by the same datacheck rules.

Small fixtures are in `tests/data/internal_extensions/`; the regression tests verify both strands, all restoration classes, unchanged middle exons/CDS, stats alignment, potential-internal preservation, input protection and sorting-error cleanup.
