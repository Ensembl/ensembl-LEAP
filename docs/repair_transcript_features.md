# Repair retained transcript features

`assets/repair_transcript_features.py` is an optional, standalone repair asset. It rebuilds retained models from the published raw annotation while preserving verified outward transcript endpoints. It is not called by Nextflow, manual-review preparation or the datachecks.

## Standard LEAP output

Run once for each corresponding annotation/stats pair:

```bash
python3 assets/repair_transcript_features.py \
  /path/to/run/raw_annotation/current_annotation_21.gff3 \
  /path/to/run/extended_transcripts/21_positive_extended.gff \
  /path/to/new_repair_directory \
  --extension-stats /path/to/run/extended_transcripts_stats/21_positive_stats.tsv \
  --leap-output-dir /path/to/run
```

Inputs may be plain or gzipped. The output directory must not already exist, including as a symbolic link. Source files are never overwritten. Python's standard library and the existing sibling datacheck assets are sufficient; no database or external sorting command is needed.

## Processed external output

Use the processed GFF and its corresponding processed stats, with the original LEAP run supplying candidate context:

```bash
python3 assets/repair_transcript_features.py \
  /path/to/run/raw_annotation/current_annotation_21.gff3 \
  /path/to/processed/external_extensions/gff/21_positive_extended.gff \
  /path/to/new_processed_repair_directory \
  --extension-stats /path/to/processed/external_extensions/stats/21_positive_stats.tsv \
  --leap-output-dir /path/to/run \
  --stage processed
```

Use the actual processed filenames if they differ. Processed output remains gene-free and retains its existing, unambiguous `transcript_class`. Internal classification does not imply blacklisting.

## Repair rules

- Match exact versioned transcript IDs and chromosome. Preserve raw IDs and Parent relationships, including distinct GFF `_PAR_Y` IDs and `gene_id` attributes.
- The retained set consists of models represented by any transcript-associated record in the supplied output. It must match the supplied stats membership exactly. Do not restore wholly omitted transcripts, add new models or change eligibility.
- Restore all raw transcript metadata and children: exons, every CDS/codon segment and phase, UTRs, IDs, parent links, exon numbers, source, score and unfamiliar feature types. Remove unexpected output features. Raw structure must itself be valid for every retained model.
- Keep the existing transcript row's outer coordinates. Never choose a new extension or move an endpoint inward. Rebuild the lowest-coordinate exon start, highest-coordinate exon end and corresponding terminal UTR outer boundaries to those coordinates. All splice-facing boundaries stay raw on both strands.
- When the transcript row is missing but children remain, use stats coordinates only if the recorded location and strand-aware extension lengths agree with raw. Candidate verification is still required for retained growth.
- Verify each extended end against the original run's `fixed_extensions/` context. Follow the writer's **last row per splice site** before validating gene membership and support. Standard output must match that selected endpoint exactly. Processed output may retain an outward endpoint contained within the selected candidate. Do not fall back to an earlier candidate. Missing or ambiguous context stops repair. Context is optional only when no extension survives.
- Recompute stats locations, strand-aware lengths and transcript `leap_tag`. Retained extensions use the selected candidate's support/count; equivalent repeated support/count entries are preserved. Reverted or unextended ends use `none / 0`. An existing blacklisted 5′ end stays `blacklisted / blacklisted` with `LEAP_blacklisted` and must remain unextended. Preserve a matching existing tag ordering; generated child tags are removed.
- Standard output includes genes containing the retained repaired transcripts. Preserve raw gene IDs/metadata where available, infer a gene record from raw transcript metadata if absent, and recalculate bounds across retained models. Report gene changes separately.

Shortened transcript rows, duplicate transcript rows, conflicting classifications, mismatched membership/gene information and unverifiable endpoints stop the entire repair. Malformed syntax, invalid selected raw models and invalid required context also stop processing. Unrelated raw models excluded from the supplied output are not rebuilt or checked for structural correctness.

## Outputs and exit codes

On success:

- `repaired.gff3`: repaired annotation, deterministically coordinate sorted, with a GFF3 header. Record/attribute order and original headers are not preserved.
- `repaired_stats.tsv`: corresponding corrected stats in the standard LEAP schema.
- `changes.csv`: one row per retained transcript, with unchanged/repaired status, original output and repaired bounds, full added/removed records, metadata changes, stats changes and generated attributes.
- `gene_changes.json`: before/after gene rows.
- `validation.csv`: results of complete transcript preservation, stats/GFF consistency and structural-integrity checks.
- `summary.json`: counts, input paths, stage, outcome and diagnostics.

Exit `0` means repaired copies passed those three checks. Exit `1` means at least one model could not be repaired unambiguously; exit `2` means an input or execution error. Failed attempts retain diagnostics when a new report directory can be created, but publish **no repaired GFF or stats**, even if some models were individually repairable. An existing output directory is rejected without modification.

This asset restores features; it does not resolve unsafe overlaps, validate reference identity, rerun extension selection or prove completeness. Run the full datacheck suite on the repaired copies before using them for release. Preserve the original reports so the initial defects remain visible.

## Tests

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s tests -p test_repair_transcript_features.py -v
```

Fixtures cover healthy models, missing/altered features and metadata, missing transcript rows, both strands, single exons, multipart codons, gzip, stats/tags, blacklisting, processed endpoints, candidate order, raw errors, membership, PAR links and overwrite protection.
