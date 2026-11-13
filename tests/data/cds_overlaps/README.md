# Original UTR overlap test datasets

These early fixtures retain the original within-gene cases. Their cross-gene protein-coding targets do not contain both UTR types, so the LEAP-aligned check now reports those contacts as `allowed_leap_overlap` for ineligible partial models. `expected.json` is authoritative. For complete cross-gene rule examples with original-run context, use [../leap_cross_gene/README.md](../leap_cross_gene/README.md).

The CLI now requires `--leap-output-dir`. The automated runner supplies minimal synthetic context using `tests/leap_test_context.py`. For immediately runnable CLI commands use the newer examples linked above.

Each folder contains `raw.gff3` and `extended.gff3`. Coordinates are deliberately small so you can inspect them directly. `source` is the extending transcript in `geneA`. The target transcripts are present only in raw, which also checks that filtering a target out of extended does not hide it.

Run the examples through the automated runner, which supplies each minimal original-run context:

```bash
python3 -m unittest discover -s tests
```

For standalone CLI commands with saved context, use the [LEAP rule examples](../leap_cross_gene/README.md).

For the first example, the source's raw first exon is `200–250` and extends to `150–250`. Its added UTR is `150–199`. The target's first exon is `140–250`, containing CDS `150–180`. Both exons end at the real internal splice boundary `250`, so the added overlap is allowed. In the failing example, the target exon ends at `240` instead.

| Dataset | Expected | Reason |
| --- | --- | --- |
| `same_gene_shared_site` | PASS, allowed row | Matching real internal exon end |
| `same_gene_wrong_site` | FAIL | Different exon ends |
| `same_gene_mixed` | FAIL | Allowed `190–199`, forbidden `160–170` |
| `single_exon_no_shared_site` | FAIL | Target has no internal splice site |
| `target_outer_boundary_not_splice_site` | FAIL | Matching end is the target's outer boundary |
| `same_gene_existing_overlap_grows` | FAIL | Added CDS overlap still requires a shared site |
| `cross_gene_intron` | PASS, no rows | Added bases avoid all other-gene exons |
| `cross_gene_cds` | PASS, audit | Partial coding target is ineligible for LEAP clipping |
| `cross_gene_utr_only` | PASS, audit | No eligible clipping or terminal guard |
| `existing_overlap_unchanged` | PASS, no rows | All shared bases were already in raw UTR |
| `existing_overlap_grows` | PASS, audit | Growth contacts an ineligible partial target |
| `opposite_strand` | PASS, no rows | Opposite-strand targets are excluded |
| `one_base_overlap` | PASS, audit | Partial target; one base is still reported |
| `adjacent_no_overlap` | PASS, no rows | Adjacent intervals share no bases |
| `other_isoform_exonic` | PASS, audit | These partial isoforms are ineligible clipping targets |
| `self_overlap_ignored` | PASS, no rows | A transcript is not compared with itself |

The `_mirrored` cases reflect coordinates around 1000 and flip strand, exercising rightward extensions on the negative strand. The `_strand_flipped` cases additionally cover leftward negative-strand and rightward positive-strand extensions. All variants retain their original expected outcome. `expected.json` records the expected exit code and report status for every dataset.

A header-only CSV means no added forbidden overlap and no shared-site exemption needed reporting. `allowed_shared_splice_site` rows also pass. Exit `2` means invalid input or an I/O error, rather than a biological failure.

## Current-annotation support examples

The `annotation_*`, `reads_same_gene_still_fails` and `mixed_annotation_reads_not_exempt` folders also include `stats.tsv`. The automated runner supplies these through `--extension-stats` alongside the required run context. In the same-gene cases, pure annotation support passes despite differing splice sites; omitting those stats fails because the provenance cannot be established from GFF3 alone.

| Dataset | Expected with its stats | Reason |
| --- | --- | --- |
| `annotation_same_gene_allowed` | PASS | Annotation support exempts the new CDS overlap |
| `annotation_cross_gene_cds_allowed` | PASS | Exception also covers other genes |
| `annotation_cross_gene_utr_allowed` | PASS | Partial target is audited; unsafe terminal contacts would still fail |
| `annotation_left_negative_allowed` | PASS | Leftward negative-strand extension uses 3′ support |
| `annotation_right_negative_allowed` | PASS | Rightward negative-strand extension uses 5′ support |
| `annotation_right_positive_allowed` | PASS | Rightward positive-strand extension uses 3′ support |
| `reads_same_gene_still_fails` | FAIL | Read support grants no additional allowance |
| `mixed_annotation_reads_not_exempt` | FAIL | Support is not purely current annotation |
| `annotation_and_reads_opposite_ends` | FAIL | Annotation-supported 5′ end passes; read-supported 3′ end fails |
| `annotation_and_shared_site_allowed` | PASS | One end passes via annotation, the other via shared splice site |

Support categories, rather than read counts, control this exception. The automated dataset runner uses each folder's `stats.tsv` when present. Stale gene/location/extension lengths, malformed stats, and duplicate chromosome/transcript entries are separately tested as input errors.
