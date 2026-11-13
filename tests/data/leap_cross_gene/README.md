# LEAP-aligned cross-gene datacheck examples

Each directory contains raw/extended GFF3 and a minimal original-run `leap_run/` context. Use the context bundled with that example; chromosome/gene names and coordinates are synthetic except in the seven-partial regression.

```bash
# Fails a full-length cross-gene CDS guard, with remaining exon contacts reported as audit allowances.
python3 assets/check_cds_overlaps.py \
  tests/data/leap_cross_gene/full_length_cds/raw.gff3 \
  tests/data/leap_cross_gene/full_length_cds/extended.gff3 /tmp/full_length_cds.csv \
  --leap-output-dir tests/data/leap_cross_gene/full_length_cds/leap_run

# Preserves the seven real within-gene partial-CDS failures after internal processing.
python3 assets/check_cds_overlaps.py \
  tests/data/leap_cross_gene/seven_partial_cds_failures/raw.gff3 \
  tests/data/leap_cross_gene/seven_partial_cds_failures/extended.gff3 /tmp/partial_cds.csv \
  --leap-output-dir tests/data/leap_cross_gene/seven_partial_cds_failures/leap_run \
  --extension-stats tests/data/leap_cross_gene/seven_partial_cds_failures/stats.tsv
```

| Example | Exit | Rule |
| --- | --- | --- |
| full_length_cds | 1 | Protect eligible multi-exon CDS; audit other bases |
| noncoding_exon | 1 | Protect exons from non-protein-coding transcript types |
| single_exon | 1 | Single-exon CDS and terminal safety both apply independently |
| unknown_terminal_pair | 1 | FE–LE contact with no shared group or recorded pair |
| natural_terminal_pair | 0 | Same contact permitted by recorded ordered pair |
| annotation_unsafe_terminal | 1 | Annotation support bypasses clipping, never terminal safety |
| bed_touch | 0 | BED endpoints touch without positive-length terminal overlap |
| seven_partial_cds_failures | 1 | Seven unsafe read-driven hits on two partial raw CDS models |
| grouped_annotation_partial_cds | 0 | SCNN1A contact lies within the verified original grouped boundary |
| ddx41_duplicate_candidates | 0 | Coincident annotation/read records use the final row, matching the writer |
| shared_first_exon_different_genes | 0 | Raw ADORA3 and TMIGD3 models share their first exon despite different gene IDs |

The last fixture is extracted from the supplied `leap_correct` raw annotation and its internally processed output, with corrected stats and the original relevant fixed-candidate metadata. Its targets are `ENST00000534890.1` on chromosome 12 and `ENST00000585355.5` on chromosome 17, each missing a 3′ UTR. All seven failures must remain even though production LEAP excludes these partial CDS models from cross-gene protection.

`tests/test_leap_cross_gene_rules.py` also covers both strands, single-CDS-gene flank rules, grouped original boundaries, ordered zip versus Cartesian/undirected matching, shared groups, contacts not tested by LEAP, missing/ambiguous context, and clipping boundary inclusion. It executes production LEAP's clipping and terminal predicates for parity without changing those scripts or requiring pandas for the datacheck.

The extra within-gene CDS restriction applies only to **partial** raw protein-coding targets. Full-length CDS and non-coding exons from the source gene follow the same production LEAP clipping predicate as other genes. `intentional_partial_cds_check` marks the extra rule in CSV; the console separates partial-CDS failures from LEAP-rule failures. The raw grouped annotation boundary and fixed BED support control LEAP clipping; stats do not supply an additional bypass for eligible features.

The ADORA3/TMIGD3 fixture is extracted from the supplied raw and processed annotations, with original fixed-candidate metadata and corrected stats. Both records explicitly label the overlapping negative-strand exon `exon_number=1`. This case passes as an existing shared-terminal annotation contact; it is not treated as evidence that new first–last contacts may bypass natural-overlap checks.

The differential clipping matrix tests 144 combinations of strand, extension end, gene relationship, candidate support, protected feature and boundary contact against the actual `stop_cds_overlap.py` function. Terminal tests execute the actual `overlap_safety_checker.py`. The seven partial-CDS failures remain explicit regression expectations.

The `grouped_annotation_partial_cds` fixture is a minimal extraction of the newer human chromosome 12 SCNN1A case. Source `ENST00000396966.6` ends at 6374803 in raw annotation; eligible donor `ENST00000360168.8` shares its splice site and extends the raw group to 6375224. The selected read endpoint is 6375258. The 34-base contact with partial `ENST00000536788.1` at `12:6374804-6374837` is reported as `allowed_grouped_annotation`, with the donor and original boundary recorded. Removing that donor eliminates the allowance even if the saved candidate baseline is unchanged. The seven unsafe-read examples above continue to fail.

The `ddx41_duplicate_candidates` fixture retains two eligible DDX41 raw models, the extended source `ENST00000330503.12`, its stats and original fixed-candidate rows. At negative-strand LE splice site 177511927, annotation and read-labelled rows both span `177511569-177511927`. The last read-labelled row supplies support and original-boundary metadata. Reversing the rows selects annotation metadata instead. These are valid ordered candidates, rather than ambiguous input. Synthetic regressions confirm both strands and terminal ends, unsafe read-last outcomes, ordered gene groups and rejection of an incompatible last row without falling back to earlier entries.
