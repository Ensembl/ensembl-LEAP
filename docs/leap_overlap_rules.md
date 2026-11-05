# LEAP overlap-rule contract

`assets/check_cds_overlaps.py` validates the overlap decisions in the current LEAP code. It does not add a blanket prohibition on entering another gene's exon or use overlapping gene spans as a natural-overlap whitelist.

The sole deliberate extra overlap restriction is the within-gene partial-CDS check. It retains the existing shared-splice and validated annotation-support allowances, and keeps the seven known processed partial-CDS findings as failures. CSV rows identify it with `check_rule=within_gene_partial_cds` and `intentional_partial_cds_check=true`. The console separates `leap_rule_failures` from `partial_cds_failures`.

This document maps the original full-length selection predicates. Production now
also applies a protection-only partial-reference stage after that selection;
see [partial protection](partial_transcript_processing.md). It cannot create or
lengthen endpoints. The independent datacheck is unchanged and its per-isoform
partial-CDS rule can still flag grouped annotation unification. The
`partial_protection/` audit records the new production restrictions.

## Production sources

| Datacheck decision | LEAP reference |
| --- | --- |
| Full-length protein-coding eligibility requires both UTR types; other transcript types are non-coding | `bin/filter_out_partials.py`, `subworkflows/local/process_annotation/main.nf` |
| Protect eligible multi-exon CDS, single-exon CDS and non-coding exons | `bin/stop_cds_overlap.py` |
| Clipping includes the original grouped outer boundary and applies within and between genes | `overlap_region_with_feature` in `bin/stop_cds_overlap.py` |
| Fixed candidate field 4 of `.` bypasses clipping; transcript stats cannot change that decision | `bin/stop_cds_overlap.py`, `bin/build_chromosome_gff_and_stats.py` |
| Single-CDS genes in different-gene natural pairs have conditional flank protection | `split_gff`, `filter_exons_by_side`, and the singleton gene-field condition in `bin/stop_cds_overlap.py` |
| Test FE–LE and terminal–single-exon contacts, with positive-length BED overlap | `modules/local/intersect_first_and_last_exon_extensions/main.nf` |
| Safe terminal contacts share a candidate gene ID or match a zipped ordered FE/LE pair in the recorded list | `bin/overlap_safety_checker.py` |
| Recorded natural FE/LE pairs come from raw terminal intersections | `bin/find_natural_overlaps.py` |
| Gene IDs are grouped by terminal splice site in their recorded order | `bin/add_gene_and_transcript_id_to_USS.py` |

Clipping and terminal safety are independent. Annotation support bypasses clipping, but an unsafe terminal contact still fails. CDS from eligible full-length models are checked through LEAP clipping, including same-gene targets; they do not receive the additional per-isoform partial-CDS rule.

## Original context and processed output

`--leap-output-dir` must refer to the original run. Its `natural_overlaps/` and `fixed_extensions/` provide ordered groups, natural pairs, grouped original boundaries and candidate support. Missing or ambiguous required context is an error. The current GFF determines terminal geometry, so reverted internal extensions are not reinstated. Use corrected processed stats with processed GFFs.

For protection-only runs, use the run root containing final protected
`fixed_extensions/`, not `fixed_extensions_before_partial_protection/`.

A shorter isoform can gain bases already present in the grouped raw annotation. Those bases are outside LEAP's read-extension clipping guard, even when the final candidate is labelled read-supported. This is different from newly extending beyond the grouped raw boundary into protected sequence. An existing overlapping gene span alone never bypasses the terminal-safety predicate.

## Shared first exons across gene IDs

The raw ADORA3 (`ENST00000241356.5`) and TMIGD3 (`ENST00001035321.1`) records both explicitly have negative-strand `exon_number=1` starting at `111503005`, with different outer boundaries. They have different gene IDs but a shared terminal annotation region. The fixture `tests/data/leap_cross_gene/shared_first_exon_different_genes` preserves this example and its original candidate metadata. Its passing result does not grant an exemption to an unknown FE–LE pair.

## Regression validation

Run `python3 -m unittest discover -s tests`.

The tests execute LEAP's actual clipping function against the datacheck over 144 combinations of strand, extension end, same/different gene, candidate support, protected feature and boundary position. They also execute LEAP's terminal checker to validate ordered pairs, shared groups and single-exon field placement. Additional tests cover special flanks, missing/ambiguous context, annotation support failing terminal safety, within-gene non-coding guards, the shared ADORA3/TMIGD3 annotation and all seven retained partial-CDS failures.

Production LEAP scripts and internal-extension processing are unchanged by this datacheck refinement. This validates overlap policy; use the separate boundary and internal-exon datachecks for transcript shortening and splice-structure changes.
