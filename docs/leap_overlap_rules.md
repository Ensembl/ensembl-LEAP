# LEAP overlap-rule contract

`assets/check_cds_overlaps.py` validates the overlap decisions in the current LEAP code. It does not add a blanket prohibition on entering another gene's exon or use overlapping gene spans as a natural-overlap whitelist.

The sole deliberate extra overlap restriction is the within-gene partial-CDS check. It retains the existing shared-splice and validated annotation-support allowances and reports sequence within a verified original grouped annotation boundary as allowed. Unexempted read-driven contacts beyond that boundary remain failures, including the seven saved unsafe-read regression examples. CSV rows identify it with `check_rule=within_gene_partial_cds` and `intentional_partial_cds_check=true`. The console separates `leap_rule_failures` from `partial_cds_failures`.

This document maps the original full-length selection predicates. Production now
also applies a protection-only partial-reference stage after that selection;
see [partial protection](partial_transcript_processing.md). It cannot create or
lengthen endpoints. The independent datacheck reports pre-existing grouped annotation unification
as an allowance and continues to fail new read-driven partial-CDS overlaps. The
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

Candidate selection follows the final writer: for each chromosome, strand, FE/LE kind and internal splice site, the **last row in the BED file wins**. Selection occurs before gene-membership or geometric validation. Earlier rows cannot supply a more permissive support label, an alternative gene group, or a containing interval when the selected last record is incompatible. Support, original boundary and ordered gene IDs are taken together from that single record, matching the stats/GFF consistency check. Missing records, duplicate context files and incompatible selected records remain input errors. Diagnostics include the splice site and, for a selected record, its source file and line.

Changed ends must fit within the selected candidate; exact endpoint equality is not required because internal-extension processing can restore part of an extension. Unchanged terminal exons may extend beyond a shorter candidate because the writer never shrinks their original annotation.

The DDX41 case (`ENSG00000183258.18`, chromosome 5, negative-strand LE splice site 177511927) has two rows spanning `177511569-177511927`. One uses annotation support (`.`), followed by a read-labelled row whose original boundary is also 177511569. Both coordinates are valid. The writer deterministically retains the final read-labelled row; different support labels alone do not make the context ambiguous. The minimal fixture is `tests/data/leap_cross_gene/ddx41_duplicate_candidates`.

For protection-only runs, use the run root containing final protected
`fixed_extensions/`, not `fixed_extensions_before_partial_protection/`.

A shorter isoform can gain bases already present in the grouped raw annotation. Those bases are outside LEAP's read-extension clipping guard, even when the final candidate is labelled read-supported. This is different from newly extending beyond the grouped raw boundary into protected sequence. An existing overlapping gene span alone never bypasses the terminal-safety predicate.

## Partial-CDS grouped annotation allowance

The within-gene partial-CDS check intersects each newly added isoform UTR with the original grouped annotation extent before applying the existing annotation-support and shared-splice exceptions. Bases inside the grouped boundary are allowed, including the boundary base itself. Bases beyond it remain eligible to fail. A mixed contact reports both portions without allowing the annotation-covered part to excuse the remaining read-driven overlap. No cross-gene clipping or terminal-safety predicate changes.

Eligible raw full-length, multi-exon terminal exons must verify the annotation extent at the same chromosome, strand, FE/LE kind and internal splice site, within the candidate's recorded gene group. Only these original endpoint-generating models may contribute; partial models, noncoding models, single-exon models and matching internal exons cannot supply a new allowance. The raw extreme is capped by the fixed BED's recorded original boundary. When annotation fallback has replaced the alternatives with `.`, eligible raw terminal donors recover the existing extent. Missing or ambiguous run context still errors; unsupported candidate metadata alone cannot grant an allowance.

CSV rows retain `check_rule=within_gene_partial_cds` and add `grouped_annotation_allowed_intervals`, `grouped_annotation_boundaries` and `grouped_annotation_donors`. Wholly annotation-covered contacts have status `allowed_grouped_annotation`; mixed contacts retain `new_cds_overlap` and the exact residual `failing_intervals`. The JBrowse builder routes grouped allowances to a dedicated labelled track.

For SCNN1A `ENST00000396966.6`, the isoform's raw end is 6374803, the original grouped endpoint is 6375224 and the selected end is 6375258. Its contact at `12:6374804-6374837` is therefore allowed as existing grouped annotation, supported by raw `ENST00000360168.8`. Read sequence beyond 6375224 receives no grouped-annotation exemption.

## Shared first exons across gene IDs

The raw ADORA3 (`ENST00000241356.5`) and TMIGD3 (`ENST00001035321.1`) records both explicitly have negative-strand `exon_number=1` starting at `111503005`, with different outer boundaries. They have different gene IDs but a shared terminal annotation region. The fixture `tests/data/leap_cross_gene/shared_first_exon_different_genes` preserves this example and its original candidate metadata. Its passing result does not grant an exemption to an unknown FE–LE pair.

## Regression validation

Run `python3 -m unittest discover -s tests`.

The tests execute LEAP's actual clipping function against the datacheck over 144 combinations of strand, extension end, same/different gene, candidate support, protected feature and boundary position. They also execute LEAP's terminal checker to validate ordered pairs, shared groups and single-exon field placement. Additional tests cover last-row precedence, incompatible final candidates, special flanks, missing/ambiguous files, annotation support failing terminal safety, within-gene non-coding guards, the shared ADORA3/TMIGD3 annotation and all seven retained partial-CDS failures.

Production LEAP scripts and internal-extension processing are unchanged by this datacheck refinement. This validates overlap policy; use the separate boundary and internal-exon datachecks for transcript shortening and splice-structure changes.
