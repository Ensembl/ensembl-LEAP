# Partials protect extensions without supplying endpoints

Partial protein-coding models are reference obstacles only. The normal LEAP
pipeline uses its original full-length annotation for exon partitioning,
terminal extraction, shared-splice-site grouping, annotation endpoints, gene-ID
groups, natural-overlap discovery and overlap selection. Partial models cannot
supply a longer endpoint or change those existing allowances and classifications.
The combined `_PC.gff` remains a reference/compatibility output of classification;
it is not consumed by endpoint generation.

## Protection contract

`bin/protect_partial_extensions.py` runs **after** the original unsafe-overlap
fixer, before the final GFF/stats writer. Running it at that point prevents a
blocked extension from enabling a longer extension at another end.

- Annotation-only candidates and single-exon context rows are unchanged.
- Only read-supported sequence beyond the original grouped annotation boundary
  is checked against additional partial protection. Existing contacts inside the
  original annotation are not newly rejected by these added checks.
- All same-chromosome, same-strand partial CDS are protected, including
  single-exon CDS. The CDS guard includes the original annotation boundary,
  matching the existing clipping predicate's boundary-touch behavior.
- The existing single-CDS-gene eligibility is computed from full-length models
  and original natural pairs. Applicable partial exon flanks are added to its
  protections without changing that eligibility. Partials cannot disable an
  existing restriction by adding a distinct CDS interval.
- Opposite terminal-exon and single-exon contacts with partials use the original
  shared-gene and ordered natural-pair allowances. Terminal contacts require
  positive-length BED overlap, as in the existing terminal checks. Partials
  do not add natural pairs or enlarge candidate gene groups. This does not
  introduce a blanket rule protecting every multi-exon UTR exon.
- A candidate can stay unchanged, choose an existing shorter read endpoint,
  or return to its original annotation fallback. Longer alternatives left in
  metadata by the original terminal fixer are excluded from consideration.
  No endpoint can extend farther than the original selected candidate.
- Unchanged candidate records are copied exactly. Changed candidates preserve
  their splice site, gene group and splice-support read count. Annotation
  fallback is marked as annotation-supported. The final writer computes lengths
  and support normally; a reverted end with no remaining extension gets
  `none` support and zero reads. The opposite end is unchanged.

Partials remain excluded from released GFF/stats. The final writer retains its
required `--partial-annotation` argument as a defensive exclusion of exact
chromosome/strand/versioned transcript keys, even though its normal source input
now contains only the original full-length models.

## Command and outputs

```bash
python3 bin/protect_partial_extensions.py \
  selected_FE.bed selected_LE.bed \
  current_annotation_FT_strand.gff partial_annotation_PT.gff \
  original_natural_overlaps.txt 12_positive
```

The command writes `12_positive_FE_protected.bed`,
`12_positive_LE_protected.bed` and `12_positive_partial_protection.tsv`.
The audit identifies changed candidate records, selected and protected endpoints,
annotation fallback, and blocking partial transcript/feature coordinates.
Plain and gzipped inputs are accepted, including empty partial files. Invalid
records and existing output paths are rejected. Exit codes are 0 on success
(including successfully shortened candidates) and 2 for input/execution errors.
This is production protection logic; the independent datachecks are unchanged.

Nextflow publishes original selected candidates in
`fixed_extensions_before_partial_protection/`, final protected candidates in
`fixed_extensions/`, and audit TSVs in `partial_protection/`. Final GFF and stats
use the protected candidates. Use a fresh output directory when testing this
revision so older candidate files do not mix with the new final context.

## Validation

Run `python3 -m unittest discover -s tests` with pandas and bedtools available.
Tests cover original endpoint grouping, annotation-only preservation, both
strands and transcript ends, shorter reads/fallback, stale longer alternatives,
partial CDS, unique flank exons, partial single-exon models without CDS,
original natural-pair/gene-group ordering, strand/chromosome separation,
boundary contacts, gzip, input errors, feature exclusion and final stats.

On 4 November 2026, all **85 tests** passed. A real two-strand Nextflow smoke
run exercised classification, full-length partitioning, strand splitting,
partial protection and final writing, and compiled the full subworkflow imports.

An exact comparison of the post-selection stages used original `leap_correct`
chromosome 12 candidates. Rebuilt baseline GFF/stats reproduced the original
output, and regenerated natural-pair sets matched the original on both strands.
All 15,918 source transcripts were retained; 15,914 have identical normalized
features and stats. No endpoint extends farther than the original output, and
no model is shortened relative to raw. CDS, exon counts and splice-facing
boundaries are unchanged throughout.

Only four source transcripts have shorter ends:

| Gene | Transcript | Original endpoint | Protected endpoint | Reason |
| --- | --- | ---: | ---: | --- |
| PTHLH | ENST00000539239.5 | 27970507 | 27969983 | Partial CDS blocks the longer read candidate |
| USP15 | ENST00001130171.1 | 62327962 | 62326231 | Partial CDS requires original annotation fallback |
| B4GALNT3 | ENST00001155883.1 | 555556 | 555431 | Partial CDS touches the original annotation boundary |
| B4GALNT3 | ENST00001155887.1 | 555556 | 555431 | Same protected splice group |

TSPAN9 ENST00000407263.2 stays at its original start, 3201177, with no added
5′ extension. The previous combined-context implementation's 383 longer
transcripts are no longer generated. Existing annotation unification at MAPKAPK5
is preserved. The original per-isoform datacheck reported one SCNN1A contact
inside the grouped annotation boundary, rather than in new read-driven sequence.
The 13 November datacheck refinement reports that contact as allowed grouped
annotation, while retaining failures beyond the original boundary. Production
partial protection is unchanged.

This validation uses saved **final** candidates to test the new post-selection
stage exactly. It does not rerun the original BAM/database stages. The separate
Nextflow smoke checks the restored full-length channel wiring; a full chromosome
12 run from the same reads remains the end-to-end confirmation. Local validation
used bedtools 2.31.1; the production container was not changed.
