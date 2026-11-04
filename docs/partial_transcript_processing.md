# Partial transcripts as LEAP processing context

LEAP classifies protein-coding models with both UTR types as full length and
models missing either UTR type as partial. Classification is unchanged. The
combined `_PC.gff` now carries both classes through the existing single-/multi-exon
partitioning, stranded annotation, terminal-exon extraction, splice-site grouping,
natural-overlap discovery and overlap decisions. Non-protein-coding annotation
continues through its existing separate protection channel.

Keeping whole partial models here makes their exons and CDS available to the
same decisions that use full-length models. It does not introduce a new blanket
exon-blocking rule, a different grouped boundary, or transcript-specific candidate
endpoints. Special single-CDS-gene classifications and natural-overlap pairs can
change because they now see the complete coding reference. The predicates that
use those data are unchanged.

The final writer receives the classified partial annotation using required
`--partial-annotation`. An empty file is valid; plain and gzipped GFF are
supported. The writer excludes each exact chromosome/strand/versioned transcript
key and its feature records before applying either extension end. Final GFF,
stats and gene bounds therefore contain only the same class of full-length
multi-exon source models as before. Partials can appear in intermediate candidate
and fixed-extension files, but are not released as extended models.

The independent datachecks do not repair output. Their existing extra
per-isoform partial-CDS restriction can still flag annotation-based unification
that the production grouped-boundary policy permits. Report those findings
separately; a zero count in that stricter check is not the acceptance criterion
for preserving LEAP's production rules.

## Validation

Run the regression suite with Python, pandas and bedtools available:

```bash
python3 -m unittest discover -s tests
```

`tests/test_partial_transcript_processing.py` checks classification, retained
partial features, actual CDS candidate selection, special flank inputs, both
strands, final feature/stat/gene exclusion, exact identifiers and input errors.
It also checks that a distinct partial CDS interval changes the existing
single-CDS-gene classification, rather than freezing eligibility to full-length
models only.
The bedtools-dependent cases check partial terminal and single-exon processing,
including a partial exon without CDS blocking an unknown terminal contact and a
partial model contributing a recorded natural-overlap pair.

Saved final candidate BEDs can support a bounded replay of the overlap and output
stages, but cannot recover every read candidate removed by the original run or
read evidence for newly represented partial-only splice sites. Such a replay
must be labelled accordingly. A complete validation requires rerunning LEAP
from the original reads with the expanded coding annotation.

## Validation recorded on 4 November 2026

All 68 automated tests passed, including 13 partial-processing cases. A Nextflow
smoke run executed the changed classification, partitioning, strand split and
final-writer wiring on both strands. With no partials, both GFF and stats were
byte-identical to the original writer in the two-strand fixture.

Bounded saved-endpoint replays retained the original source transcript sets,
CDS coordinates, exon counts and internal splice sites in all three supplied
datasets. No model became shorter than raw, no partial model entered final
output, and internal-extension processing dropped no source models.

| Dataset | Source models | Fixed candidates rechecked | Further production CDS/exon clipping needed |
| --- | ---: | ---: | ---: |
| Older human | 234,296 | 134,824 | 0 |
| Older mouse | 161,526 | 106,607 | 0 |
| Newer human (`leap_correct`) | 306,857 | 160,450 | 0 |

The unchanged independent datacheck reported these unique source transcripts
with its stricter partial-CDS findings:

| Dataset | Original output | Replayed output | Replayed output after internal-extension processing |
| --- | ---: | ---: | ---: |
| Older human | 670 | 367 | 180 |
| Older mouse | 143 | 82 | 51 |
| Newer human | 208 | 154 | 6 |

None of the remaining partial-CDS finding intervals intersects a read-supported
region beyond its grouped annotation boundary. They remain reported rather
than weakening the independent check.

Older human also has 66 additional `single_cds_gene_flank` checker rows at
CGGBP1 and SCO2. Each gene has one distinct CDS interval among full-length
models and two when a truncated partial-CDS interval is included. Production
therefore no longer applies its single-CDS-gene flank condition, while the
unchanged checker reconstructs that condition from full-length models only.
The candidates pass the actual production predicate with combined coding
context. This is an eligibility consequence of supplying additional reference
models, not a change to the predicate.

Replay used bedtools 2.31.1 and assumed overlap penalty 10 from the documented
example; original run parameter files were unavailable. The production bedtools
container remains unchanged. Partial-only splice groups used annotation-only
evidence. These results do not replace a complete rerun from the original reads.
