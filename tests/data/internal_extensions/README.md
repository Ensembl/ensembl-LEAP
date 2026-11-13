# Internal-extension processing fixtures

Each folder is a small LEAP output directory with raw annotation, extended annotation and extension stats. All cases include:

- `anchor.1`, defining the gene's extreme CDS segments;
- `source.1`, with 50 added bases at each end and an unchanged middle exon;
- `utr_only.1`, without CDS, which is intentionally excluded from processed output.

The source's 5′ extension is supported by `current_annotation` with read count 7; its 3′ extension is supported by `reads` with read count 3. This makes wrong-end resets on the negative strand visible.

| Folder prefix | Source behavior |
| --- | --- |
| `external` | Keep both extensions and both original support categories/counts |
| `internal_left` | Restore the genomic left terminal exon from raw |
| `internal_right` | Restore the genomic right terminal exon from raw |
| `internal_both` | Restore both ends and keep a zero-extension stats row |

Every case has positive- and negative-strand versions. On the positive strand, left means 5′; on the negative strand, left means 3′. `expected.json` records the exact processed start/end coordinates for each case.

Run one case from the repository root:

```bash
python3 assets/process_internal_extensions.py tests/data/internal_extensions/internal_left_negative /tmp/leap_internal_left_negative_review
```

Use a new output directory for each run. Compare `external_extensions/gff/1_negative_extended.gff` and `external_extensions/stats/1_negative_stats.tsv` against the input. The 3′ extension should be zero with support `none` and read count `0`, while the 5′ extension remains 50 with its original annotation support and read count 7. The original extended source records are retained under `external_extensions/potential_internal/`.

Run the regression tests with:

```bash
python3 -m unittest discover -s tests -p 'test_process_internal_extensions.py'
```
