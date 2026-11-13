# Duplicate transcript example

`cox14/` contains two real COX14 transcripts from the published raw and final
chromosome 12 annotations in `/Users/rmerritt/temp/leap_new_real_run/`.
The source datasets were read without modification; these are subsets containing
the two models and their gene record. Other transcripts at this gene are omitted.

| Transcript | Raw bounds | Final bounds |
| --- | --- | --- |
| ENST00000550487.6 | 12:50112236-50120453 | 12:50112029-50120463 |
| ENST00001138219.1 | 12:50112029-50120463 | 12:50112029-50120463 |

Both are on the positive strand. Their final exon intervals are
`50112029-50112301;50120036-50120463`. Both have the same CDS interval
`50120044-50120217`, phase `0`. Their raw CDS structures already matched, but
their raw exon boundaries differed. Extension brings the first transcript to the
second transcript's annotation boundaries, producing a duplicate group by the
exon-and-CDS definition. This does not establish that the extension was unsafe.

```bash
python3 assets/check_duplicate_transcripts.py tests/data/duplicate_transcripts/cox14/raw.gff3 raw_duplicates.csv
# Exit 0: no duplicate group.
python3 assets/check_duplicate_transcripts.py tests/data/duplicate_transcripts/cox14/extended.gff3 extended_duplicates.csv
# Exit 1: one group containing two transcripts.
```

Use the individual files rather than the fixture directory: raw and final are
two stages of the same models, not independent output partitions.
