# LEAP impact statistics and dashboard

Every LEAP run now writes `impact/index.html`, alongside machine-readable statistics. Open the HTML file directly in a browser, or share the **whole `impact/` directory** so its CSV downloads remain available. Charts and filters work offline, without installing a web server. The print button provides a compact view for saving as PDF.

The report runs independently of manual-review preparation. With datachecks enabled, it includes the completed datacheck summary and finishes before the final gate. A failing datacheck does not suppress the impact report. When checks are explicitly disabled, statistics still run and the dashboard states `NOT_RUN` for validation. An incomplete statistics report is published with diagnostics and fails the run after review preparation completes.

## What the report measures

Filter by chromosome, strand and result view:

- **Original LEAP output**: every retained output transcript, including zero-extension models.
- **External output**: the result of applying the existing internal-extension processor in memory. Internal ends are restored; models without CDS are excluded. Zero-extension retained models remain in the denominator.
- **Internal additions only**: only positive additions removed by that processor. A transcript with one internal and one external extended end appears in both result views, with each view counting its corresponding added sequence. Internal classification without positive growth does not count as an internal extension.

The report reuses the existing processor's classification, first-seen CDS tie behaviour and support handling. It does not change extension selection, modify annotations, run the repair asset or publish a new processed GFF. Manual-review preparation continues to supply the actual processed annotation files. Datacheck outcomes apply to the original pipeline output; the calculated external view is not a separately validated release.

Headline measures include extended transcripts, affected genes, 5′/3′ extension counts, added bases, support, blacklist counts, length distributions and positive-extension quantiles. The dashboard highlights:

- **Newly covered genomic bases**: all extension intervals merged across transcripts and both strands, minus **every raw GFF feature interval on either strand**. This includes gene and transcript spans, exons, CDS, UTRs, codons and unfamiliar feature types. A genomic position counts once per result view, regardless of how many transcripts cover it or which evidence supports their ends.
- **New exonic bases**: the same extension union minus raw **exons** on either strand. This includes new exon coverage within a previously annotated gene or transcript span, which is excluded from the stricter “no previous feature” measure.
- **Unique extension bases**: all extension intervals merged by chromosome across transcripts and both strands, including positions already annotated in another model.
- **New read-supported genomic bases**: the read-supported subset outside every raw feature on either strand, merged across transcripts and strands. A position counts once even if several read-supported extensions cover it.
- **Added bases across transcripts**, retained in the exports: the sum of each transcript's additions. A shared 100-base addition in ten isoforms contributes 1,000 bases here and 100 bases to the union measure. It is not used as the headline newly covered base count.

All intervals are inclusive. Added intervals exclude the original boundary base: a start moving from 100 to 90 adds bases 90–99. Negative-strand 5′ ends correspond to the high-coordinate boundary. **New genomic coverage measures merge across strands**, so the same chromosome position counts once even if positive and negative transcripts both extend over it. Chromosomes remain distinct. A strand filter selects which extensions contribute; the strict raw-feature exclusion still considers both strands. An additional same-strand raw-feature comparison is exported as `new_same_strand_feature_bases`.

The older `unique_added_bases`, `unique_read_supported_bases` and `novel_read_supported_bases` export columns retain their per-strand definitions for compatibility; opposite strands count separately there. `unique_genomic_extension_bases`, `newly_covered_bases` and `newly_exonic_bases` provide the across-strand genomic counts. Gene groups use chromosome, strand and exact `gene_id`; transcript counts preserve exact versioned IDs and chromosome/strand copies, including PAR copies. Chromosome filters describe the supplied run scope, not a genome-wide denominator if only a subset was processed.

Each selected end has one evidence source: `reads` or `current_annotation`. A transcript may have read support at one end and annotation support at the other; these ends are counted separately. `none` and `blacklisted` describe unextended ends. Repeated identical labels on an end do not increase its base count; conflicting read/annotation labels on one end are an input error. Read counts are retained verbatim in the transcript CSV; partition/gene tables include the minimum and maximum recorded read-support counts. These values are **not summed into a number of unique reads or molecules**.

Length mean, median, p90, p95 and maximum describe **positive extensions only**, separately for each end and stage. Percentiles use linear interpolation. Run-wide statistics are computed from individual end lengths rather than averaging chromosome medians. Interactive histograms and headline counts follow the filters; `extension_lengths.csv` contains the exact run-wide quantiles.

This report measures annotation growth. It does not measure read coverage, gene expression, new transcript discovery, new coding sequence, functional benefit or assembly identity. A positive impact number is not evidence that validation passed. Missing CDS models are counted as processor exclusions, separately from internal additions. External and internal **unique** base unions may overlap between isoforms, so their union totals must not be added together to reconstruct the original unique total.

## Files

| File | Contents |
| --- | --- |
| `index.html` | Offline interactive dashboard, original/external/internal comparison, largest examples, datacheck outcomes and counting rules |
| `summary.json` | Complete aggregate metrics, scope, definitions, provenance, diagnostics, plot data and supplied datacheck results |
| `summary.txt` | Short written totals for sharing or preparing slides |
| `overall_statistics.csv` | Run-wide totals for each of the three views |
| `genomic_coverage.csv` | Original/external/internal unique genomic extension, newly covered and new exonic bases by chromosome and positive/negative/both views; `both` is recomputed as a union rather than a sum |
| `newly_covered_intervals.csv` | Disjoint inclusive intervals outside every raw feature, once per stage and chromosome; the base column sums to the strict run-wide newly covered total |
| `transcripts.csv` | One row per original output transcript: exact ID, gene, strand, raw/original/external bounds, `chromosome:start-end`, support, read counts, internal/external lengths and exclusion flag |
| `genes.csv` | Gene impact by stage, chromosome and strand |
| `partition_statistics.csv` | Metrics by stage, chromosome and strand, including raw/eligible denominators |
| `support.csv` | Transcript-end counts by support category, end and stage |
| `extension_lengths.csv` | Run-wide length summaries by end and stage |
| `diagnostics.csv` | Missing/ambiguous inputs, inconsistent GFF/stats or processing errors |

HTML embeds aggregate data and up to 20 largest examples per partition/view. The example table shows the largest 50 matching those stored examples; it is not a search across every transcript. Full transcript and gene tables are in the CSVs. Free-text CSV cells beginning with spreadsheet formula characters receive a leading apostrophe for safe spreadsheet viewing.

## Rebuild from an existing run

```bash
python3 assets/build_leap_impact_report.py \
  --leap-output-dir /path/to/leap_run \
  --chromosomes 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y \
  --output /path/to/new_impact_report \
  --label 'Human LEAP impact'
```

`--strand positive|negative|both` defaults to both. Use the exact chromosome names in that run. One raw annotation and one GFF/stats pair per selected chromosome/strand are required. Plain and gzipped pipeline files are supported; ambiguous plain/gzip copies are rejected. Empty final GFFs and header-only stats are accepted as zero models, while completeness remains the responsibility of the completeness datacheck.

The command automatically uses `datachecks/summary/summary.json` if present. Otherwise validation is labelled `NOT_RUN`. Supply `--datacheck-summary /path/to/summary.json` to use a matching summary saved elsewhere; only supply results from the same original run and scope. A failed datacheck remains visible but does not prevent rebuilding descriptive statistics.

The output directory must be new; existing directories and linked destinations are rejected. Source files remain untouched. Exit `0` means all selected statistics partitions were analysed, not that biological validation passed. Exit `2` means missing, ambiguous, unreadable or inconsistent required inputs. A readable partial report is retained with `INPUT_ERROR` and explicit completed/expected counts. Statistics from failed partitions are omitted, never treated as zero-success results.

## Pipeline wiring and tests

`LEAP_IMPACT_REPORT` consumes declared scope, raw GFF and final GFF/stats channels. It does not discover asynchronously published files. The original production and datacheck modules remain responsible for biological decisions. The reporting process stages only the relevant annotation/statistics helpers and the HTML template; database connection files are not inputs.

The datacheck gate also requires a complete impact summary. If biological checks are disabled, `LEAP_IMPACT_GATE` performs that report-completeness gate after publication and review preparation. This keeps report errors visible in both execution modes.

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s tests -p test_build_leap_impact_report.py -v
LEAP_NEXTFLOW=/path/to/nextflow python3 -m unittest discover -s tests -p test_impact_report_nextflow.py -v
```

Tests cover union versus transcript totals, merging across opposite strands, distinct chromosomes, all-feature and exon-only subtraction, inclusive boundaries, gene spans, unstranded/unknown features, same-strand comparisons, distinct support at opposite ends, conflicting support on one end, internal/external classification, no-CDS exclusions, zero growth, blacklisting, gzip, empty partitions, invalid/missing/ambiguous inputs, output protection, safe HTML and real workflow publication before gating. The existing datacheck/manual-review smoke tests also assert dashboard publication when checks fail.
