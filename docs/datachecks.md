# Automatic LEAP datachecks

Datachecks run by default after the normal extension outputs have been generated. Configure `datacheck_reference_lengths` with a FASTA `.fai` or two-column sequence-length file for the same reference assembly. Plain and gzipped files are accepted; chromosome names must match exactly. LEAP rejects missing/unreadable reference paths before starting production tasks. `run_datachecks = false` explicitly opts out.

## What runs

| Task alias | Check |
| --- | --- |
| `CHECK_BOUNDARIES` | No original transcript boundary moves inward |
| `CHECK_INTERNAL_EXONS` | Exon count and all splice-facing boundaries are preserved |
| `CHECK_PRESERVATION` | Complete transcript feature and metadata preservation |
| `CHECK_STATS_GFF` | Final coordinates, lengths, support, read counts and tags agree, including fixed-candidate provenance |
| `CHECK_STRUCTURE` | GFF parent graph and feature containment |
| `CHECK_COMPLETENESS` | Every eligible raw transcript appears in GFF and stats |
| `CHECK_PARTIAL_PROTECTION` | Only permitted candidate changes occur, with matching blocker evidence |
| `CHECK_REFERENCE_COORDINATES` | Raw and final GFF features fit exact reference sequence lengths |
| `CHECK_CDS_OVERLAPS` | LEAP overlap rules plus the existing strict within-gene partial-CDS diagnostic |
| `CHECK_DUPLICATE_TRANSCRIPTS` | Non-failing GFF retention/removal recommendations for identical exon/CDS models |

The suite does not change the standalone datachecks' comparison rules. Partial-CDS contacts within a verified existing grouped annotation boundary are reported as allowed; unexempted new read-driven contacts beyond it fail. Automatic checking does not run the standalone internal-extension processor or imply that its release-stage output has been checked.

## Channel wiring and missing partitions

`PROCESS_ANNOTATION` additionally emits its published raw annotation. `EXTEND_TRANSCRIPTS` additionally emits a tuple containing partition ID, GFF, stats, protected FE/LE candidates and protection audit. `OVERLAP_READS` exposes that tuple and the original pre-protection FE/LE candidates. Extension processing is unchanged; annotation export now uses the bundled Python exporter and a declared database-secrets JSON input.

`DATACHECKS` receives those channels, natural pairs and the expected chromosome/strand scope derived independently from the reads chromosome sample sheet. Both strands are expected for every selected chromosome. Channels are collected with `toList`, which retains empty streams, and matched against that expected scope. Missing outputs are represented by empty input lists, so completeness still runs for an absent partition. Unexpected and duplicate scope records are retained as summary input errors.

Each task stages its own declared inputs into a private pipeline-shaped `context/` directory. It does not read `params.output_directory`, wait for publish operations, modify production files or infer input scope from whichever outputs happen to exist. Staged check scripts and the adapter are explicit task inputs so changes invalidate their task caches. The same reusable process is included under ten named aliases, with at most four concurrent tasks per alias. Resource labels use existing pipeline configuration; execution requires Python 3 and no additional Python packages.

## Results and failure handling

Every check task retains a CSV, stdout, stderr and a uniquely named status JSON. Ordinary checker exits `0`, `1` and `2` are recorded as `PASS`, `FAIL` and `INPUT_ERROR`; duplicate-transcript exit `1` is instead recorded as non-failing `REVIEW` with gate exit `0`; the adapter task itself completes so other checks can finish. Unexpected checker errors, abnormal exits and missing required CSV reports are recorded as input/execution errors. An error placeholder CSV points to diagnostics when the standalone script could not produce a biological report.

When manual-review preparation is enabled, the final gate also waits for its completed bundle. A datacheck failure therefore does not cancel review preparation before its tracks and launch files are available. Review processing consumes copies and does not change the standard LEAP annotation evaluated by these checks.

The summary requires exactly one result per expected partition/check. Missing, duplicated, unrecognized or malformed result records cannot silently pass. Reports are published under:

```text
datachecks/
  1_positive_boundary/
    report.csv
    stdout.txt
    stderr.txt
    1_positive_boundary.status.json
  ...
  summary/
    summary.html
    summary.csv
    summary.json
    verdict.json
```

The separate final gate runs only after all check results and the combined summary exist. It returns `0` when everything passes or is explicitly inapplicable, `1` for detected violations, and `2` when any input/execution error prevents complete assessment. Nextflow reports an unsuccessful workflow for either nonzero gate result. Reports remain available for review and `-resume` can reuse completed validation tasks. Upstream production/task-staging failures still follow Nextflow's normal error handling; this mechanism delays failure caused by datacheck findings, rather than overriding production failures.

Some standalone scripts require at least one final transcript. For an entirely empty final GFF, the adapter first runs scoped completeness. Only a passing preflight with zero eligible transcripts, zero final models and zero stats rows permits `NOT_APPLICABLE` for boundary, internal-exon, preservation, stats/GFF, structural and overlap comparisons. The preflight CSV and logs are retained. Completeness, partial protection, reference coordinates and duplicate recommendations always run. Nonempty gene-only files, missing eligible models or malformed inputs are not treated as valid empty output.

## Validation

Run all Python tests:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s tests
```

The real Nextflow smoke tests use tiny published-layout fixtures and the actual `DATACHECKS` subworkflow; no database, reads processing or containers are needed. They run when `nextflow` is on PATH or an executable is supplied:

```bash
LEAP_NEXTFLOW=/path/to/nextflow python3 -m unittest discover -s tests -p test_datachecks_nextflow.py -v
```

They verify all ten checks on both strands, retained reports after both biological failure and input error, missing-partition completeness and valid zero-eligible output. The standalone adapter tests also cover checker crashes, missing reports, summary completeness and verdict ordering. These tests exercise real task execution and published outputs rather than only inspecting module text.

## Duplicate recommendations

The tenth module runs `check_duplicate_transcripts.py` with `--recommend-removals`
on each final chromosome/strand GFF. Its report recommends GFF-only retention or
removal from MANE, canonical extended and blacklist tags, followed by numeric
stable ID. Multiple MANE models stay protected. Unknown/tied numeric IDs require
manual review. Recommendations apply across genes and do not modify annotation.

`REVIEW` is accepted only for this checker with a duplicate-finding exit of `1`.
A duplicate-only run completes successfully with overall `REVIEW`; missing or
malformed input, a missing CSV, and biological failures in other checks still fail
the gate after all reports publish. Empty final GFFs pass this audit with an empty
CSV. Duplicate reports do not establish redundancy in the saved database, where
the uploader preserves original MANE coordinates.

## Impact report

The automatic workflow also writes `impact/index.html` and CSV/JSON impact statistics from declared annotation/stats channels. Reporting waits for the completed datacheck summary; the final gate waits for both impact reporting and manual-review preparation. Biological failures retain the dashboard. Missing or inconsistent reporting inputs produce an incomplete dashboard with diagnostics and fail the final gate. See [impact statistics](impact_dashboard.md) for definitions and the standalone rebuild command.
