# Annotation exporter and connection JSON

LEAP bundles the modern Python exporter from the db2gxf project at commit
`dff0ec308a9fd738e29470923a43036ef7bf942e`. Source files are under `bin/gxf_dump/`,
with `bin/gxf_dump.py` as the command-line entry point. The upstream MIT notice is
in `licenses/db2gxf-LICENSE`. No external checkout is required on the cluster.

LEAP-specific changes expose `MANE_select` and `MANE_plus_clinical` database remarks
as `MANE_Select` and `MANE_Plus_Clinical`, alongside dedicated MANE attributes.
Canonical extended and the three blacklist-related remarks are normalized to their
standard tags. Recognition uses complete known values, not free-text substring
matches. Hidden notes cannot independently grant MANE status. The Y-PAR
compatibility helper projects complete default-assembly counterpart genes
onto Y using the database PAR mappings, preserves exon/CDS phases and traversal,
and supplies `_PAR_Y` GFF feature IDs while retaining stable `gene_id` and
`transcript_id` metadata. Stored copies are not duplicated. Genes crossing a
projection boundary are excluded from Y projection with a diagnostic, matching
the complete-gene transfer seen in the supplied legacy X/Y exports. Other export
selection, sequence reconstruction, feature coordinates and ordering follow the
pinned exporter. Live CLI/task diagnostics suppress driver messages that might
contain connection values.

## Migration

Copy `examples/db_secrets.example.json` to a private location, fill in the same
connection used by LEAP previously, and launch with:

```bash
nextflow run main.nf -profile hpc_slurm \
  --db_secrets /private/path/db_secrets.json \
  --datacheck_reference_lengths /path/to/reference.fa.fai \
  --reads /path/to/reads.bam --reads_index /path/to/reads.bam.csi \
  --cage /path/to/cage.bed --polya /path/to/polya.bed \
  --output_directory /path/to/leap_results
```

The old `--host`, `--port`, `--user`, `--pass` and `--dbname` parameters are rejected.
The JSON needs `host`, `port` (1–65535), `user`, `database` and a string `password`
(which can be empty). `dbname` and `pass` are accepted aliases, and `dialect`
defaults to `mysql`. Conflicting aliases are errors. File/schema validation occurs
before production tasks are launched. JSON values are not interpolated into task
commands. Only the selected source file is staged as a declared task input;
connection files are not published with annotation or reports. Keep populated
connection files out of Git; the supplied example contains placeholders only.

The annotation process declares a Conda environment with Python 3.11, SQLAlchemy
2 and PyMySQL 1. When Conda is disabled, provide Python 3.10+ and equivalent
`SQLAlchemy>=2,<3` and `PyMySQL>=1.1,<2` packages. Other LEAP helper requirements
remain. Annotation extraction no longer requires the Perl exporter or MySQL CLI;
the separate database uploader still requires Otter/Perl.

## Preserved annotation processing

Each task exports its explicitly selected chromosome, both strands, in GFF3 with
versioned identifiers and `--vega`. It retains exporter traversal and exon-rank
order. Do not coordinate-sort this file before LEAP's terminal-exon extraction.
Canonical extended models remain available as annotation and extension models.

The Python wrapper then performs the existing whole-gene readthrough removal and
whole-model OTT filtering, retaining OTT exon identifiers attached to retained
models. The blacklist still selects current protein-coding transcripts with
upstream-ATG, CAGE-supported-TSS or RP-supported-TIS remarks. Attribute type is
resolved by `remark` code instead of a hard-coded numeric ID. Queries use the same
JSON connection and selected default-assembly sequence-region IDs. Outputs retain
`current_annotation_<chromosome>.gff3` and the two strand-specific blacklist names.
The published raw annotation remains feature-only, as in the previous wrapper.

## Validation

The supplied chromosome 21 comparison contains 1,104 genes, 9,495 transcripts and
158,962 feature records in each export. Coordinates, phases, membership, traversal
order and non-tag attributes match. The modern file exposes 135 canonical extended
models; the other differences are tag ordering and repeated tag normalization.
Neither supplied file exposes MANE status, so these historical files alone cannot
validate the newly added MANE-remark mapping. A fresh chromosome 21 export with
the bundled exporter was therefore obtained on 13 November 2026. All 158,962
records match the supplied Perl export in traversal order, coordinates, phases,
membership and non-tag attributes. The fresh export contains 214 MANE Select
transcripts and one MANE Plus Clinical transcript; their tags agree with the
database's dedicated attributes and exact MANE remarks, with zero discrepancies.

To repeat this validation, export chromosome 21 with the bundled version, using
the same database snapshot where possible:

```bash
python3 bin/gxf_dump.py --db-secrets /private/path/db_secrets.json \
  --chromosome 21 --vega --format gff3 --output modern_chr21_with_mane.gff3
```

Live readthrough/blacklist queries and whole-model filtering were also checked.
Both chromosome 21 exports produce the same 158,059 retained records, in the same
order apart from tags; the queries identify eight readthrough genes and 36
blacklisted transcripts. The supplied unfiltered Perl export and the fresh export
both contain 11 gene rows without transcripts. Ten are removed by OTT filtering;
one existing empty gene, `ENSG00000214955.5`, remains in both filtered copies. This
input finding is recorded separately and is not repaired by the exporter upgrade.

LEAP preprocessing of the supplied exports produces identical full-length,
partial, coding and noncoding feature collections apart from tags. All four
terminal unique-splice-site tables are byte-identical. Local tests additionally
exercise SQL extraction against synthetic SQLite schemas, tag interpretation,
filtering, connection validation and credential-safe errors.

The additional supplied X/Y validation compares 10,230 projected PAR feature
records exactly, including phases, IDs, parent links and metadata. The legacy Y
file has 11 additional records from two genes absent from the supplied filtered X
file; they cannot be reconstructed from that filtered X input. The production
exporter queries the unfiltered counterpart annotation before applying LEAP's
filters. A fresh production Y export includes 64 PAR genes and 640 PAR transcripts;
its structural and reference-coordinate checks pass. The older Y run is a
different annotation baseline: 57 transcript IDs differ in each direction, while
the 583 shared PAR transcript IDs retain identical locations. It is not evidence
of exact parity with that older snapshot. An independent current-database
manifest agrees with all 640 fresh PAR transcript IDs, gene memberships,
coordinates and strands, with no missing or mismatched models. The fresh
chromosome 21 comparison is the matched export-parity validation.

The regression suite ran 291 tests: 288 passed and three existing bedtools tests
were skipped because bedtools was unavailable. Nine real Nextflow integration
tests cover annotation staging/preflight and all ten automatic datachecks,
including publication before gating, non-failing duplicate reviews and failing
input errors. A complete read-processing/extension run was not performed as part
of this export validation.
