# Uploading LEAP annotations

`assets/update_LEAPed_annotation.pl` is a patched standalone copy of the supplied
Otter uploader. It retains the existing command interface and database connection
mechanism; it is not invoked by Nextflow or by the datachecks.

The script requires the same Ensembl, Vega and Otter Perl environment and configured
Otter dataset as the original uploader. Test a real dataset with a dry run first:

```bash
perl assets/update_LEAPed_annotation.pl \
  -in /path/to/extended_annotation.gff3 \
  -dataset YOUR_OTTER_DATASET -chr 12 \
  -out /path/to/review/mane_upload -dry_run
```

Plain and gzipped GFF3 input are supported as before. Database persistence is enabled
only with `-write`; the default is dry run, and `-dry_run` overrides `-write`.
The existing `_input_check.txt`, `_edits.log` and write-mode `_postedit_check.txt`
outputs remain. Review `PLAN:` entries in the edits log for target IDs, source IDs,
actions and final endpoints.

## Existing canonical extended transcripts

MANE status is read from database `remark` values `MANE_select` and
`MANE_plus_clinical`. Canonical extended models are recognized through the
`Ensembl_canonical_extended` database remark, also accepting its lowercase spelling.
These attributes do not have to occur in the LEAP GFF.

For a MANE extension, the uploader searches the same gene for an existing canonical
extended transcript with the same chromosome, strand, exon count, every protected
exon boundary, genomic CDS boundaries, and exon phases/end phases. Only the
lowest-coordinate exon start and highest-coordinate exon end may differ. This
preserves the full internal exon structure and splice-facing terminal boundaries,
and distinguishes different coding models even when they share an intron chain.
Coordinates are normalized to the genome so differing slice origins do not affect
matching. Single-exon models can match, with both outer ends permitted to differ.

| Situation | Action |
| --- | --- |
| Matching existing model already covers both requested ends | `already_covered`: keep it; create no transcript |
| Matching existing model needs further extension | `update_canonical_extended`: reuse its transcript/protein IDs and metadata |
| Existing model is farther at one end and this run is farther at the other | Keep the smaller start and larger end; update the existing model |
| No matching tagged model | `create_canonical_extended`: copy the MANE model, extend it, clear copy IDs and add the canonical extended remark |
| Multiple matching tagged models | Stop before editing that gene, reporting their IDs for review |

Different protected exon or CDS structures are not candidates for reuse. The
original MANE transcript retains its identity, structure and MANE remarks. Changed
terminal exons of canonical extended models are detached before editing so shared
exon objects do not change other models. Existing LEAP remarks are retained without
adding repeated instances of the same remark.

The whole gene is planned before edits. If both a canonical transcript and its MANE
source occur in the input, their selected extensions are combined into one update
using the farther endpoint at each end. Equivalent MANE sources requiring a new
canonical model are combined into one new copy. Processing order cannot undo a
longer selected endpoint. Repeating the same MANE upload reuses the existing result.
The script does not delete or consolidate pre-existing duplicate canonical models.

As in the supplied uploader, source models must pass the initial database/GFF
coordinate, terminal-exon and LEAP-tag checks to be scheduled for extension. Those
eligibility checks and the non-MANE extension behavior remain. The new planner also
rejects invalid coordinates and source shortening. Already completed gene writes
are not rolled back if a later gene has an error.

Post-upload checks compare the saved model with the planned endpoint union and
protected structure, verify canonical/MANE status, and verify that each original
MANE model stayed unchanged. They identify reused models by stable ID and newly
created copies by structure and final coordinates. Mismatch messages are diagnostic
log findings, as in the original uploader; they do not roll back database writes.

## Validation

```bash
python3 -m unittest discover -s tests -p test_update_leaped_annotation.py -v
```

The tests run the actual planner and editing routine with mock transcript, exon,
translation and Otter objects. They cover both strands, either/both ends, unchanged
or farther existing models, endpoint unions, repeated uploads, multiple MANE
sources, canonical input rows, CDS/exon incompatibilities, phase differences,
ambiguous matches, single-exon CDS, shared exon objects, slice offsets and
post-upload diagnostics. A separate compile-only test stubs the unavailable Otter
imports. Neither test connects to or writes a database. Real Otter persistence,
stable-ID assignment and server validation still require an environment-matched
dry run and subsequent database validation.
