import csv
import gzip
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from datacheck_fixtures import write_run
from test_check_transcript_preservation import annotation, extend

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / 'assets/repair_transcript_features.py'
sys.path.insert(0, str(ROOT / 'assets'))
spec = importlib.util.spec_from_file_location('repair_asset', SCRIPT)
repair = importlib.util.module_from_spec(spec); spec.loader.exec_module(repair)


class RepairTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(); self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.run = write_run(self.root / 'leap')
        self.raw = self.run / 'raw_annotation/current_annotation_1.gff3'
        self.extended = self.run / 'extended_transcripts/1_positive_extended.gff'
        self.stats = self.run / 'extended_transcripts_stats/1_positive_stats.tsv'
        # The production writer puts leap_tag on transcript rows, not children.
        for path in (self.extended, self.run / 'extended_transcripts/1_negative_extended.gff'):
            path.write_text(''.join(line if '\ttranscript\t' in line else line.replace(';leap_tag=LEAP_5,LEAP_3', '')
                                    for line in path.read_text().splitlines(True)))

    def invoke(self, expected=0, stage='leap', context=True):
        output = self.root / ('result_' + str(len(list(self.root.glob('result_*')))))
        command = [sys.executable, str(SCRIPT), str(self.raw), str(self.extended), str(output),
                   '--extension-stats', str(self.stats), '--stage', stage]
        if context: command += ['--leap-output-dir', str(self.run)]
        before = {p: p.read_bytes() for p in (self.raw, self.extended, self.stats)}
        result = subprocess.run(command, capture_output=True, text=True)
        self.assertEqual(result.returncode, expected, result.stdout + result.stderr)
        self.assertNotIn('Traceback', result.stderr)
        self.assertEqual({p: p.read_bytes() for p in before}, before)
        return output

    def test_healthy_models_are_unchanged_and_all_checks_pass(self):
        output = self.invoke()
        with (output / 'changes.csv').open() as handle: rows = list(csv.DictReader(handle))
        self.assertEqual(rows[0]['status'], 'unchanged')
        with (output / 'validation.csv').open() as handle:
            self.assertEqual([r['status'] for r in csv.DictReader(handle)], ['PASS'] * 3)

    def test_rebuilds_internal_exons_cds_codons_utrs_metadata_and_unknown_features(self):
        extra = '1\ttest\tcustom_feature\t220\t225\t4\t+\t.\tID=custom;Parent=t_plus.1;gene_id=g_plus;transcript_id=t_plus.1;note=retain\n'
        self.raw.write_text(self.raw.read_text() + extra)
        records = []
        for record in repair.read_records(self.extended):
            f = record.fields
            if f[2] in ('CDS', 'start_codon', 'stop_codon', 'three_prime_UTR'): continue
            if f[2] == 'exon' and f[3] == '200':
                record = repair.changed(record, 201, attrs=dict(record.attrs, exon_number='99', source_note='bad'))
            if f[2] == 'transcript': record = repair.changed(record, attrs=dict(record.attrs, transcript_name='wrong'))
            records.append(record)
        records.append(repair.Record(('1', 'bad', 'unwanted', '210', '215', '.', '+', '.'),
            dict(ID='extra', Parent='t_plus.1', gene_id='g_plus', transcript_id='t_plus.1')))
        repair.write_gff(self.extended, records)
        output = self.invoke()
        actual = list(repair.read_records(output / 'repaired.gff3'))
        self.assertIn('custom_feature', {r.fields[2] for r in actual})
        self.assertNotIn('unwanted', {r.fields[2] for r in actual})
        self.assertNotIn('source_note=bad', (output / 'repaired.gff3').read_text())
        transcript = next(r for r in actual if r.fields[2] == 'transcript')
        self.assertEqual(transcript.fields[3:5], ('90', '360'))
        self.assertEqual(json.loads((output / 'summary.json').read_text())['counts'], {'repaired': 1})

    def test_missing_terminal_exons_and_transcript_row_can_be_reconstructed(self):
        for missing_tx in (False, True):
            original = self.extended.read_text()
            records = [r for r in repair.read_records(self.extended)
                       if r.fields[2] != 'exon' and not (missing_tx and r.fields[2] == 'transcript')]
            repair.write_gff(self.extended, records)
            output = self.invoke()
            transcript = next(r for r in repair.read_records(output / 'repaired.gff3') if r.fields[2] == 'transcript')
            self.assertEqual(transcript.fields[3:5], ('90', '360'))
            self.extended.write_text(original)

    def test_both_strands_single_exon_multipart_codons_and_gzip(self):
        for strand, label, offset in (('+', 'positive', 0), ('-', 'negative', 400)):
            self.extended = self.run / 'extended_transcripts' / ('1_' + label + '_extended.gff')
            self.stats = self.run / 'extended_transcripts_stats' / ('1_' + label + '_stats.tsv')
            records = list(repair.read_records(self.extended))
            # Remove all CDS/codon records; raw holds the authoritative parts.
            repair.write_gff(self.extended, [r for r in records if r.fields[2] not in ('CDS', 'start_codon', 'stop_codon')])
            self.invoke()
        self.extended = self.run / 'extended_transcripts/1_positive_extended.gff'
        self.stats = self.run / 'extended_transcripts_stats/1_positive_stats.tsv'
        raw = annotation('t_plus.1', single=True).replace('g1', 'g_plus')
        raw += '1\ttest\tstart_codon\t128\t129\t.\t+\t1\tID=split_codon;Parent=t_plus.1;gene_id=g_plus;transcript_id=t_plus.1\n'
        self.raw.write_text(raw)
        self.extended.write_text(extend(raw))
        compressed = self.root / 'raw.gff3.gz'
        with gzip.open(compressed, 'wt') as handle: handle.write(self.raw.read_text())
        self.raw = compressed
        # Single-exon splice-facing boundaries differ from the multi-exon fixture.
        (self.run / 'fixed_extensions/1_positive_FE_protected.bed').write_text('1\t90\t350\t100\t0\t+\t.\t.\t7\tg_plus\n')
        (self.run / 'fixed_extensions/1_positive_LE_protected.bed').write_text('1\t100\t360\t350\t0\t+\t.\t.\t7\tg_plus\n')
        output = self.invoke()
        self.assertEqual(sum(r.fields[2] == 'start_codon' for r in repair.read_records(output / 'repaired.gff3')), 2)

    def test_stats_and_tags_are_recomputed_from_verified_endpoints(self):
        with self.stats.open() as handle: reader = csv.DictReader(handle, delimiter='\t'); fields = reader.fieldnames; rows = list(reader)
        rows[0].update(five_prime_extension='999', five_prime_support='current_annotation', five_prime_read_count='0')
        repair.write_csv(self.stats, fields, rows, '\t')
        output = self.invoke()
        with (output / 'repaired_stats.tsv').open() as handle: rows = list(csv.DictReader(handle, delimiter='\t'))
        self.assertEqual(rows[0]['five_prime_extension'], '10')
        self.assertEqual(rows[0]['five_prime_support'], 'reads')
        self.assertEqual(rows[0]['five_prime_read_count'], '7')

    def test_processed_restored_end_retains_blacklist_and_other_end_support(self):
        records = []
        for r in repair.read_records(self.extended):
            if r.fields[2] in ('transcript', 'exon', 'five_prime_UTR') and r.fields[3] == '90': r = repair.changed(r, 100)
            if r.fields[2] == 'transcript': r = repair.changed(r, attrs=dict(r.attrs, leap_tag='LEAP_blacklisted,LEAP_3', transcript_class='internal_left'))
            if r.fields[2] != 'gene': records.append(r)
        repair.write_gff(self.extended, records)
        output = self.invoke(stage='processed')
        self.assertFalse(any(r.fields[2] == 'gene' for r in repair.read_records(output / 'repaired.gff3')))
        with (output / 'repaired_stats.tsv').open() as handle: row = next(csv.DictReader(handle, delimiter='\t'))
        self.assertEqual(row['five_prime_support'], 'blacklisted')
        self.assertEqual(row['five_prime_read_count'], 'blacklisted')
        self.assertEqual(row['three_prime_read_count'], '7')
        self.assertEqual(row['five_prime_extension'], '0')

    def test_unverified_inward_or_ambiguous_boundaries_stop_without_fixed_outputs(self):
        output = self.invoke(expected=1, context=False)
        self.assertFalse((output / 'repaired.gff3').exists())
        self.assertTrue((output / 'changes.csv').is_file())
        original = self.extended.read_text()
        records = list(repair.read_records(self.extended))
        transcript = next(r for r in records if r.fields[2] == 'transcript')
        repair.write_gff(self.extended, records + [transcript])
        output = self.invoke(expected=1)
        self.assertFalse((output / 'repaired_stats.tsv').exists())
        self.extended.write_text(original)
        repair.write_gff(self.extended, [repair.changed(r, 101) if r.fields[2] == 'transcript' else r for r in records])
        self.invoke(expected=1)

    def test_last_candidate_is_authoritative_and_never_falls_back(self):
        path = self.run / 'fixed_extensions/1_positive_FE_protected.bed'
        path.write_text(path.read_text() + '1\t95\t150\t100\t0\t+\t.\t.\t2\tg_plus\n')
        output = self.invoke(expected=1)
        self.assertIn('last candidate', (output / 'changes.csv').read_text())

    def test_no_growth_needs_no_context_and_clears_support_without_restoring_excluded_models(self):
        records = [r for r in repair.read_records(self.raw)
                   if r.fields[6] == '+']
        repair.write_gff(self.extended, records)
        output = self.invoke(context=False)
        with (output / 'repaired_stats.tsv').open() as handle:
            row = next(csv.DictReader(handle, delimiter='\t'))
        for side in ('five_prime', 'three_prime'):
            self.assertEqual(row[side + '_extension'], '0')
            self.assertEqual(row[side + '_support'], 'none')
            self.assertEqual(row[side + '_read_count'], '0')
        self.assertNotIn('t_minus.1', (output / 'repaired.gff3').read_text())

    def test_processed_contained_endpoint_and_annotation_support(self):
        for path in self.run.glob('fixed_extensions/*positive*'):
            path.write_text(path.read_text().replace('\t100\t0\t', '\t.\t0\t')
                            .replace('\t350\t0\t', '\t.\t0\t').replace('\t7\t', '\t0\t'))
        records = [repair.changed(r, 95) if r.fields[3] == '90' else r
                   for r in repair.read_records(self.extended)]
        repair.write_gff(self.extended, records)
        self.invoke(expected=1)
        output = self.invoke(stage='processed')
        with (output / 'repaired_stats.tsv').open() as handle:
            row = next(csv.DictReader(handle, delimiter='\t'))
        self.assertEqual(row['five_prime_extension'], '5')
        self.assertEqual(row['five_prime_support'], 'current_annotation')
        self.assertEqual(row['five_prime_read_count'], '0')

    def test_conflicting_class_and_mismatched_membership_are_not_repaired(self):
        records = [repair.changed(r, attrs=dict(r.attrs, transcript_class=(
            'external' if r.fields[2] == 'transcript' else 'internal'))) for r in repair.read_records(self.extended)]
        repair.write_gff(self.extended, records)
        output = self.invoke(expected=1)
        self.assertIn('Conflicting transcript_class', (output / 'changes.csv').read_text())
        records = [r for r in records if r.fields[2] == 'gene']
        repair.write_gff(self.extended, records)
        self.invoke(expected=2)

    def test_ambiguous_context_and_malformed_inputs_produce_diagnostics(self):
        path = self.run / 'fixed_extensions/1_positive_FE_protected.bed'
        duplicate = path.with_name('other_1_positive_FE_protected.bed')
        duplicate.write_bytes(path.read_bytes())
        output = self.invoke(expected=2)
        self.assertIn('exactly one', (output / 'summary.json').read_text())
        duplicate.unlink()
        self.raw.write_text(self.raw.read_text() + 'not a GFF record\n')
        output = self.invoke(expected=2)
        self.assertIn('nine GFF3 columns', (output / 'summary.json').read_text())

    def test_membership_version_raw_errors_and_overwrite_protection(self):
        original = self.extended.read_text()
        self.extended.write_text(original.replace('t_plus.1', 't_plus.2'))
        output = self.invoke(expected=2)
        self.assertFalse((output / 'repaired.gff3').exists())
        self.extended.write_text(original)
        raw_records = list(repair.read_records(self.raw))
        exon = next(r for r in raw_records if r.fields[2] == 'exon')
        repair.write_gff(self.raw, raw_records + [exon])
        self.invoke(expected=2)
        result = subprocess.run([sys.executable, str(SCRIPT), str(self.raw), str(self.extended), str(self.run),
            '--extension-stats', str(self.stats)], capture_output=True, text=True)
        self.assertEqual(result.returncode, 2)
        self.assertIn('new directory', result.stderr)

    def test_par_y_ids_and_gene_relationships_are_preserved(self):
        for path in (self.raw, self.extended):
            records = []
            for r in repair.read_records(path):
                attrs = dict(r.attrs)
                if attrs.get('ID') == 'g_plus': attrs['ID'] = 'g_plus_PAR_Y'
                if attrs.get('Parent') == 'g_plus': attrs['Parent'] = 'g_plus_PAR_Y'
                records.append(repair.changed(r, attrs=attrs))
            repair.write_gff(path, records)
        output = self.invoke()
        gene = next(r for r in repair.read_records(output / 'repaired.gff3') if r.fields[2] == 'gene')
        self.assertEqual(gene.attrs['ID'], 'g_plus_PAR_Y')
        self.assertEqual(gene.attrs['gene_id'], 'g_plus')


if __name__ == '__main__': unittest.main()
