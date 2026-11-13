import csv
import gzip
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from test_partial_transcript_processing import model

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / 'assets/check_partial_protection.py'
PRODUCTION = ROOT / 'bin/protect_partial_extensions.py'


class PartialProtectionAuditTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.output = self.root / 'report.csv'
        for folder in ('raw_annotation','fixed_extensions_before_partial_protection','fixed_extensions','partial_protection','natural_overlaps'):
            (self.root / folder).mkdir()

    def setup_partition(self, partial='', full=None, strand='+', natural='', fe=None, le=None):
        full = model('source.1', strand=strand) if full is None else full
        raw = self.root / 'raw_annotation/current_annotation_1.gff3'
        raw.write_text(full + partial)
        label = 'positive' if strand == '+' else 'negative'
        originals = []
        for kind, value in (('FE', fe), ('LE', le)):
            left = (kind == 'FE') == (strand == '+')
            if value is None:
                value = ('1\t150\t250\t185/200\t1\t{}\t.\t.\t3\tA\n' if left else
                         '1\t400\t500\t470/450\t1\t{}\t.\t.\t4\tA\n').format(strand)
            path = self.root / 'fixed_extensions_before_partial_protection' / ('test_1_' + label + '_' + kind + '_SCF.bed')
            path.write_text(value)
            originals.append(path)
        pairs = self.root / 'natural_overlaps' / ('1_' + label + '_NO.txt')
        pairs.write_text(natural)
        # Fixture generation uses production only to exercise real output formats;
        # endpoint expectations in the tests below are asserted independently.
        full_path, partial_path = self.root / 'full.gff', self.root / 'partial.gff'
        full_path.write_text(full)
        partial_path.write_text(partial)
        prefix = self.root / 'generated'
        result = subprocess.run([sys.executable, str(PRODUCTION), *map(str, originals), str(full_path), str(partial_path), str(pairs), str(prefix)], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        after = []
        for kind in ('FE','LE'):
            path = self.root / 'fixed_extensions' / ('1_' + label + '_' + kind + '_protected.bed')
            (self.root / ('generated_' + kind + '_protected.bed')).replace(path)
            after.append(path)
        audit = self.root / 'partial_protection' / ('1_' + label + '_partial_protection.tsv')
        (self.root / 'generated_partial_protection.tsv').replace(audit)
        return originals, after, audit, label

    def check(self, label='positive', extra=()):
        result = subprocess.run([sys.executable, str(SCRIPT), str(self.root), str(self.output), '--chromosomes', '1', '--strand', label] + list(extra), capture_output=True, text=True)
        rows = []
        if result.returncode != 2 and self.output.exists():
            with self.output.open() as handle:
                rows = list(csv.DictReader(handle))
        return result, rows

    def test_empty_partials_unchanged_and_annotation_context(self):
        before, after, _, label = self.setup_partition(fe='1\t150\t250\t.\t0\t+\t.\t.\t3\tA\n',
                    le='1\t400\t450\tsingle_exon\t0\t+\tB\t.\t.\t0\n')
        result, rows = self.check(label)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual([r['status'] for r in rows], ['unchanged','unchanged'])
        self.assertIn('candidate_records: 2', result.stdout)
        after[0].write_text(after[0].read_text().replace('\t150\t', '\t140\t'))
        result, rows = self.check(label)
        self.assertEqual(result.returncode, 1)
        self.assertEqual(rows[0]['lengthened_endpoint'], 'True')

    def test_shorter_and_fallback_both_ends_both_strands(self):
        for strand in ('+','-'):
            for kind in ('FE','LE'):
                for fallback in (False, True):
                    with self.subTest(strand=strand, kind=kind, fallback=fallback):
                        # reset generated files, keeping the temporary root
                        for folder in self.root.iterdir():
                            if folder.is_dir():
                                for p in folder.iterdir(): p.unlink()
                        left = (kind == 'FE') == (strand == '+')
                        exons = ((100,195 if fallback else 170),(300,350)) if left else ((480,520),(600,650))
                        cds = ((160,195 if fallback else 170),) if left else ((450 if fallback else 485,495),)
                        if not left and fallback:
                            exons = ((450,520),(600,650))
                        partial = model('partial.1', strand=strand, exons=exons, cds=cds, five=False, three=False)
                        _, _, _, label = self.setup_partition(partial=partial, strand=strand)
                        result, rows = self.check(label)
                        self.assertEqual(result.returncode, 0, result.stderr)
                        row = next(r for r in rows if r['end'] == kind)
                        self.assertEqual(row['status'], 'annotation_fallback' if fallback else 'shorter_read_candidate')
                        self.assertEqual(row['expected_endpoint'], ('200' if left else '450') if fallback else ('185' if left else '470'))
                        self.assertIn('partial_cds:partial.1', row['blocker_evidence'])

    def test_growth_invented_endpoints_and_metadata_changes(self):
        _, after, _, label = self.setup_partition()
        original = after[0].read_text()
        for before, replacement, flag in (('\t150\t', '\t140\t', 'lengthened_endpoint'),
                    ('\t150\t', '\t175\t', 'invented_endpoint'),
                    ('\t250\t', '\t249\t', 'protected_fields_changed'),
                    ('\t3\tA\n', '\t4\tA\n', 'protected_fields_changed'),
                    ('185/200', '200', 'alternatives_mismatch'),
                    ('\tA\n','\tA,B\n','protected_fields_changed')):
            after[0].write_text(original.replace(before, replacement))
            result, rows = self.check(label)
            self.assertEqual(result.returncode, 1, result.stderr)
            self.assertEqual(rows[0][flag], 'True')
        after[0].write_bytes(original.replace('\n','\r\n').encode())
        result, rows = self.check(label)
        self.assertEqual(result.returncode, 1)
        self.assertEqual(rows[0]['unchanged_record_changed'], 'True')

    def test_duplicate_multiplicity_and_order(self):
        row = '1\t150\t250\t185/200\t1\t+\t.\t.\t3\tA\n'
        _, after, _, label = self.setup_partition(fe=row + row)
        result, rows = self.check(label)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(len(rows), 3)
        after[0].write_text(row)
        result, rows = self.check(label)
        self.assertEqual(result.returncode, 1)
        self.assertTrue(any(r['missing_candidate'] == 'True' for r in rows))
        after[0].write_text(row * 3)
        result, rows = self.check(label)
        self.assertEqual(result.returncode, 1)
        self.assertTrue(any(r['unexpected_candidate'] == 'True' for r in rows))
        before = next((self.root/'fixed_extensions_before_partial_protection').glob('*_FE_*'))
        other = row.replace('\t150\t','\t160\t')
        before.write_text(row + other)
        after[0].write_text(other + row)
        result, rows = self.check(label)
        self.assertEqual(result.returncode, 1)
        self.assertTrue(any(r['candidate_order_mismatch'] == 'True' for r in rows))

    def test_audit_missing_wrong_extra_and_multiplicity(self):
        partial = model('partial.1', exons=((100,170),(300,350)), cds=((160,170),), five=False, three=False)
        _, _, audit, label = self.setup_partition(partial=partial)
        original = audit.read_text()
        header, data = original.splitlines(True)
        for content, status in ((header,'missing_audit'), (original + data,'unexpected_audit'),
                                (original.replace('partial_cds:partial.1', 'partial_cds:wrong'), 'missing_audit')):
            audit.write_text(content)
            result, rows = self.check(label)
            self.assertEqual(result.returncode, 1, result.stderr)
            self.assertTrue(any(r['status'] == status for r in rows))
        audit.write_text(original)
        before = next((self.root/'fixed_extensions_before_partial_protection').glob('*_FE_*'))
        before.write_text(before.read_text()*2)
        after = self.root/'fixed_extensions/1_positive_FE_protected.bed'
        after.write_text(after.read_text()*2)
        audit.write_text(original + data)
        result, rows = self.check(label)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(sum(r['status']=='shorter_read_candidate' for r in rows), 2)

    def test_terminal_and_single_exon_contacts_allowances(self):
        # Non-CDS terminal obstacles use ordered pairs and positive-length contact.
        for partial, pairs, expected in (
            (model('p', 'B', exons=((50,100),(160,180)), cds=(), five=False, three=False), '', 'shorter_read_candidate'),
            (model('p', 'B', exons=((50,100),(160,180)), cds=(), five=False, three=False), 'A\tB\n', 'unchanged'),
            (model('p', 'B', exons=((50,100),(160,180)), cds=(), five=False, three=False), 'B\tA\n', 'shorter_read_candidate'),
            (model('p', 'A', exons=((50,100),(160,180)), cds=(), five=False, three=False), '', 'unchanged'),
            (model('p', 'B', exons=((160,180),), cds=(), five=False, three=False), '', 'shorter_read_candidate'),
            (model('p', 'B', exons=((100,150),), cds=(), five=False, three=False), '', 'unchanged')):
            for folder in self.root.iterdir():
                if folder.is_dir():
                    for p in folder.iterdir(): p.unlink()
            _, _, _, label = self.setup_partition(partial=partial, natural=pairs)
            result, rows = self.check(label)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(rows[0]['status'], expected)

    def test_special_flank_protection_and_boundary_touch(self):
        # Gene A has one full-length CDS and is in a natural-overlap pair.
        full = model('source.1', cds=((220,240),))
        partial = model('partial.1', exons=((160,170),(300,350)), cds=((320,330),), five=False, three=False)
        _, _, _, label = self.setup_partition(partial=partial, full=full, natural='A\tB\n')
        result, rows = self.check(label)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('partial_single_cds_flank', rows[0]['blocker_evidence'])
        # CDS contact at the baseline is inclusive and forces fallback.
        raw = self.root/'raw_annotation/current_annotation_1.gff3'
        raw.write_text(model('source.1') + model('partial.1', exons=((200,210),), cds=((200,205),), five=False, three=False))
        result, rows = self.check(label)
        self.assertEqual(result.returncode, 1)
        self.assertEqual(rows[0]['expected_endpoint'], '200')

    def test_stale_longer_options_annotation_fallback_and_original_groups(self):
        partial = model('p', exons=((100,195),(300,350)), cds=((160,195),), five=False, three=False)
        before, after, _, label = self.setup_partition(partial=partial,
                fe='1\t150\t250\t100/185/200\t2\t+\t.\t.\t3\tA\n')
        result, rows = self.check(label)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(rows[0]['expected_endpoint'], '200')
        self.assertEqual(rows[0]['actual_candidate'].split('\t')[8], '3')
        after[0].write_text(before[0].read_text().replace('\t150\t','\t100\t'))
        result, rows = self.check(label)
        self.assertEqual(result.returncode, 1)
        self.assertEqual(rows[0]['lengthened_endpoint'], 'True')

    def test_gzip_and_audit_record_order(self):
        partial = model('p', exons=((100,170),(480,520)), cds=((160,170),(485,495)), five=False, three=False)
        _, _, audit, label = self.setup_partition(partial=partial)
        lines = audit.read_text().splitlines(True)
        audit.write_text(lines[0] + ''.join(reversed(lines[1:])))
        for folder in self.root.iterdir():
            if folder.is_dir():
                for path in list(folder.iterdir()):
                    with gzip.open(str(path)+'.gz','wb') as handle: handle.write(path.read_bytes())
                    path.unlink()
        result, rows = self.check(label)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(sum(r['status']=='shorter_read_candidate' for r in rows), 2)

    def test_context_errors_and_invalid_records(self):
        _, after, audit, label = self.setup_partition()
        original = after[0].read_text()
        self.output.write_text('existing\n')
        for invalid in ('bad\n', original.replace('\t150\t','\tbad\t'), original.replace('\t+\t','\t-\t')):
            after[0].write_text(invalid)
            result, _ = self.check(label)
            self.assertEqual(result.returncode, 2, result.stderr)
            self.assertEqual(self.output.read_text(), 'existing\n')
        after[0].write_text(original)
        audit.write_text('wrong header\n')
        result, _ = self.check(label)
        self.assertEqual(result.returncode, 2)
        audit.unlink()
        result, _ = self.check(label)
        self.assertEqual(result.returncode, 2)
        self.assertIn('expected exactly one', result.stderr)

    def test_ambiguous_inputs_and_overwrite_links(self):
        before, after, _, label = self.setup_partition()
        copy = Path(str(before[0])+'.gz')
        with gzip.open(copy,'wb') as handle: handle.write(before[0].read_bytes())
        result, _ = self.check(label)
        self.assertEqual(result.returncode, 2)
        copy.unlink()
        for target in (before[0],after[0],self.root/'raw_annotation/current_annotation_1.gff3'):
            for kind in ('direct','symlink','hardlink'):
                self.output = self.root/('report-'+kind)
                self.output.unlink(missing_ok=True)
                if kind == 'direct': self.output = target
                elif kind == 'symlink': self.output.symlink_to(target)
                else: os.link(target, self.output)
                result, _ = self.check(label)
                self.assertEqual(result.returncode, 2, result.stderr)
                self.assertIn('overwrite', result.stderr)

    def test_ordered_gene_groups_and_context_separation(self):
        partial = model('p', 'B', exons=((50,100),(160,180)), cds=(), five=False, three=False)
        cases = [('A,C', 'C\tB\n', partial, 'shorter_read_candidate'),
                 ('A,C', 'A\tB\n', partial, 'unchanged'),
                 ('C,A', 'A\tB\n', partial, 'shorter_read_candidate'),
                 ('A,B', '', partial, 'unchanged'),
                 ('A', '', partial.replace('1\ttest', '2\ttest'), 'unchanged'),
                 ('A', '', partial.replace('\t+\t', '\t-\t'), 'unchanged')]
        for genes, pairs, context, expected in cases:
            for folder in self.root.iterdir():
                if folder.is_dir():
                    for p in folder.iterdir(): p.unlink()
            _, _, _, label = self.setup_partition(partial=context, natural=pairs,
                fe='1\t150\t250\t185/200\t1\t+\t.\t.\t3\t' + genes + '\n')
            result, rows = self.check(label)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(rows[0]['status'], expected)

    def test_annotation_support_and_nonoutward_selection_are_preserved(self):
        partial = model('p', exons=((100,210),(300,350)), cds=((150,205),), five=False, three=False)
        for candidate in ('1\t150\t250\t.\t0\t+\t.\t.\t3\tA\n',
                          '1\t200\t250\t150/200\t1\t+\t.\t.\t3\tA\n'):
            for folder in self.root.iterdir():
                if folder.is_dir():
                    for p in folder.iterdir(): p.unlink()
            _, _, _, label = self.setup_partition(partial=partial, fe=candidate)
            result, rows = self.check(label)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(rows[0]['status'], 'unchanged')

    def test_malformed_raw_and_natural_context_and_missing_output(self):
        _, after, _, label = self.setup_partition()
        raw = self.root/'raw_annotation/current_annotation_1.gff3'
        original = raw.read_text()
        for content in (original + original.splitlines(True)[0],
                        original.replace('transcript_type=protein_coding;', ''),
                        original.replace('\texon\t', '\tother\t')):
            raw.write_text(content)
            result, _ = self.check(label)
            self.assertEqual(result.returncode, 2, result.stderr)
        raw.write_text(original)
        pairs = self.root/'natural_overlaps/1_positive_NO.txt'
        pairs.write_text('only_one_gene\n')
        result, _ = self.check(label)
        self.assertEqual(result.returncode, 2)
        pairs.write_text('')
        after[0].unlink()
        result, _ = self.check(label)
        self.assertEqual(result.returncode, 2)
        self.assertIn('expected exactly one', result.stderr)
