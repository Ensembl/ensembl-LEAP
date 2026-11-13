import csv
import gzip
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / 'assets/check_reference_coordinates.py'


def feature(chrom='1', start=1, end=100, kind='exon', attrs='ID=f;transcript_id=t.1;gene_id=g'):
    return '{}\ttest\t{}\t{}\t{}\t.\t+\t.\t{}\n'.format(chrom, kind, start, end, attrs)


class ReferenceCoordinatesTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.annotation = self.root / 'annotation.gff'
        self.reference = self.root / 'reference.sizes'
        self.output = self.root / 'report.csv'

    def put(self, path, value):
        with (gzip.open if str(path).endswith('.gz') else open)(path, 'wt') as handle:
            handle.write(value)

    def check(self, annotation=None, reference=None):
        if annotation is not None:
            self.put(self.annotation, annotation)
        if reference is not None:
            self.put(self.reference, reference)
        result = subprocess.run([sys.executable, str(SCRIPT), str(self.annotation), str(self.reference), str(self.output)], capture_output=True, text=True)
        rows = []
        if result.returncode != 2 and self.output.exists():
            with self.output.open() as handle:
                rows = list(csv.DictReader(handle))
        return result, rows

    def test_boundaries_all_feature_types_and_fai(self):
        for reference in ('1\t100\n', '1\t100\t3\t60\t61\n'):
            result, rows = self.check(''.join(feature(kind=kind) for kind in ('gene','transcript','exon','CDS','unfamiliar')), reference)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(rows, [])
            self.assertIn('features: 5', result.stdout)

    def test_every_coordinate_violation_and_exact_names(self):
        for kwargs, issue in (({'chrom': 'chr1'}, 'unknown_chromosome'),
                              ({'start': 0}, 'start_before_reference'),
                              ({'start': -2, 'end': -1}, 'end_before_reference'),
                              ({'start': 90, 'end': 80}, 'reversed_interval'),
                              ({'end': 101}, 'beyond_reference_end'),
                              ({'start': 101, 'end': 110}, 'beyond_reference_end')):
            with self.subTest(kwargs=kwargs):
                result, rows = self.check(feature(**kwargs), '1\t100\n')
                self.assertEqual(result.returncode, 1, result.stderr)
                self.assertIn(issue, rows[0]['status'])
                self.assertEqual(rows[0]['line'], '1')
                self.assertEqual(rows[0]['transcript_id'], 't.1')
                self.assertIn('violating_features: 1', result.stdout)
        result, rows = self.check(feature(chrom='Y', end=101), '1\t100\nY\t100\n')
        self.assertEqual(rows[0]['reference_length'], '100')

    def test_empty_header_only_and_embedded_fasta(self):
        for value in ('', '##gff-version 3\n', feature() + '##FASTA\n>1\nACGT\n'):
            result, rows = self.check(value, '1\t100\n')
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(rows, [])
        result, _ = self.check(feature(attrs='.'), '1\t100\n')
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_directory_gzip_and_ambiguity(self):
        self.annotation = self.root / 'gffs'
        self.annotation.mkdir()
        self.reference = self.root / 'reference.fai.gz'
        self.put(self.annotation / '1_positive.gff.gz', feature())
        self.put(self.annotation / '1_negative.gff3', feature(end=101))
        result, rows = self.check(reference='1\t100\t3\t60\t61\n')
        self.assertEqual(result.returncode, 1, result.stderr)
        self.assertEqual(len(rows), 1)
        self.assertIn('files: 2', result.stdout)
        self.put(self.annotation / '1_positive.gff', feature())
        result, _ = self.check()
        self.assertEqual(result.returncode, 2)
        self.assertIn('ambiguous', result.stderr)

    def test_invalid_lengths(self):
        for reference in ('', '1\t0\n', '1\t-1\n', '1\tbad\n', '1\t100\n1\t100\n',
                          '1\t100\t3\n', '1\t100\t3\t0\t1\n', '1\t100\t3\t60\t59\n',
                          '1\t100\n2\t100\t3\t60\t61\n'):
            result, _ = self.check(feature(), reference)
            self.assertEqual(result.returncode, 2, result.stderr)

    def test_malformed_gff_preserves_existing_report(self):
        self.output.write_text('existing report\n')
        for invalid in ('1\tbad\n', feature(start='x'), feature().replace('ID=f;', 'ID=f;ID=g;'),
                        feature().replace('\t+\t', '\tinvalid\t'), feature().replace('\t.\t+\t', '\tbad\t+\t')):
            result, _ = self.check(feature(end=101) + invalid, '1\t100\n')
            self.assertEqual(result.returncode, 2, result.stderr)
            self.assertEqual(self.output.read_text(), 'existing report\n')

    def test_missing_files_corrupt_gzip_and_empty_directory(self):
        result, _ = self.check(reference='1\t100\n')
        self.assertEqual(result.returncode, 2)
        self.annotation = self.root / 'empty'
        self.annotation.mkdir()
        result, _ = self.check()
        self.assertEqual(result.returncode, 2)
        self.annotation = self.root / 'bad.gff.gz'
        self.annotation.write_bytes(bytes.fromhex('1f8b080000000000000007'))
        result, _ = self.check()
        self.assertEqual(result.returncode, 2)
        self.assertNotIn('Traceback', result.stderr)

    def test_overwrite_protection_including_links(self):
        self.put(self.annotation, feature())
        self.put(self.reference, '1\t100\n')
        for target in (self.annotation, self.reference):
            for kind in ('direct', 'symlink', 'hardlink'):
                self.output = self.root / ('report-' + kind)
                self.output.unlink(missing_ok=True)
                if kind == 'direct':
                    self.output = target
                elif kind == 'symlink':
                    self.output.symlink_to(target)
                else:
                    os.link(target, self.output)
                result, _ = self.check()
                self.assertEqual(result.returncode, 2, result.stderr)
                self.assertIn('overwrite', result.stderr)
