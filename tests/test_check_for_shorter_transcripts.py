"""Run with python3 -m unittest discover -s tests."""
import csv
import gzip
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / 'assets' / 'check_for_shorter_transcripts.py'


def transcript(identifier='t1', start=100, end=200, strand='+', chromosome='1'):
    return '{}\ttest\ttranscript\t{}\t{}\t.\t{}\t.\ttranscript_id={};\n'.format(
        chromosome, start, end, strand, identifier)


class BoundaryCheckTests(unittest.TestCase):
    def run_check(self, raw, extended, compressed=False):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            inputs = []
            for name, content in [('raw', raw), ('extended', extended)]:
                path = root / (name + ('.gff3.gz' if compressed else '.gff3'))
                opener = gzip.open if compressed else open
                with opener(path, 'wt') as handle:
                    handle.write(content)
                inputs.append(str(path))
            report = root / 'report.csv'
            result = subprocess.run([sys.executable, str(SCRIPT)] + inputs + [str(report)],
                                    capture_output=True, text=True)
            if report.exists():
                with report.open(newline='') as handle:
                    rows = list(csv.DictReader(handle))
            else:
                rows = None
            return result, rows

    def test_boundaries_on_both_strands(self):
        cases = [(100, 200, 'unchanged', False, False),
                 (90, 210, 'extended', False, False),
                 (90, 200, 'extended', False, False),
                 (100, 210, 'extended', False, False),
                 (110, 200, 'shortened', True, False),
                 (100, 190, 'shortened', False, True),
                 (110, 190, 'shortened', True, True),
                 (110, 300, 'shortened', True, False),
                 (1, 190, 'shortened', False, True)]
        for strand in ('+', '-'):
            for start, end, status, left, right in cases:
                with self.subTest(strand=strand, start=start, end=end):
                    result, rows = self.run_check(transcript(strand=strand),
                                                 transcript(start=start, end=end, strand=strand))
                    self.assertEqual(result.returncode, int(status == 'shortened'), result.stderr)
                    self.assertEqual(rows[0]['status'], status)
                    self.assertEqual(rows[0]['shorter_left'], str(left))
                    self.assertEqual(rows[0]['shorter_right'], str(right))
                    self.assertEqual(rows[0]['length_current'], '101')
                    self.assertEqual(rows[0]['length_extended'], str(end - start + 1))
                    self.assertIn('{}: 1'.format(status), result.stdout)

    def test_unmatched_exact_ids_and_gzip(self):
        result, rows = self.run_check('# header\n\n' + transcript('t1.1') + transcript('t2'),
                                      transcript('t1.2') + transcript('t2'), compressed=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual([row['status'] for row in rows],
                         ['missing_from_extended', 'missing_from_raw', 'unchanged'])
        self.assertEqual(rows[0]['start_extended'], '')
        self.assertEqual(rows[1]['start_current'], '')
        self.assertIn('PASS: 0 failing transcript(s)', result.stdout)
        self.assertIn('missing_from_extended: 1', result.stdout)
        self.assertIn('missing_from_raw: 1', result.stdout)

    def test_location_mismatches(self):
        for kwargs, flag in [({'chromosome': '2'}, 'chromosome_mismatch'),
                             ({'strand': '-'}, 'strand_mismatch')]:
            with self.subTest(flag=flag):
                result, rows = self.run_check(transcript(), transcript(**kwargs))
                self.assertEqual(result.returncode, 1)
                self.assertEqual(rows[0]['status'], 'location_mismatch')
                self.assertEqual(rows[0][flag], 'True')
                self.assertEqual(rows[0]['shorter_left'], '')
                self.assertIn('location_mismatch: 1', result.stdout)

    def test_invalid_inputs(self):
        cases = [('', 'no transcript records'),
                 ('# comment\n', 'no transcript records'),
                 (transcript().replace('transcript_id=t1;', 'gene_id=g1;'), 'missing transcript_id'),
                 (transcript() * 2, 'duplicate transcript_id'),
                 (transcript(start=0), 'invalid transcript coordinates'),
                 (transcript(start=201), 'invalid transcript coordinates'),
                 (transcript(start='x'), 'must be integers'),
                 (transcript().replace('transcript_id=t1;', 'broken'), 'malformed'),
                 (transcript().replace('transcript_id=t1;', 'transcript_id=t1;transcript_id=t2;'), 'duplicate attribute'),
                 ('1\ttest\ttranscript\n', '9 GFF3 columns'),
                 (transcript(strand='x'), 'invalid chromosome or strand')]
        for invalid, message in cases:
            for raw, new in [(invalid, transcript()), (transcript(), invalid)]:
                with self.subTest(message=message, raw=raw):
                    result, rows = self.run_check(raw, new)
                    self.assertEqual(result.returncode, 2)
                    self.assertIn(message, result.stderr)
                    self.assertIsNone(rows)

    def test_output_errors_and_input_protection(self):
        with tempfile.TemporaryDirectory() as directory:
            raw = Path(directory) / 'raw.gff3'
            raw.write_text(transcript())
            for output in (str(raw), directory + '/missing/report.csv'):
                result = subprocess.run([sys.executable, str(SCRIPT), str(raw), str(raw), output],
                                        capture_output=True, text=True)
                self.assertEqual(result.returncode, 2)
                self.assertIn('ERROR:', result.stderr)
                self.assertEqual(raw.read_text(), transcript())

    def test_missing_file(self):
        with tempfile.TemporaryDirectory() as directory:
            result = subprocess.run([sys.executable, str(SCRIPT), directory + '/missing',
                                     directory + '/also_missing', directory + '/out.csv'],
                                    capture_output=True, text=True)
            self.assertEqual(result.returncode, 2)
            self.assertIn('ERROR:', result.stderr)


if __name__ == '__main__':
    unittest.main()
