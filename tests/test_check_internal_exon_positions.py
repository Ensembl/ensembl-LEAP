import csv
import gzip
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / 'assets/check_internal_exon_positions.py'
EXONS = [(100, 150), (200, 250), (300, 350)]


def annotation(exons=None, identifier='t1', chromosome='1', strand='+', reverse=False, exon_id='old'):
    exons = EXONS if exons is None else exons
    def feature(kind, start, end, extra=''):
        return '{}\ttest\t{}\t{}\t{}\t.\t{}\t.\ttranscript_id={};{}\n'.format(
            chromosome, kind, start, end, strand, identifier, extra)
    rows = [feature('transcript', 1, 1000)] + [
        feature('exon', start, end, 'exon_id={};'.format(exon_id)) for start, end in exons]
    return ''.join(reversed(rows) if reverse else rows)


class InternalExonTests(unittest.TestCase):
    def run_check(self, raw, extended, compressed=False):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            paths = []
            for name, content in [('raw', raw), ('extended', extended)]:
                path = root / (name + ('.gff3.gz' if compressed else '.gff3'))
                with (gzip.open if compressed else open)(path, 'wt') as handle:
                    handle.write(content)
                paths.append(str(path))
            output = root / 'report.csv'
            result = subprocess.run([sys.executable, str(SCRIPT)] + paths + [str(output)],
                                    capture_output=True, text=True)
            rows = None
            if output.exists():
                with output.open(newline='') as handle:
                    rows = list(csv.DictReader(handle))
            return result, rows

    def test_unchanged_and_outer_boundaries(self):
        for strand in ('+', '-'):
            for exons in [EXONS, [(90, 150), (200, 250), (300, 400)],
                          [(110, 150), (200, 250), (300, 340)]]:
                with self.subTest(strand=strand, exons=exons):
                    result, rows = self.run_check(annotation(strand=strand),
                        annotation(exons, strand=strand, reverse=True, exon_id='new'), compressed=True)
                    self.assertEqual(result.returncode, 0, result.stderr)
                    self.assertEqual(rows[0]['status'], 'unchanged_internal_structure')
                    self.assertEqual(rows[0]['exon_count_current'], '3')
                    self.assertEqual(rows[0]['exon_coordinates_extended'],
                                     ';'.join('{}-{}'.format(*exon) for exon in exons))
                    self.assertIn('PASS: 0 failing transcript(s)', result.stdout)

    def test_each_protected_boundary(self):
        for strand in ('+', '-'):
            for index, offset in [(0, 1), (1, 0), (1, 1), (2, 0)]:
                exons = [list(exon) for exon in EXONS]
                exons[index][offset] += 1
                with self.subTest(strand=strand, index=index, offset=offset):
                    result, rows = self.run_check(annotation(strand=strand), annotation(exons, strand=strand))
                    self.assertEqual(result.returncode, 1, result.stderr)
                    self.assertEqual(rows[0]['status'], 'changed_boundaries')
                    detail = 'exon_{}_{}:{}->{}'.format(index + 1, ('start', 'end')[offset],
                                                       EXONS[index][offset], exons[index][offset])
                    self.assertEqual(rows[0]['changed_boundary_details'], detail)
                    self.assertIn('flag_changed_boundaries: 1', result.stdout)

    def test_count_changes_and_combined_flags(self):
        for exons in [EXONS[:2], EXONS + [(400, 450)]]:
            result, rows = self.run_check(annotation(), annotation(exons))
            self.assertEqual(result.returncode, 1)
            self.assertEqual(rows[0]['status'], 'exon_count_mismatch')
            self.assertEqual(rows[0]['changed_boundaries'], '')
        result, rows = self.run_check(annotation(), annotation(EXONS[:2], chromosome='2', strand='-'))
        self.assertEqual(result.returncode, 1)
        for flag in ('chromosome_mismatch', 'strand_mismatch', 'exon_count_mismatch'):
            self.assertEqual(rows[0][flag], 'True')
        self.assertEqual(rows[0]['status'], 'location_mismatch')
        result, rows = self.run_check(annotation(), annotation([(100, 151)] + EXONS[1:], chromosome='2'))
        self.assertEqual(result.returncode, 1)
        self.assertEqual(rows[0]['changed_boundaries'], 'True')

    def test_single_and_two_exons(self):
        for strand in ('+', '-'):
            result, rows = self.run_check(annotation([(100, 150)], strand=strand),
                                          annotation([(90, 160)], strand=strand))
            self.assertEqual(result.returncode, 0, result.stderr)
            result, rows = self.run_check(annotation(EXONS[:2], strand=strand),
                                          annotation([(90, 150), (200, 260)], strand=strand))
            self.assertEqual(result.returncode, 0, result.stderr)
            result, rows = self.run_check(annotation(EXONS[:2], strand=strand),
                                          annotation([(90, 151), (201, 260)], strand=strand))
            self.assertEqual(result.returncode, 1, result.stderr)
            self.assertEqual(rows[0]['changed_boundary_details'], 'exon_1_end:150->151;exon_2_start:200->201')

    def test_missing_ids(self):
        result, rows = self.run_check(annotation(identifier='t1.1') + annotation(identifier='shared'),
                                      annotation(identifier='t1.2') + annotation(identifier='shared'))
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual([row['status'] for row in rows],
                         ['unchanged_internal_structure', 'missing_from_extended', 'missing_from_raw'])
        self.assertEqual(rows[1]['exon_coordinates_extended'], '')
        self.assertIn('missing_from_extended: 1', result.stdout)
        self.assertIn('missing_from_raw: 1', result.stdout)

    def test_invalid_inputs(self):
        valid = annotation()
        lines = valid.splitlines(keepends=True)
        invalids = [('', 'no transcript'), ('# comment\n', 'no transcript'),
                    (valid + lines[0], 'duplicate transcript'),
                    (valid.replace('transcript_id=t1;', 'gene_id=g1;'), 'missing transcript_id'),
                    (valid.replace('transcript_id=t1;', 'broken;'), 'malformed'),
                    (valid.replace('\t100\t150\t', '\tx\t150\t'), 'integers'),
                    (annotation([(0, 150)]), 'invalid coordinates'),
                    (annotation([(151, 150)]), 'invalid coordinates'),
                    (annotation([(100, 150), (150, 250)]), 'overlapping'),
                    (annotation([(100, 150), (100, 150)]), 'duplicate exons'),
                    (annotation([(100, 1001)]), 'outside transcript'),
                    (lines[0], 'no exons'), (''.join(lines[1:]), 'orphan exon'),
                    (lines[0] + ''.join(lines[1:]).replace('1\ttest', '2\ttest'), 'inconsistent exon'),
                    (valid.replace('\t+\t', '\tx\t'), 'invalid chromosome or strand'),
                    ('1\ttest\texon\n', '9 GFF3 columns')]
        for invalid, message in invalids:
            for raw, new in [(invalid, valid), (valid, invalid)]:
                with self.subTest(message=message, raw=raw):
                    result, rows = self.run_check(raw, new)
                    self.assertEqual(result.returncode, 2, result.stderr)
                    self.assertIn(message, result.stderr)
                    self.assertIsNone(rows)

    def test_file_errors_and_output_protection(self):
        with tempfile.TemporaryDirectory() as directory:
            raw = Path(directory) / 'raw.gff3'
            raw.write_text(annotation())
            alias = Path(directory) / 'alias.gff3'
            alias.symlink_to(raw)
            for output in (raw, alias, Path(directory) / 'absent/report.csv'):
                result = subprocess.run([sys.executable, str(SCRIPT), str(raw), str(raw), str(output)],
                                        capture_output=True, text=True)
                self.assertEqual(result.returncode, 2)
                self.assertEqual(raw.read_text(), annotation())
            result = subprocess.run([sys.executable, str(SCRIPT), directory + '/missing',
                                     str(raw), directory + '/out.csv'], capture_output=True, text=True)
            self.assertEqual(result.returncode, 2)


if __name__ == '__main__':
    unittest.main()
