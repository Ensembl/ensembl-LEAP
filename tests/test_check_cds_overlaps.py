"""CLI regression tests using the inspectable GFF3 datasets."""

import csv
import gzip
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / 'assets/check_cds_overlaps.py'
DATA = ROOT / 'tests/data/cds_overlaps'
spec = importlib.util.spec_from_file_location('cds_check', SCRIPT)
check = importlib.util.module_from_spec(spec)
spec.loader.exec_module(check)


from leap_test_context import write_context


class CdsOverlapTests(unittest.TestCase):
    def run_check(self, raw, extended, compressed=False, stats=None):
        with tempfile.TemporaryDirectory() as directory:
            paths = []
            for name, content in [('raw', raw), ('extended', extended)]:
                path = Path(directory) / (
                    name + ('.gff3.gz' if compressed else '.gff3')
                )
                with (gzip.open if compressed else open)(path, 'wt') as handle:
                    handle.write(content)
                paths.append(str(path))
            output = Path(directory) / 'report.csv'
            context = write_context(directory, *paths, stats=stats)
            options = ['--leap-output-dir', str(context)]
            if stats is not None:
                stats_path = Path(directory) / 'stats.tsv'
                stats_path.write_text(stats)
                options += ['--extension-stats', str(stats_path)]
            result = subprocess.run(
                [sys.executable, str(SCRIPT)] + paths + [str(output)] + options,
                capture_output=True,
                text=True,
            )
            rows = None
            if output.exists():
                with output.open(newline='') as handle:
                    rows = list(csv.DictReader(handle))
            return result, rows

    def dataset(self, name):
        return tuple(
            (DATA / name / filename).read_text()
            for filename in ('raw.gff3', 'extended.gff3')
        )

    def test_all_datasets(self):
        expected = json.loads((DATA / 'expected.json').read_text())
        for name, expectation in expected.items():
            with self.subTest(name=name):
                stats_path = DATA / name / 'stats.tsv'
                result, rows = self.run_check(
                    *self.dataset(name),
                    stats=stats_path.read_text() if stats_path.exists() else None
                )
                self.assertEqual(
                    result.returncode, expectation['exit_code'], result.stderr
                )
                if expectation['overlap_status']:
                    self.assertEqual(
                        {row['status'] for row in rows}, {expectation['overlap_status']}
                    )
                else:
                    self.assertEqual(rows, [])
                self.assertIn(
                    'FAIL:' if expectation['exit_code'] else 'PASS:', result.stdout
                )

    def test_report_evidence_and_mixed_outcomes(self):
        result, rows = self.run_check(*self.dataset('same_gene_mixed'), compressed=True)
        self.assertEqual(result.returncode, 1)
        self.assertEqual(len(rows), 1)
        row = rows[0]
        self.assertEqual(row['allowed_intervals'], '190-199')
        self.assertEqual(row['allowed_bases'], '10')
        self.assertEqual(row['failing_intervals'], '160-170')
        self.assertEqual(row['failing_bases'], '11')
        self.assertEqual(row['added_overlap_bases'], '21')
        self.assertEqual(row['shared_splice_sites'], 'end:250')
        self.assertIn('new_cds_overlap: 1', result.stdout)

    def test_annotation_support_is_per_end_and_reports_both_reasons(self):
        for name, status in [
            ('annotation_and_reads_opposite_ends', 'new_cds_overlap'),
            ('annotation_and_shared_site_allowed', 'allowed_multiple_reasons'),
        ]:
            with self.subTest(name=name):
                stats = (DATA / name / 'stats.tsv').read_text()
                result, rows = self.run_check(*self.dataset(name), stats=stats)
                self.assertEqual(rows[0]['status'], status)
                self.assertEqual(rows[0]['annotation_allowed_intervals'], '160-180')
                self.assertEqual(rows[0]['annotation_extension_ends'], 'five_prime')
                if status == 'new_cds_overlap':
                    self.assertEqual(result.returncode, 1)
                    self.assertEqual(rows[0]['failing_intervals'], '460-470')
                else:
                    self.assertEqual(result.returncode, 0)
                    self.assertEqual(
                        rows[0]['splice_site_allowed_intervals'], '460-470'
                    )

    def test_missing_or_stale_provenance_does_not_grant_exemption(self):
        name = 'annotation_same_gene_allowed'
        raw, extended = self.dataset(name)
        stats = (DATA / name / 'stats.tsv').read_text()
        result, rows = self.run_check(raw, extended)
        self.assertEqual(result.returncode, 1)
        self.assertIn('disabled', result.stdout)
        result, rows = self.run_check(
            raw, extended, stats=stats.replace('\tsource\t', '\tother\t')
        )
        self.assertEqual(result.returncode, 1)
        invalids = [
            (stats.replace('\t50\t', '\t51\t'), 'length differs'),
            (stats.replace('geneA', 'wrong_gene'), 'gene/location differs'),
            (stats.replace('1:150-450', '1:151-450'), 'gene/location differs'),
            (stats + stats.splitlines()[1] + '\n', 'duplicate stats'),
            (
                stats.replace('current_annotation', 'unknown'),
                'unknown extension support',
            ),
            ('wrong\theader\n', 'missing extension stats columns'),
        ]
        for invalid, message in invalids:
            with self.subTest(message=message):
                result, rows = self.run_check(raw, extended, stats=invalid)
                self.assertEqual(result.returncode, 2, result.stderr)
                self.assertIn(message, result.stderr)
                self.assertIsNone(rows)

    def test_stats_directory_and_gzip(self):
        name = 'annotation_cross_gene_cds_allowed'
        with tempfile.TemporaryDirectory() as directory:
            stats_dir = Path(directory)
            with gzip.open(stats_dir / 'extension.tsv.gz', 'wt') as handle:
                original = (DATA / name / 'stats.tsv').read_text()
                # The same transcript ID on a different chromosome is legitimate.
                other = original.splitlines()[1].replace('1:150-450', 'Y:150-450')
                handle.write(original + other + '\n')
            stats = check.read_extension_stats(stats_dir)
            self.assertEqual(
                stats[('1', 'source')]['five_prime_support'], 'current_annotation'
            )
            self.assertIn(('Y', 'source'), stats)
            stats_dir.joinpath('extension.tsv').write_text(
                (DATA / name / 'stats.tsv').read_text()
            )
            with self.assertRaisesRegex(ValueError, 'duplicate stats'):
                check.read_extension_stats(stats_dir)

    def test_duplicate_intervals_and_unordered_records(self):
        raw, extended = self.dataset('same_gene_shared_site')
        baseline, baseline_rows = self.run_check(raw, extended)
        raw = '\n'.join(reversed(raw.splitlines())) + '\n' + raw
        # Duplicate interval features are fine; transcript records are not duplicated.
        raw = (
            '\n'.join(
                line
                for i, line in enumerate(raw.splitlines())
                if '\ttranscript\t' not in line or i < len(raw.splitlines()) // 2
            )
            + '\n'
        )
        result, rows = self.run_check(
            raw, '\n'.join(reversed(extended.splitlines())) + '\n'
        )
        self.assertEqual(result.returncode, baseline.returncode, result.stderr)
        self.assertEqual(rows, baseline_rows)

    def test_missing_source_and_missing_feature_types(self):
        raw, extended = self.dataset('cross_gene_cds')
        result, rows = self.run_check(
            raw, extended.replace('transcript_id=source;', 'transcript_id=unknown;')
        )
        self.assertEqual(result.returncode, 0)
        self.assertEqual(rows[0]['status'], 'missing_from_raw')
        self.assertIn('missing_from_raw: 1', result.stdout)
        result, rows = self.run_check(
            raw,
            '\n'.join(line for line in extended.splitlines() if '_UTR\t' not in line)
            + '\n',
        )
        self.assertEqual(result.returncode, 0)
        self.assertEqual(rows, [])
        # Unrelated transcript-only models, such as the old OTT artifact, need no exons.
        raw += '1\ttest\ttranscript\t900\t950\t.\t+\t.\ttranscript_id=OTT;gene_id=other;transcript_type=protein_coding;\n'
        result, rows = self.run_check(raw, extended)
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_interval_helpers_and_index(self):
        self.assertEqual(
            list(check.subtract([(1, 20)], [(3, 5), (10, 15)])),
            [(1, 2), (6, 9), (16, 20)],
        )
        intervals = [(1, 100, 'long'), (20, 21, 'short'), (101, 110, 'next')]
        index = check.IntervalIndex(intervals)
        self.assertEqual(list(index.query(50, 50)), [intervals[0]])
        self.assertEqual(list(index.query(111, 120)), [])

    def test_invalid_inputs(self):
        raw, extended = self.dataset('same_gene_shared_site')
        invalids = [
            ('', 'no transcript-associated'),
            ('1\ttest\tCDS\n', '9 GFF3 columns'),
            (raw.replace('gene_id=geneA;', ''), 'missing gene_id'),
            (raw.replace('transcript_id=source;', 'broken;'), 'malformed'),
            (raw.replace('\t200\t250\t', '\t0\t250\t'), 'invalid coordinates'),
            (raw.replace('\t200\t250\t', '\tx\t250\t'), 'integers'),
            (raw.replace('\t+\t', '\t.\t'), 'explicit + or -'),
            (
                '\n'.join(line for line in raw.splitlines() if '\texon\t' not in line)
                + '\n',
                'containing exon',
            ),
        ]
        for invalid, message in invalids:
            with self.subTest(message=message):
                result, rows = self.run_check(invalid, extended)
                self.assertEqual(result.returncode, 2)
                self.assertIn(message, result.stderr)
                self.assertIsNone(rows)

    def test_output_protection_and_io_errors(self):
        with tempfile.TemporaryDirectory() as directory:
            raw = Path(directory) / 'raw.gff3'
            content = self.dataset('same_gene_shared_site')[0]
            raw.write_text(content)
            for source, output in [
                (raw, raw),
                (raw, Path(directory) / 'absent/report.csv'),
                (Path(directory) / 'missing.gff3', Path(directory) / 'out.csv'),
            ]:
                result = subprocess.run(
                    [sys.executable, str(SCRIPT), str(source), str(raw), str(output)],
                    capture_output=True,
                    text=True,
                )
                self.assertEqual(result.returncode, 2)
                self.assertEqual(raw.read_text(), content)


if __name__ == '__main__':
    unittest.main()
