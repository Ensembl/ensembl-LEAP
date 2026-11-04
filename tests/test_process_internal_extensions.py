"""Tests against small, inspectable LEAP directory fixtures."""
import csv
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / 'assets/process_internal_extensions.py'
DATA = ROOT / 'tests/data/internal_extensions'

def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

process = load('process_internal', SCRIPT)
check = load('overlap_check', ROOT / 'assets/check_cds_overlaps.py')


@unittest.skipUnless(shutil.which('sort'), 'requires sort')
class InternalProcessingTests(unittest.TestCase):
    def test_all_fixture_strands_and_stats(self):
        for name, expected in json.loads((DATA / 'expected.json').read_text()).items():
            with self.subTest(name=name), tempfile.TemporaryDirectory() as directory:
                source, output = DATA / name, Path(directory) / 'processed'
                before = {str(p): p.read_bytes() for p in source.rglob('*') if p.is_file()}
                result = subprocess.run([sys.executable, str(SCRIPT), str(source), str(output)], capture_output=True, text=True)
                self.assertEqual(result.returncode, 0, result.stderr)
                folder = output / 'external_extensions'
                strand = expected['strand']
                path = folder / 'gff' / ('1_' + strand + '_extended.gff')
                raw = check.read_annotation(source / 'raw_annotation/current_annotation_1.gff3')
                extended = check.read_annotation(path)
                stats = check.read_extension_stats(folder / 'stats')
                check.validate_extension_stats(raw, extended, stats)
                self.assertEqual(set(extended), {'anchor.1', 'source.1'})
                source_model = extended['source.1']
                self.assertEqual((min(a for a,b in source_model['exon']), max(b for a,b in source_model['exon'])),
                                 (expected['start'], expected['end']))
                self.assertEqual(source_model['exon'][1], raw['source.1']['exon'][1])
                self.assertEqual(source_model['CDS'], raw['source.1']['CDS'])
                self.assertEqual(len(source_model['exon']), 3)
                row = stats[('1', 'source.1')]
                for side, restored in [('left', expected['restore_left']), ('right', expected['restore_right'])]:
                    prime = 'five_prime' if (side == 'left') == (strand == 'positive') else 'three_prime'
                    self.assertEqual(row[prime + '_extension'], 0 if restored else 50)
                    self.assertEqual(row[prime + '_support'], 'none' if restored else ('current_annotation' if prime == 'five_prime' else 'reads'))
                    self.assertEqual(row[prime + '_read_count'], '0' if restored else ('7' if prime == 'five_prime' else '3'))
                potential = folder / 'potential_internal' / ('1_' + strand + '_extended.gff')
                original = source / 'extended_transcripts' / path.name
                source_lines = [line for line in original.read_text().splitlines() if 'transcript_id=source.1;' in line]
                potential_lines = [line for line in potential.read_text().splitlines() if 'transcript_id=source.1;' in line]
                if expected['restore_left'] or expected['restore_right']:
                    self.assertCountEqual(source_lines, potential_lines)
                else:
                    self.assertEqual(potential_lines, [])
                self.assertTrue('source.1' in extended)  # Both-restored stats rows are retained.
                for path, content in before.items():
                    self.assertEqual(Path(path).read_bytes(), content)

    def test_classification_retains_boundary_or_rule(self):
        source = DATA / 'external_positive/extended_transcripts/1_positive_extended.gff'
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'test.gff'
            # Shares the earliest CDS start but not its end; LEAP2GO treats this as a match.
            text = source.read_text().replace('\t120\t150\t.\t+\t.\ttranscript_id=source.1;', '\t120\t160\t.\t+\t.\ttranscript_id=source.1;')
            path.write_text(text)
            models, dropped = process.scan_extended(path)
            self.assertEqual(models[('1','source.1')]['class'], 'external')
            self.assertIn(('1','utr_only.1'), dropped)

    def test_sort_failure_aborts_and_cleans_staging(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / 'processed'
            with patch.object(process.subprocess, 'run', side_effect=subprocess.CalledProcessError(1, 'sort')):
                with self.assertRaises(subprocess.CalledProcessError):
                    process.run(DATA / 'internal_left_negative', output)
            self.assertFalse(output.exists())
            self.assertFalse(list(Path(directory).glob('.processed-build-*')))

    def test_existing_output_and_stale_stats_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / 'existing'
            output.mkdir()
            with self.assertRaisesRegex(ValueError, 'already exists'):
                process.run(DATA / 'internal_left_positive', output)
            source = Path(directory) / 'input'
            shutil.copytree(DATA / 'internal_left_positive', source)
            stats = source / 'extended_transcripts_stats/1_positive_stats.tsv'
            stats.write_text(stats.read_text().replace('1:250-950', '1:251-950'))
            with self.assertRaisesRegex(ValueError, 'gene/location differs'):
                process.run(source, Path(directory) / 'bad')
            self.assertFalse((Path(directory) / 'bad').exists())
            extended = source / 'extended_transcripts/1_positive_extended.gff'
            shutil.copy2(extended, extended.with_suffix('.gff3'))
            with self.assertRaisesRegex(ValueError, 'Duplicate extended annotation basenames'):
                process.run(source, Path(directory) / 'duplicate')

    def test_single_exon_restores_whole_exon_and_both_stats_ends(self):
        with tempfile.TemporaryDirectory() as directory:
            source, output = Path(directory) / 'input', Path(directory) / 'processed'
            shutil.copytree(DATA / 'internal_both_positive', source)
            for relative, bounds in [('raw_annotation/current_annotation_1.gff3', (300, 400)),
                                     ('extended_transcripts/1_positive_extended.gff', (250, 450))]:
                path = source / relative
                rows = []
                for fields, attrs, context in process.records(path):
                    if attrs.get('transcript_id') == 'source.1':
                        if fields[2] == 'transcript':
                            fields[3:5] = list(map(str, bounds))
                        elif attrs.get('exon_number') != '1':
                            continue
                        else:
                            if int(fields[3]) in (250, 300):
                                fields[3] = str(bounds[0])
                            if int(fields[4]) == 400:
                                fields[4] = str(bounds[1])
                    rows.append(process.serialize(fields))
                path.write_text(''.join(rows))
            stats_path = source / 'extended_transcripts_stats/1_positive_stats.tsv'
            stats_path.write_text(stats_path.read_text().replace('1:250-850', '1:250-450'))
            process.run(source, output)
            folder = output / 'external_extensions'
            raw = check.read_annotation(source / 'raw_annotation/current_annotation_1.gff3')
            after = check.read_annotation(folder / 'gff/1_positive_extended.gff')
            stats = check.read_extension_stats(folder / 'stats')
            check.validate_extension_stats(raw, after, stats)
            self.assertEqual(after['source.1']['exon'], [(300, 400)])
            self.assertEqual(stats[('1', 'source.1')]['five_prime_extension'], 0)
            self.assertEqual(stats[('1', 'source.1')]['three_prime_extension'], 0)


if __name__ == '__main__':
    unittest.main()
