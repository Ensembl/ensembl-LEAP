from collections import Counter
import csv
import gzip
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from read_evidence_fixtures import alignments, sam_records, write_bam
from test_prepare_manual_review import fixture

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('review_reads', ROOT / 'bin/prepare_review_reads.py')
reads = importlib.util.module_from_spec(spec); spec.loader.exec_module(reads)
SCOPES = [('1', 'positive'), ('1', 'negative'), ('2', 'positive'), ('2', 'negative')]


@unittest.skipUnless(shutil.which('samtools'), 'requires samtools')
class ReviewReadTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(); self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def partitions(self, rows=None):
        bam, index = write_bam(self.root / 'source.bam', rows)
        # Ensure the explicit index works even without the usual sidecar name.
        renamed = index.with_name('custom source index.csi'); index.rename(renamed)
        parts = reads.partition_source(bam, renamed, ['1', '2'], self.root / 'parts')
        return bam, renamed, parts

    def test_record_multiplicity_flags_tags_and_strands_are_preserved(self):
        bam, index, parts = self.partitions()
        before = {p: p.read_bytes() for p in (bam, index)}
        output = self.root / 'review'
        reads.prepare_reads(parts, SCOPES, output)
        expected = [r for r in sam_records(bam, '-F', '4') if r.split('\t')[2] in ('1', '2')]
        def normalized(rows):
            # Samtools may move RG among optional tags; tag values stay intact.
            return Counter(tuple(r.split('\t')[:11]) + tuple(sorted(r.split('\t')[11:])) for r in rows)
        self.assertEqual(normalized(sam_records(output / 'unsplit.bam')), normalized(expected))
        positive = sam_records(output / 'positive.bam'); negative = sam_records(output / 'negative.bam')
        self.assertTrue(all(not int(r.split('\t')[1]) & 16 for r in positive))
        self.assertTrue(all(int(r.split('\t')[1]) & 16 for r in negative))
        self.assertEqual(normalized(positive + negative), normalized(expected))
        summary = json.loads((output / 'read-summary.json').read_text())
        self.assertEqual(summary['counts'], dict(positive=3, negative=3, unsplit=6))
        self.assertEqual({p: p.read_bytes() for p in (bam, index)}, before)
        self.assertIn('@RG\tID:sample\tSM:fixture', reads.command('view', '-H', output / 'unsplit.bam'))
        for view in reads.VIEWS:
            with gzip.open(output / (view + '.bam.csi'), 'rb') as handle:
                self.assertEqual(handle.read(4), b'CSI\x01')

    def test_empty_views_and_restricted_strands(self):
        _, _, parts = self.partitions([r for r in alignments() if not int(r.split('\t')[1]) & 16])
        reads.prepare_reads(parts, SCOPES, self.root / 'positive_only')
        self.assertEqual(sam_records(self.root / 'positive_only/negative.bam'), [])
        reads.prepare_reads(parts, [('1', 'negative'), ('2', 'negative')], self.root / 'empty')
        summary = json.loads((self.root / 'empty/read-summary.json').read_text())
        self.assertEqual(summary['counts'], dict(positive=0, negative=0, unsplit=0))
        self.assertTrue(all((self.root / 'empty' / (v + '.bam.csi')).is_file() for v in reads.VIEWS))

    def test_identical_alignment_records_are_not_collapsed(self):
        row = alignments()[0]
        bam, index = write_bam(self.root / 'duplicates.bam', [row, row])
        reads.prepare_reads([dict(chromosome='1', bam=bam, index=index)],
            [('1', 'positive'), ('1', 'negative')], self.root / 'duplicates')
        self.assertEqual(len(sam_records(self.root / 'duplicates/unsplit.bam')), 2)

    def test_csi_supports_coordinates_beyond_bai_limits(self):
        row = alignments()[0].split('\t'); row[3] = '600000001'
        bam, index = write_bam(self.root / 'long.bam', ['\t'.join(row)], {'1': 700000000})
        reads.prepare_reads([dict(chromosome='1', bam=bam, index=index)],
            [('1', 'positive'), ('1', 'negative')], self.root / 'long_review')
        count = reads.command('view', '-c', '-X', self.root / 'long_review/positive.bam',
            self.root / 'long_review/positive.bam.csi', '1:600000001-600000010')
        self.assertEqual(int(count), 1)

    def test_scope_missing_duplicate_headers_indexes_and_overwrite_fail(self):
        bam, index, parts = self.partitions()
        for supplied in (parts[:1], parts + parts[:1]):
            with self.assertRaisesRegex(ValueError, 'Missing, duplicate or unexpected'):
                reads.prepare_reads(supplied, SCOPES, self.root / 'bad')
        with self.assertRaisesRegex(ValueError, 'outside chromosome'):
            reads.prepare_reads([dict(chromosome='1', bam=bam, index=index)], [('1', 'positive')], self.root / 'bad')
        for kwargs in (dict(group='different'), dict(lengths={'1': 3000, '2': 2000, 'MT': 2000})):
            other, other_index = write_bam(self.root / 'other.bam', [alignments()[4]], **kwargs)
            with self.assertRaisesRegex(ValueError, 'Incompatible BAM'):
                reads.prepare_reads([parts[0], dict(chromosome='2', bam=other, index=other_index)], SCOPES, self.root / 'bad')
        bad_index = self.root / 'bad.csi'; bad_index.write_bytes(b'invalid')
        with self.assertRaisesRegex(ValueError, 'Samtools view failed'):
            reads.prepare_reads([dict(parts[0], index=bad_index), parts[1]], SCOPES, self.root / 'bad')
        with self.assertRaisesRegex(ValueError, 'Missing or unreadable'):
            reads.prepare_reads([dict(parts[0], index=self.root / 'missing'), parts[1]], SCOPES, self.root / 'bad')
        invalid = self.root / 'invalid.bam'; invalid.write_bytes(b'not a BAM')
        with self.assertRaisesRegex(ValueError, 'Samtools quickcheck failed'):
            reads.prepare_reads([dict(parts[0], bam=invalid), parts[1]], SCOPES, self.root / 'bad')
        with self.assertRaisesRegex(ValueError, 'new directory'):
            reads.prepare_reads(parts, SCOPES, self.root)
        self.assertFalse((self.root / 'bad').exists())

    def test_failed_sort_removes_staging_and_does_not_publish(self):
        _, _, parts = self.partitions()
        run = reads.subprocess.run
        def failure(command, *args, **kwargs):
            if command[:2] == ['samtools', 'sort']:
                return subprocess.CompletedProcess(command, 1, stderr='forced sorting failure')
            return run(command, *args, **kwargs)
        with patch.object(reads.subprocess, 'run', side_effect=failure):
            with self.assertRaisesRegex(ValueError, 'filtering/sorting failed'):
                reads.prepare_reads(parts, SCOPES, self.root / 'bad')
        self.assertFalse((self.root / 'bad').exists())
        self.assertFalse(list(self.root.glob('.bad-*')))

    def test_standalone_sheets_publish_alignments_and_keep_annotation_only_compatibility(self):
        source = fixture(self.root / 'leap')
        bam, index = write_bam(self.root / "evidence's reads.bam")
        output = self.root / 'bundle'
        reference = self.root / 'reference.fa'
        reference.write_text(''.join('>{}\n{}\n'.format(c, 'A' * 2000) for c in ('1', '2', 'MT')))
        pipeline = self.root / 'jbrowser'; pipeline.mkdir(); (pipeline / 'main.nf').touch()
        result = subprocess.run([sys.executable, str(ROOT / 'bin/prepare_manual_review.py'), 'run', str(source), str(output),
            '--chromosomes', '1', '--reads', str(bam), '--reads-index', str(index),
            '--reference-genome', str(reference), '--genome-id', 'test', '--species', 'test',
            '--jbrowser-pipeline', str(pipeline), '--public-base-url', 'https://example.org/root',
            '--ftp-target-dir', '/example/ftp/root', '--hub-email', 'test@example.org'], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        with (output / 'samplesheet.csv').open() as handle:
            rows = [r for r in csv.DictReader(handle) if r['type'] == 'bam']
        self.assertEqual(len(rows), 3)
        self.assertTrue(all(r['reference_file'] == str(reference.resolve()) and r['data_category'] == 'alignment' for r in rows))
        self.assertTrue(all(Path(r['file']).is_file() for r in rows))
        with (output / 'samplesheet_offline.csv').open() as handle:
            self.assertTrue(all(r['publication'] == 'offline' for r in csv.DictReader(handle)))
        summary = json.loads((output / 'review-summary.json').read_text())
        self.assertEqual(summary['read_evidence']['counts']['unsplit'], 4)
        self.assertEqual(len(summary['tracks']), 15)
        self.assertEqual(summary['read_evidence']['original_input']['bam'], str(bam.resolve()))
        bins = self.root / 'executables'; bins.mkdir()
        captured = self.root / 'jbrowser-arguments.json'
        fake = bins / 'nextflow'
        fake.write_text('#!' + sys.executable + '\nimport json,sys\nfrom pathlib import Path\nPath(' +
            repr(str(captured)) + ').write_text(json.dumps(sys.argv[1:]))\n')
        fake.chmod(0o755)
        env = dict(os.environ, PATH=str(bins) + os.pathsep + os.environ['PATH'])
        for mode in ('', '_offline'):
            result = subprocess.run(['bash', str(output / ('run_jbrowser' + mode + '.sbatch'))],
                                    env=env, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            arguments = json.loads(captured.read_text())
            self.assertEqual(arguments[arguments.index('--input') + 1], str(output / ('samplesheet' + mode + '.csv')))
            if mode:
                self.assertNotIn('--public_base_url', arguments)
        with (output / 'samplesheet_offline.csv').open() as handle:
            self.assertTrue(all(r['publication'] == 'offline' for r in csv.DictReader(handle)))
        with (output / 'samplesheet.csv').open() as handle:
            fields = list(csv.DictReader(handle).fieldnames)
        with (output / 'samplesheet.csv').open() as handle:
            rows = list(csv.DictReader(handle))
        for row in rows:
            if row['type'] == 'bam': row['reference_file'] = str(self.root / 'missing-reference.fa')
        with (output / 'samplesheet.csv').open('w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=fields); writer.writeheader(); writer.writerows(rows)
        captured.unlink()
        result = subprocess.run(['bash', str(output / 'run_jbrowser_offline.sbatch')], env=env, capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('Alignment reference is missing', result.stderr)
        self.assertFalse(captured.exists())
        result = subprocess.run([sys.executable, str(ROOT / 'bin/prepare_manual_review.py'), 'run', str(source),
            str(self.root / 'invalid'), '--chromosomes', '1', '--reads', str(bam)], capture_output=True, text=True)
        self.assertEqual(result.returncode, 2)
        self.assertIn('supplied together', result.stderr)

    @unittest.skipUnless(os.environ.get('LEAP_JBROWSER_PIPELINE'), 'set LEAP_JBROWSER_PIPELINE for JBrowser CRAM validation')
    def test_jbrowser_alignment_route_builds_cram_and_crai(self):
        _, _, parts = self.partitions()
        evidence = self.root / 'evidence'; reads.prepare_reads(parts, SCOPES, evidence)
        reference = self.root / 'reference.fa'
        reference.write_text(''.join('>{}\n{}\n'.format(c, 'A' * 2000) for c in ('1', '2', 'MT')))
        helper = Path(os.environ['LEAP_JBROWSER_PIPELINE']) / 'bin/alignment_cram_lib.py'
        for view in reads.VIEWS:
            cram = self.root / (view + '.cram')
            result = subprocess.run([sys.executable, str(helper), 'validate', '--input', str(evidence / (view + '.bam')),
                '--sample', view, '--reference-file', str(reference), '--validation-tsv', str(self.root / (view + '.tsv')),
                '--mqc-json', str(self.root / (view + '-validation.json')), '--status-file', str(self.root / (view + '.status'))],
                capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual((self.root / (view + '.status')).read_text().strip(), 'PASS')
            result = subprocess.run([sys.executable, str(helper), 'prepare', '--input', str(evidence / (view + '.bam')),
                '--sample', view, '--reference-file', str(reference), '--output-cram', str(cram),
                '--output-crai', str(cram) + '.crai', '--mqc-json', str(self.root / (view + '.json'))], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertTrue(Path(str(cram) + '.crai').is_file())
            self.assertEqual(int(reads.command('view', '-c', '-T', reference, cram)), 6 if view == 'unsplit' else 3)


if __name__ == '__main__':
    unittest.main()
