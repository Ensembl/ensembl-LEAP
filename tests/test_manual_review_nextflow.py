import csv
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

from datacheck_fixtures import write_run
from read_evidence_fixtures import alignments, sam_records, write_bam

ROOT = Path(__file__).resolve().parents[1]
NEXTFLOW = os.environ.get('LEAP_NEXTFLOW') or shutil.which('nextflow')


@unittest.skipUnless(NEXTFLOW and shutil.which('samtools'), 'requires Nextflow and samtools')
class ManualReviewNextflowTests(unittest.TestCase):
    def run_workflow(self, scenario):
        temp = tempfile.TemporaryDirectory(); self.addCleanup(temp.cleanup)
        root = Path(temp.name)
        inputs = write_run(root / 'inputs', empty_negative=scenario == 'empty')
        (inputs/'CAGE.bed').write_text('chr1\t80\t90\tcage_positive\t2000\t+\t82\t83\t255,0,0\nchr1\t760\t765\tcage_negative\t20\t-\t761\t762\t255,0,0\n')
        (inputs/'polyA.bed').write_text('1\t360\t361\tpolyA_positive\t0.25\t+\t1.0\t2\t0.25\tIG\tAATAAA\n1\t480\t481\tpolyA_negative\t0.5\t-\t2.0\t3\t0.5\tIG\tAATAAA\n')
        rows = [r for r in alignments() if r.split('\t')[2] == '1' and
                not (scenario == 'empty' and int(r.split('\t')[1]) & 16)]
        write_bam(inputs / 'reads/1.bam', rows)
        if scenario == 'failed_check':
            (inputs / 'partial_protection/1_positive_partial_protection.tsv').write_text('malformed\n')
        config = root / 'local.config'
        config.write_text("process.executor='local'\nprocess.cpus=1\nprocess.memory='1 GB'\nexecutor.queueSize=6\nsingularity.enabled=false\nconda.enabled=false\n")
        command = [NEXTFLOW, '-C', str(config), 'run', str(ROOT / 'tests/nextflow/manual_review.nf'),
            '-ansi-log', 'false', '-work-dir', str(root / 'work'), '--fixture_dir', str(inputs),
            '--output_directory', str(root / 'published')]
        if scenario == 'missing': command += ['--missing_negative', 'true']
        if scenario == 'failed_check': command += ['--with_checks', 'true']
        if scenario == 'no_reference': command += ['--no_reference', 'true']
        if scenario == 'missing_reads': command += ['--missing_reads', 'true']
        env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', NXF_OFFLINE='true', NXF_DISABLE_CHECK_LATEST='true')
        result = subprocess.run(command, cwd=root, env=env, capture_output=True, text=True, timeout=180)
        return result, root

    def test_all_modules_publish_twelve_tracks_and_sbatch(self):
        result, root = self.run_workflow('pass')
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        output = root / 'published/manual_review'
        self.assertEqual(len(list((output / 'annotations').glob('*/*.gff3'))), 12)
        self.assertTrue((output / 'run_jbrowser.sbatch').is_file())
        self.assertTrue((output / 'run_jbrowser_offline.sbatch').is_file())
        self.assertTrue((output / 'samplesheet_offline.csv').is_file())
        self.assertTrue((output / 'spreadsheets/transcripts.csv').is_file())
        with (output / 'samplesheet.csv').open() as handle:
            spreadsheet = next(r for r in csv.DictReader(handle) if r['type']=='spreadsheet')
        self.assertEqual(spreadsheet['data_category'], 'leap_review')
        summary = json.loads((output / 'review-summary.json').read_text())
        self.assertEqual(len(summary['scopes']), 2)
        self.assertEqual(summary['settings']['published_root'], str(output))
        self.assertEqual(summary['settings']['species'], 'hg38')
        for sheet in ('samplesheet.csv', 'samplesheet_offline.csv'):
            with (output / sheet).open() as handle:
                rows = list(csv.DictReader(handle))
            self.assertTrue(all(r['species']=='hg38' and r['genome_id']=='hg38' for r in rows))
        self.assertEqual(summary['read_evidence']['counts'], dict(positive=2, negative=2, unsplit=4))
        self.assertTrue((output / 'reads/unsplit.bam.csi').is_file())
        self.assertEqual(len(list((output/'sites').glob('*/*.bed'))),6)
        with (output/'samplesheet.csv').open() as handle:
            sites = [r for r in csv.DictReader(handle) if r['type']=='regions']
        self.assertEqual(len(sites),6)
        self.assertEqual({r['region_format'] for r in sites},{'cageBed9','polyAtlas'})
        self.assertTrue((output/'region_presets/cageBed9.tsv').is_file())
        self.assertTrue((output/'install_review_region_presets.py').is_file())

    def test_empty_partition_is_retained_and_missing_partition_fails(self):
        result, root = self.run_workflow('empty')
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual((root / 'published/manual_review/annotations/leaped/negative.gff3').read_text(), '##gff-version 3\n')
        self.assertEqual(sam_records(root / 'published/manual_review/reads/negative.bam'), [])
        self.assertTrue((root / 'published/manual_review/reads/negative.bam.csi').is_file())
        with (root / 'published/manual_review/samplesheet.csv').open() as handle:
            self.assertNotIn('LEAP_reads_negative', {r['sample'] for r in csv.DictReader(handle)})
        result, root = self.run_workflow('missing')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('manual-review annotation partitions', result.stdout + result.stderr)

    def test_failed_datacheck_waits_for_review_bundle(self):
        result, root = self.run_workflow('failed_check')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('DATACHECK_GATE', result.stdout + result.stderr)
        self.assertTrue((root / 'published/manual_review/samplesheet.csv').is_file(), result.stdout + result.stderr)
        self.assertTrue((root / 'published/datachecks/summary/summary.json').is_file())
        self.assertTrue((root / 'published/manual_review/reads/unsplit.bam.csi').is_file())
        self.assertTrue((root / 'published/impact/index.html').is_file())
        self.assertTrue((root/'published/manual_review/sites/cage/unsplit.bed').is_file())
        with (root/'published/manual_review/spreadsheets/transcripts.csv').open() as handle:
            reader = csv.DictReader(handle); fields = reader.fieldnames; review_rows = list(reader)
        with (root/'published/impact/transcripts.csv').open() as handle:
            impact_rows = list(csv.DictReader(handle))
        self.assertNotIn('raw_start', fields)
        self.assertEqual(review_rows, [{k:row[k] for k in fields} for row in impact_rows])

    def test_without_reference_lengths_still_writes_editable_bundle(self):
        result, root = self.run_workflow('no_reference')
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        output = root / 'published/manual_review'
        self.assertTrue((output / 'run_jbrowser.sbatch').is_file())
        self.assertFalse((output / 'reference.chrom.sizes').exists())

    def test_missing_read_partition_does_not_silently_make_annotation_only_bundle(self):
        result, root = self.run_workflow('missing_reads')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('manual-review chromosome BAM/index pairs', result.stdout + result.stderr)
        self.assertFalse((root / 'published/manual_review').exists())


if __name__ == '__main__':
    unittest.main()
