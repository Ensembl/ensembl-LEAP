import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

from datacheck_fixtures import write_run

ROOT = Path(__file__).resolve().parents[1]
NEXTFLOW = os.environ.get('LEAP_NEXTFLOW') or shutil.which('nextflow')


@unittest.skipUnless(NEXTFLOW, 'requires Nextflow')
class ManualReviewNextflowTests(unittest.TestCase):
    def run_workflow(self, scenario):
        temp = tempfile.TemporaryDirectory(); self.addCleanup(temp.cleanup)
        root = Path(temp.name)
        inputs = write_run(root / 'inputs', empty_negative=scenario == 'empty')
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
        summary = json.loads((output / 'review-summary.json').read_text())
        self.assertEqual(len(summary['scopes']), 2)
        self.assertEqual(summary['settings']['published_root'], str(output))

    def test_empty_partition_is_retained_and_missing_partition_fails(self):
        result, root = self.run_workflow('empty')
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual((root / 'published/manual_review/annotations/leaped/negative.gff3').read_text(), '##gff-version 3\n')
        result, root = self.run_workflow('missing')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('manual-review annotation partitions', result.stdout + result.stderr)

    def test_failed_datacheck_waits_for_review_bundle(self):
        result, root = self.run_workflow('failed_check')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('DATACHECK_GATE', result.stdout + result.stderr)
        self.assertTrue((root / 'published/manual_review/samplesheet.csv').is_file(), result.stdout + result.stderr)
        self.assertTrue((root / 'published/datachecks/summary/summary.json').is_file())

    def test_without_reference_lengths_still_writes_editable_bundle(self):
        result, root = self.run_workflow('no_reference')
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        output = root / 'published/manual_review'
        self.assertTrue((output / 'run_jbrowser.sbatch').is_file())
        self.assertFalse((output / 'reference.chrom.sizes').exists())


if __name__ == '__main__':
    unittest.main()
