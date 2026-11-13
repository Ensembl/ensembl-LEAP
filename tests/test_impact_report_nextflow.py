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
class ImpactNextflowTests(unittest.TestCase):
    def run_workflow(self, invalid=False):
        temp = tempfile.TemporaryDirectory(); self.addCleanup(temp.cleanup)
        root = Path(temp.name)
        run = write_run(root / 'inputs')
        if invalid:
            (run / 'extended_transcripts_stats/1_positive_stats.tsv').write_text('malformed\n')
        config = root / 'local.config'
        config.write_text("process.executor='local'\nprocess.cpus=1\nprocess.memory='1 GB'\nsingularity.enabled=false\nconda.enabled=false\n")
        result = subprocess.run([NEXTFLOW, '-C', str(config), 'run', str(ROOT / 'tests/nextflow/impact_report.nf'),
            '-ansi-log', 'false', '-work-dir', str(root / 'work'), '--fixture_dir', str(run),
            '--output_directory', str(root / 'published')], cwd=root,
            env=dict(os.environ, PYTHONDONTWRITEBYTECODE='1', NXF_OFFLINE='true', NXF_DISABLE_CHECK_LATEST='true'),
            capture_output=True, text=True, timeout=120)
        self.assertTrue((root / 'published/impact/index.html').is_file(), result.stdout + result.stderr)
        return result, json.loads((root / 'published/impact/summary.json').read_text())

    def test_stats_run_without_manual_review_or_datachecks(self):
        result, summary = self.run_workflow()
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(summary['datachecks']['status'], 'NOT_RUN')
        self.assertEqual(summary['totals']['external']['added_bases'], 40)

    def test_bad_stats_publish_dashboard_before_failure_even_with_checks_disabled(self):
        result, summary = self.run_workflow(invalid=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(summary['status'], 'INPUT_ERROR')
        self.assertIn('LEAP_IMPACT_GATE', result.stdout + result.stderr)


if __name__ == '__main__': unittest.main()
