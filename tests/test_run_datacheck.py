import csv
import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from datacheck_fixtures import write_run

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT/'bin/run_datacheck.py'
SPEC = importlib.util.spec_from_file_location('datacheck_runner',SCRIPT)
RUNNER = importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(RUNNER)


class DatacheckAdapterTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.context = write_run(self.root/'context')
        self.index = 0

    def run_check(self, name, strand='positive', checkers=None):
        self.index += 1
        output = self.root/('result-'+str(self.index))
        result = subprocess.run([sys.executable,str(SCRIPT),'check','--check',name,
            '--chromosome','1','--strand',strand,'--context',str(self.context),
            '--checkers',str(checkers or ROOT/'assets'),'--reference-lengths',str(self.context/'reference.sizes'),
            '--output',str(output)],capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stderr)
        status = json.loads((output/('1_'+strand+'_'+name+'.status.json')).read_text())
        return status,output

    def test_all_ten_checks_both_strands(self):
        for strand in ('positive','negative'):
            for name in RUNNER.CHECKS:
                with self.subTest(strand=strand,check=name):
                    status,output = self.run_check(name,strand)
                    self.assertEqual(status['status'],'PASS',status)
                    self.assertEqual(status['exit_code'],0)
                    self.assertTrue((output/'report.csv').is_file())

    def test_reference_inputs_are_not_in_published_reports(self):
        status,output = self.run_check('reference_coordinates')
        self.assertEqual(status['status'],'PASS',status)
        self.assertFalse((output/'reference_inputs').exists())
        self.assertFalse(any(p.suffix in ('.gff','.gff3') for p in output.rglob('*')))
        self.assertTrue(output.with_name(output.name+'_reference_inputs').is_dir())

    def test_failures_and_input_errors_are_retained(self):
        gff = self.context/'extended_transcripts/1_positive_extended.gff'
        gff.write_text(gff.read_text().replace('\tgene\t90\t360\t','\tgene\t91\t359\t'))
        status,output = self.run_check('structure')
        self.assertEqual(status['status'],'FAIL',status)
        self.assertEqual(status['exit_code'],1)
        self.assertIn('FAIL:',(output/'stdout.txt').read_text())
        (self.context/'partial_protection/1_positive_partial_protection.tsv').unlink()
        status,output = self.run_check('partial_protection')
        self.assertEqual(status['status'],'INPUT_ERROR',status)
        self.assertIn('ERROR:',(output/'stderr.txt').read_text())
        self.assertIn('INPUT_ERROR',(output/'report.csv').read_text())

    def test_empty_partitions_require_verified_completeness(self):
        write_run(self.context,empty_negative=True)
        for name in RUNNER.CHECKS:
            status,output = self.run_check(name,'negative')
            self.assertEqual(status['status'],'NOT_APPLICABLE' if name in RUNNER.LEGACY_EMPTY else 'PASS',status)
            self.assertEqual(status['exit_code'],0)
        # Eligible missing models cannot become NOT_APPLICABLE.
        write_run(self.context)
        (self.context/'extended_transcripts/1_negative_extended.gff').write_text('')
        status,_ = self.run_check('preservation','negative')
        self.assertEqual(status['status'],'INPUT_ERROR',status)
        status,_ = self.run_check('completeness','negative')
        self.assertEqual(status['status'],'FAIL',status)

    def test_unexpected_checker_errors_and_missing_report(self):
        checkers = self.root/'checkers';checkers.mkdir()
        path = checkers/RUNNER.CHECKS['boundary']
        for source in ('raise RuntimeError("unexpected")\n','pass\n','import sys;sys.exit(137)\n'):
            path.write_text(source)
            status,output = self.run_check('boundary',checkers=checkers)
            self.assertEqual(status['status'],'INPUT_ERROR',status)
            self.assertEqual(status['exit_code'],2)
            self.assertTrue((output/'report.csv').exists())

    def summary(self, statuses, scopes=None, issues=None):
        folder = self.root/('statuses-'+str(self.index));folder.mkdir()
        for i,status in enumerate(statuses):
            (folder/(str(i)+'.status.json')).write_text(json.dumps(status))
        manifest = self.root/('expected-'+str(self.index)+'.json')
        manifest.write_text(json.dumps(dict(scopes=scopes if scopes is not None else [['1','positive']],issues=issues or [])))
        output = self.root/('summary-'+str(self.index));self.index+=1
        result = subprocess.run([sys.executable,str(SCRIPT),'summarize','--expected',str(manifest),
            '--statuses',str(folder),'--output',str(output)],capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stderr)
        return json.loads((output/'summary.json').read_text()),output

    def valid_statuses(self):
        return [dict(chromosome='1',strand='positive',check=name,status='PASS',exit_code=0,
            report='1_positive_'+name+'/report.csv',stdout='1_positive_'+name+'/stdout.txt',
            stderr='1_positive_'+name+'/stderr.txt',message='Pass',elapsed_seconds=1) for name in RUNNER.CHECKS]

    def test_summary_and_final_gate_only_fail_after_reports_exist(self):
        statuses = self.valid_statuses()
        statuses[0].update(status='FAIL',exit_code=1,message='<script>unsafe</script>')
        summary,output = self.summary(statuses)
        self.assertEqual(summary['exit_code'],1)
        self.assertEqual(summary['counts'],{'FAIL':1,'PASS':9})
        self.assertIn('&lt;script&gt;', (output/'summary.html').read_text())
        self.assertNotIn('<script>', (output/'summary.html').read_text())
        result = subprocess.run([sys.executable,str(SCRIPT),'gate','--verdict',str(output/'verdict.json')],cwd=self.root,capture_output=True,text=True)
        self.assertEqual(result.returncode,1)
        self.assertTrue((output/'summary.csv').is_file())
        statuses[0].update(status='INPUT_ERROR',exit_code=2)
        summary,output = self.summary(statuses)
        self.assertEqual(summary['exit_code'],2)
        result = subprocess.run([sys.executable,str(SCRIPT),'gate','--verdict',str(output/'verdict.json')],cwd=self.root,capture_output=True,text=True)
        self.assertEqual(result.returncode,2)
        summary,output = self.summary(self.valid_statuses())
        result = subprocess.run([sys.executable,str(SCRIPT),'gate','--verdict',str(output/'verdict.json')],cwd=self.root,capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stderr)
        self.assertTrue((self.root/'datachecks.passed').exists())

    def test_missing_duplicate_extra_results_and_empty_scope(self):
        statuses = self.valid_statuses()
        for observed in (statuses[:-1],statuses+[statuses[0]],statuses+[dict(statuses[0],chromosome='2')]):
            summary,_ = self.summary(observed)
            self.assertEqual(summary['exit_code'],2)
            self.assertGreater(summary['counts']['INPUT_ERROR'],0)
        summary,_ = self.summary([],scopes=[])
        self.assertEqual(summary['exit_code'],2)
        summary,_ = self.summary(statuses,issues=['Unexpected final output chromosome 2'])
        self.assertEqual(summary['exit_code'],2)

    def test_duplicate_findings_are_review_and_do_not_fail_gate(self):
        path=self.context/'extended_transcripts/1_positive_extended.gff'
        original=path.read_text()
        copied=''.join(line.replace('t_plus.1','t_copy.1') for line in original.splitlines(True) if '\tgene\t' not in line)
        path.write_text(original+copied)
        status, output=self.run_check('duplicate_transcripts')
        self.assertEqual(status['status'],'REVIEW',status)
        self.assertEqual(status['exit_code'],0)
        self.assertEqual(status['checker_exit_code'],1)
        statuses=self.valid_statuses()
        statuses=[status if row['check']=='duplicate_transcripts' else row for row in statuses]
        summary,folder=self.summary(statuses)
        self.assertEqual(summary['status'],'REVIEW')
        self.assertEqual(summary['exit_code'],0)
        result=subprocess.run([sys.executable,str(SCRIPT),'gate','--verdict',str(folder/'verdict.json')],cwd=self.root,capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stderr)
        self.assertIn('require review',result.stdout)
        # REVIEW cannot disguise a failure in another datacheck.
        statuses[0].update(status='REVIEW',exit_code=0,checker_exit_code=1)
        summary,_=self.summary(statuses)
        self.assertEqual(summary['exit_code'],2)

    def test_duplicate_input_error_and_missing_report_are_errors(self):
        path=self.context/'extended_transcripts/1_positive_extended.gff'
        path.write_text('malformed GFF\n')
        status,_=self.run_check('duplicate_transcripts')
        self.assertEqual(status['status'],'INPUT_ERROR')
        checkers=self.root/'duplicate-checkers';checkers.mkdir()
        (checkers/RUNNER.CHECKS['duplicate_transcripts']).write_text('import sys;sys.exit(1)\n')
        status,_=self.run_check('duplicate_transcripts',checkers=checkers)
        self.assertEqual(status['status'],'INPUT_ERROR')
