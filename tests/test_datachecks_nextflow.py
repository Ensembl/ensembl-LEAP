import json
import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

from datacheck_fixtures import write_run

ROOT = Path(__file__).resolve().parents[1]
NEXTFLOW = os.environ.get('LEAP_NEXTFLOW') or shutil.which('nextflow')


@unittest.skipUnless(NEXTFLOW, 'Set LEAP_NEXTFLOW or install Nextflow for the real DSL2 smoke tests')
class NextflowDatachecksTests(unittest.TestCase):
    def run_workflow(self, scenario):
        temp = tempfile.TemporaryDirectory();self.addCleanup(temp.cleanup)
        root = Path(temp.name)
        inputs = write_run(root/'inputs',empty_negative=scenario=='empty')
        if scenario == 'fail':
            gff = inputs/'extended_transcripts/1_positive_extended.gff'
            gff.write_text(gff.read_text().replace('\tgene\t90\t360\t','\tgene\t91\t359\t'))
            (inputs/'partial_protection/1_positive_partial_protection.tsv').write_text('malformed header\n')
        if scenario == 'read_overlap':
            raw = inputs/'raw_annotation/current_annotation_1.gff3'
            text = raw.read_text().replace('\tgene\t100\t350\t', '\tgene\t90\t350\t')
            attrs = 'gene_id=g_plus;transcript_id=read_blocker.1;transcript_type=protein_coding'
            for kind,start,end,identifier,parent in (
                    ('transcript',90,98,'read_blocker.1','g_plus'),
                    ('exon',90,98,'exon:read_blocker.1','read_blocker.1'),
                    ('CDS',92,96,'CDS:read_blocker.1','read_blocker.1')):
                text += '1\ttest\t{}\t{}\t{}\t.\t+\t{}\tID={};Parent={};{}\n'.format(
                    kind,start,end,'0' if kind=='CDS' else '.',identifier,parent,attrs)
            raw.write_text(text)
        if scenario == 'duplicates':
            import csv
            for path in (inputs/'raw_annotation/current_annotation_1.gff3',inputs/'extended_transcripts/1_positive_extended.gff'):
                original=path.read_text()
                copied=''.join(line.replace('t_plus.1','t_copy.1') for line in original.splitlines(True)
                    if '\tgene\t' not in line and 'transcript_id=t_plus.1;' in line)
                path.write_text(original+copied)
            path=inputs/'extended_transcripts_stats/1_positive_stats.tsv'
            with path.open() as handle:
                reader=csv.DictReader(handle,delimiter='\t');fields=reader.fieldnames;rows=list(reader)
            rows.append(dict(rows[0],transcript_id='t_copy.1'))
            with path.open('w',newline='') as handle:
                writer=csv.DictWriter(handle,fieldnames=fields,delimiter='\t');writer.writeheader();writer.writerows(rows)
        config = root/'local.config'
        config.write_text("process.executor='local'\nprocess.cpus=1\nprocess.memory='1 GB'\nexecutor.queueSize=6\n"+
                          "singularity.enabled=false\nconda.enabled=false\n")
        command = [NEXTFLOW,'-C',str(config),'run',str(ROOT/'tests/nextflow/datachecks.nf'),
                   '-ansi-log','false','-work-dir',str(root/'work'),
                   '--fixture_dir',str(inputs),'--output_directory',str(root/'published')]
        if scenario == 'pass':
            reference_name = "reference's lengths.sizes"
            (inputs/'reference.sizes').rename(inputs/reference_name)
            command += ['--reference_name',reference_name]
        if scenario == 'missing':
            command += ['--missing_negative','true']
        env = dict(os.environ,PYTHONDONTWRITEBYTECODE='1',NXF_OFFLINE='true',NXF_DISABLE_CHECK_LATEST='true')
        result = subprocess.run(command,cwd=root,env=env,capture_output=True,text=True,timeout=180)
        (root/'console.txt').write_text(result.stdout+result.stderr)
        summary_path = root/'published/datachecks/summary/summary.json'
        self.assertTrue(summary_path.exists(),result.stdout+result.stderr)
        summary = json.loads(summary_path.read_text())
        self.assertEqual(summary['expected_checks'],20)
        self.assertEqual(summary['received_status_files'],20)
        for entry in summary['results']:
            if entry['check']=='reference_coordinates':
                report_folder = (root/'published/datachecks'/entry['report']).parent
                self.assertFalse((report_folder/'reference_inputs').exists())
                self.assertFalse(any(p.suffix in ('.gff','.gff3') for p in report_folder.rglob('*')))
            for field in ('report','stdout','stderr'):
                self.assertTrue((root/'published/datachecks'/entry[field]).is_file(),entry)
        return result,summary

    def test_all_ten_checks_run_and_publish_before_success(self):
        result,summary = self.run_workflow('pass')
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        self.assertEqual(summary['counts'],{'PASS':20})

    def test_failures_and_input_errors_finish_all_checks_before_gate_fails(self):
        result,summary = self.run_workflow('fail')
        self.assertNotEqual(result.returncode,0)
        self.assertEqual(summary['exit_code'],2)
        self.assertEqual(summary['counts'],{'PASS':18,'FAIL':1,'INPUT_ERROR':1})
        self.assertIn('DATACHECK_GATE',result.stdout+result.stderr)

    def test_new_read_driven_partial_cds_overlap_still_fails_the_workflow(self):
        result,summary = self.run_workflow('read_overlap')
        self.assertNotEqual(result.returncode,0)
        self.assertEqual(summary['exit_code'],1)
        overlap = next(r for r in summary['results'] if r['strand']=='positive' and r['check']=='cds_overlaps')
        self.assertEqual(overlap['status'],'FAIL')
        self.assertEqual(next(r for r in summary['results'] if r['check']=='partial_protection' and r['strand']=='positive')['status'],'FAIL')
        self.assertEqual(summary['counts'],{'PASS':18,'FAIL':2})

    def test_missing_partition_reaches_completeness(self):
        result,summary = self.run_workflow('missing')
        self.assertNotEqual(result.returncode,0)
        completeness = next(r for r in summary['results'] if r['strand']=='negative' and r['check']=='completeness')
        self.assertEqual(completeness['status'],'FAIL')
        self.assertEqual(summary['counts'],{'PASS':10,'INPUT_ERROR':9,'FAIL':1})

    def test_valid_empty_partition_is_explicitly_not_applicable(self):
        result,summary = self.run_workflow('empty')
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        self.assertEqual(summary['counts'],{'PASS':14,'NOT_APPLICABLE':6})

    def test_main_entrypoint_reference_guard_and_explicit_disable(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            env = dict(os.environ,NXF_OFFLINE='true',NXF_DISABLE_CHECK_LATEST='true')
            secrets=root/"database's connection.json"
            secrets.write_text(json.dumps(dict(host='unused',port=3306,user='unused',password='',database='unused')))
            def run(*arguments):
                return subprocess.run([NEXTFLOW,'run',str(ROOT/'main.nf'),'-preview','--db_secrets',str(secrets),*arguments],
                    cwd=root,env=env,capture_output=True,text=True,timeout=90)
            result = run()
            self.assertNotEqual(result.returncode,0)
            self.assertIn('--datacheck_reference_lengths',result.stdout+result.stderr)
            result = run('--run_datachecks','false')
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            (root/'reference.sizes').write_text('1\t2000\n')
            result = run('--datacheck_reference_lengths',str(root/'reference.sizes'))
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)

    def test_duplicate_review_is_published_without_failing_workflow(self):
        result,summary=self.run_workflow('duplicates')
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        self.assertEqual(summary['counts'],{'PASS':19,'REVIEW':1})
        self.assertEqual(summary['status'],'REVIEW')
        duplicate=next(r for r in summary['results'] if r['check']=='duplicate_transcripts' and r['strand']=='positive')
        self.assertEqual(duplicate['checker_exit_code'],1)
        self.assertEqual(duplicate['exit_code'],0)
