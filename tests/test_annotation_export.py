import contextlib
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'bin'))
from annotation_db import load_connection
import get_current_annotation as prepare
from gxf_dump.exporter import TranscriptRecord, rows_to_records, records_to_features
from test_check_duplicate_transcripts import annotation


class AnnotationExportTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        self.secrets=self.root/'connection.json'
        self.values=dict(host='private-host',port=3306,user='private-user',password='PRIVATE_PASSWORD',database='private-db')
        self.secrets.write_text(json.dumps(self.values))

    def test_json_validation_aliases_and_no_credential_diagnostics(self):
        self.assertEqual(load_connection(self.secrets)['database'],'private-db')
        for change in ({'port':0},{'port':True},{'port':65536},{'port':'bad'}, {'host':''},
                       {'user':[]},{'database':None},{'password':None},{'dialect':'sqlite'},
                       {'dbname':'PRIVATE_ALIAS_DATABASE'},{'pass':'PRIVATE_ALIAS_PASSWORD'}):
            self.secrets.write_text(json.dumps(dict(self.values,**change)))
            result=subprocess.run([sys.executable,str(ROOT/'bin/annotation_db.py'),str(self.secrets)],capture_output=True,text=True)
            self.assertEqual(result.returncode,2,result.stderr)
            for secret in ('PRIVATE_PASSWORD','private-user','private-db','private-host','PRIVATE_ALIAS_DATABASE','PRIVATE_ALIAS_PASSWORD'):
                self.assertNotIn(secret,result.stdout+result.stderr)
        for text in ('[1,2]', '{"password":"PRIVATE_PASSWORD",','null'):
            self.secrets.write_text(text)
            with self.assertRaises(ValueError): load_connection(self.secrets)
        self.secrets.write_text(json.dumps(dict(host='h',port='3306',user='u',dbname='db',**{'pass':''})))
        self.assertEqual(load_connection(self.secrets)['password'],'')

    def test_filter_whole_ott_and_readthrough_models_keep_ott_exon_ids(self):
        source=self.root/'source.gff3';target=self.root/'final.gff3'
        content=annotation('ENST1.1',gene='ENSG1.1')+annotation('OTTHUMT2.1',gene='ENSG1.1')+annotation('ENST3.1',gene='OTTG1.1')+annotation('ENST4.1',gene='ENSG4.1')
        content='1\ttest\tgene\t100\t350\t.\t+\t.\tID=OTTG1.1;gene_id=OTTG1.1\n'+content.replace('exon_id=e','exon_id=OTTHUME')
        source.write_text(content)
        prepare.filter_annotation(source,target,{'ENSG4'},remove_ott=True)
        text=target.read_text()
        self.assertIn('ENST1.1',text)
        self.assertIn('OTTHUME',text)
        for tid in ('OTTHUMT2.1','ENST3.1','ENST4.1'): self.assertNotIn(tid,text)
        self.assertEqual(source.read_text(),content)

    def test_workflow_outputs_order_and_credential_safe_errors(self):
        def export(args):
            self.assertIn('--vega',args)
            Path(args[args.index('--output')+1]).write_text(annotation('ENST1.1')+annotation('ENST2.1',gene='readthrough.1'))
        with patch.object(prepare.exporter,'export_main',side_effect=export), patch.object(prepare,'exclusion_ids',return_value=({'readthrough'},{'ENST1'})):
            with contextlib.redirect_stdout(io.StringIO()): prepare.run(self.secrets,'21',self.root)
        result=(self.root/'current_annotation_21.gff3').read_text()
        self.assertEqual(result,annotation('ENST1.1'))
        self.assertEqual((self.root/'do_not_extend_transcripts_positive.txt').read_text(),'stable_id\nENST1\n')
        self.assertEqual((self.root/'do_not_extend_transcripts_negative.txt').read_text(),'stable_id\nENST1\n')
        previous=result
        def fail(args):
            print('PRIVATE_PASSWORD',file=sys.stderr)
            raise RuntimeError('private-user PRIVATE_PASSWORD')
        output=io.StringIO()
        with patch.object(prepare.exporter,'export_main',side_effect=fail), contextlib.redirect_stderr(output):
            with self.assertRaisesRegex(ValueError,'Database export/filter query failed'):
                prepare.run(self.secrets,'21',self.root)
        self.assertNotIn('PRIVATE_PASSWORD',output.getvalue())
        self.assertEqual((self.root/'current_annotation_21.gff3').read_text(),previous)
        self.assertFalse(list(self.root.glob('.annotation-*')))

    def test_bad_chromosome_and_empty_filtered_annotation(self):
        with self.assertRaisesRegex(ValueError,'Invalid chromosome'): prepare.run(self.secrets,'../21',self.root)
        def export(args): Path(args[-1]).write_text('##gff-version 3\n')
        with patch.object(prepare.exporter,'export_main',side_effect=export), patch.object(prepare,'exclusion_ids',return_value=(set(),set())):
            with self.assertRaisesRegex(ValueError,'empty'): prepare.run(self.secrets,'21',self.root)
        self.assertFalse((self.root/'current_annotation_21.gff3').exists())


    def test_bundled_cli_suppresses_malformed_secret_contents(self):
        self.secrets.write_text('{"password":"PRIVATE_PASSWORD",')
        result=subprocess.run([sys.executable,str(ROOT/'bin/gxf_dump.py'),'--db-secrets',str(self.secrets),
            '--chromosome','21','--output',str(self.root/'never.gff3')],capture_output=True,text=True)
        self.assertEqual(result.returncode,2)
        self.assertNotIn('PRIVATE_PASSWORD',result.stdout+result.stderr)
        self.assertNotIn('Traceback',result.stderr)
        self.assertFalse((self.root/'never.gff3').exists())

if __name__=='__main__': unittest.main()
