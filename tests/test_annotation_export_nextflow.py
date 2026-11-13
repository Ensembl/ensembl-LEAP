import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

from test_check_duplicate_transcripts import annotation

ROOT=Path(__file__).resolve().parents[1]
NEXTFLOW=os.environ.get('LEAP_NEXTFLOW') or shutil.which('nextflow')


@unittest.skipUnless(NEXTFLOW,'requires Nextflow')
class AnnotationExportNextflowTests(unittest.TestCase):
    def test_real_module_stages_secret_path_and_publishes_only_annotation(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);helpers=root/'helpers';helpers.mkdir()
            secrets=root/"private connection's config.json"
            values=dict(host='PRIVATE_HOST_VALUE',port=3306,user='PRIVATE_USER_VALUE',password='PRIVATE_PASSWORD_VALUE',database='PRIVATE_DATABASE_VALUE')
            secrets.write_text(json.dumps(values));before=secrets.read_bytes()
            shutil.copyfile(ROOT/'bin/get_current_annotation.py',helpers/'prepare_annotation.py')
            shutil.copyfile(ROOT/'bin/annotation_db.py',helpers/'annotation_db.py')
            shutil.copytree(ROOT/'bin/gxf_dump',helpers/'gxf_dump',ignore=shutil.ignore_patterns('__pycache__'))
            (helpers/'fixture.gff3').write_text(annotation('ENST1.1'))
            (helpers/'get_current_annotation.py').write_text('''from pathlib import Path
import prepare_annotation as task

def replay(args):
    Path(args[-1]).write_text((Path(__file__).parent/'fixture.gff3').read_text())

task.exporter.export_main=replay
task.exclusion_ids=lambda secrets, chromosome: (set(), {'ENST1'})
raise SystemExit(task.main())
''')
            config=root/'local.config'
            config.write_text("process.executor='local'\nprocess.cpus=1\nprocess.memory='1 GB'\nconda.enabled=false\nsingularity.enabled=false\n")
            command=[NEXTFLOW,'-C',str(config),'run',str(ROOT/'tests/nextflow/annotation_export.nf'),
                     '-ansi-log','false','-work-dir',str(root/'work'),'--db_secrets',str(secrets),
                     '--helpers',str(helpers),'--output_directory',str(root/'published')]
            env=dict(os.environ,NXF_OFFLINE='true',NXF_DISABLE_CHECK_LATEST='true',PYTHONDONTWRITEBYTECODE='1')
            result=subprocess.run(command,cwd=root,env=env,capture_output=True,text=True,timeout=120)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            published=list((root/'published').rglob('*'))
            self.assertEqual([p.name for p in published if p.is_file()],['current_annotation_21.gff3'])
            self.assertEqual((root/'published/raw_annotation/current_annotation_21.gff3').read_text(),annotation('ENST1.1'))
            self.assertEqual(secrets.read_bytes(),before)
            scripts=list((root/'work').rglob('.command.sh'))
            self.assertEqual(len(scripts),1)
            script=scripts[0].read_text()
            self.assertIn('connection/db_secrets.json',script)
            for value in values.values():
                if isinstance(value,str): self.assertNotIn(value,script+result.stdout+result.stderr)
            self.assertEqual((scripts[0].parent/'do_not_extend_transcripts_negative.txt').read_text(),'stable_id\nENST1\n')

    def test_entrypoint_json_and_legacy_guards_do_not_echo_credentials(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);path=root/'connection.json'
            good=dict(host='h',port=3306,user='u',password='DO_NOT_PRINT_PASSWORD',database='db')
            env=dict(os.environ,NXF_OFFLINE='true',NXF_DISABLE_CHECK_LATEST='true')
            for content, extra, expected in (
                (json.dumps(good),['--pass','DO_NOT_PRINT_PASSWORD'],'Legacy database parameters'),
                ('{"password":"DO_NOT_PRINT_PASSWORD",',[],'Invalid database secrets JSON'),
                (json.dumps(dict(good,port=0)),[],'Invalid database secrets JSON'),
                ('[]',[],'Invalid database secrets JSON'),
            ):
                path.write_text(content)
                command=[NEXTFLOW,'run',str(ROOT/'main.nf'),'-preview','--run_datachecks','false','--db_secrets',str(path),*extra]
                result=subprocess.run(command,cwd=root,env=env,capture_output=True,text=True,timeout=90)
                self.assertNotEqual(result.returncode,0,result.stdout+result.stderr)
                self.assertIn(expected,result.stdout+result.stderr)
                self.assertNotIn('DO_NOT_PRINT_PASSWORD',result.stdout+result.stderr)
            result=subprocess.run([NEXTFLOW,'run',str(ROOT/'main.nf'),'-preview','--run_datachecks','false'],cwd=root,env=env,capture_output=True,text=True,timeout=90)
            self.assertNotEqual(result.returncode,0)
            self.assertIn('--db_secrets',result.stdout+result.stderr)


if __name__=='__main__': unittest.main()
