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

from test_prepare_manual_review import fixture

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('review_sites',ROOT/'bin/prepare_review_sites.py')
sites=importlib.util.module_from_spec(spec);spec.loader.exec_module(sites)


def write_sites(root):
    root.mkdir(parents=True,exist_ok=True)
    cage=root/'CAGE.bed'
    cage.write_text('chr1\t0\t10\tpeak_plus\t2972\t+\t2\t3\t255,0,0\n'
                    'chr1\t30\t40\tpeak_minus\t397\t-\t33\t34\t255,0,0\n'
                    'chr1\t50\t60\tunknown_strand\t0\t.\t53\t54\t0\n'
                    'chr2\t10\t20\tother_chromosome\t10\t+\t12\t13\t0\n')
    polya=root/'polyA.bed'
    polya.write_text('1\t20\t21\tpoly_plus\t0.2690\t+\t0.15\t4\t0.2690\tIN\tCATAAA@-18@16467\n'
                     '1\t60\t61\tpoly_minus\t0.0064\t-\t0.00\t1\t0.0064\tIG\tATTAAA@-24@77654\n')
    return dict(cage=cage,polya=polya)


@unittest.skipUnless(shutil.which('sort'),'requires local sort')
class SiteEvidenceTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        self.sources=write_sites(self.root/'inputs')
        self.scopes=[('1','positive'),('1','negative')]

    def test_bed_views_keep_coordinates_multiplicity_and_native_polya_fields(self):
        before={p:p.read_bytes() for p in self.sources.values()}
        summary=sites.prepare_sites(self.sources,self.scopes,self.root/'sites')
        self.assertEqual({p:p.read_bytes() for p in before},before)
        self.assertEqual(summary['sources']['cage']['region_format'],'cageBed9')
        self.assertEqual(summary['sources']['cage']['counts']['unsplit'],3)
        self.assertEqual(summary['sources']['cage']['counts']['unstranded'],1)
        self.assertEqual(summary['sources']['cage']['counts']['excluded'],1)
        rows=(self.root/'sites/cage/positive.bed').read_text().splitlines()
        self.assertEqual(rows,['1\t0\t10\tpeak_plus\t1000\t+\t2\t3\t255,0,0\t2972'])
        self.assertTrue((self.root/'sites/region_presets/cageBed9.tsv').is_file())
        self.assertEqual(summary['sources']['polya']['region_format'],'polyAtlas')
        self.assertEqual((self.root/'sites/polya/unsplit.bed').read_text(),self.sources['polya'].read_text())
        with (self.root/'sites/cage/source_records.tsv').open() as handle:
            rows=list(csv.DictReader(handle,delimiter='\t'))
        self.assertEqual(json.loads(rows[0]['source_fields_json'])[4],'2972')
        self.sources['cage'].write_text(self.sources['cage'].read_text()*2)
        summary=sites.prepare_sites({'cage':self.sources['cage']},self.scopes,self.root/'duplicate')
        self.assertEqual(summary['sources']['cage']['counts']['unsplit'],6)

    def test_gzip_gff_conversion_and_restricted_strand_scope(self):
        source=self.root/'CAGE.gff3.gz'
        with gzip.open(source,'wt') as handle:
            handle.write('##gff-version 3\n1\tFANTOM\tCAGE_cluster\t1\t10\t1.25\t+\t.\tID=p;Name=peak%20one\n'
                         '1\tFANTOM\tCAGE_cluster\t31\t40\t.\t-\t.\tID=n\n')
        summary=sites.prepare_sites({'cage':source},[('1','positive')],self.root/'gff')
        self.assertEqual(summary['sources']['cage']['region_format'],'cageBed6')
        self.assertEqual((self.root/'gff/cage/positive.bed').read_text(),'1\t0\t10\tpeak one\t1\t+\t1.25\n')
        self.assertEqual((self.root/'gff/cage/negative.bed').read_text(),'')

    def test_empty_source_and_single_evidence_source(self):
        self.sources['cage'].write_text('# no sites\n')
        summary=sites.prepare_sites({'cage':self.sources['cage']},self.scopes,self.root/'empty')
        self.assertEqual(summary['sources']['cage']['counts']['selected'],0)
        self.assertEqual(len(list((self.root/'empty/cage').glob('*.bed'))),3)

    def test_invalid_records_reference_scope_and_output_protection(self):
        original=self.sources['cage'].read_text()
        for bad in ('1\t-1\t5\tx\t0\t+\n','1\t5\t5\tx\t0\t+\n','1\t1\t5\tx\t0\t!\n',
                    '1\t1\t5\tx\tNaN\t+\n','invalid\n'):
            self.sources['cage'].write_text(bad)
            with self.assertRaisesRegex(ValueError,'malformed'):
                sites.prepare_sites(self.sources,self.scopes,self.root/'bad')
            self.assertFalse((self.root/'bad').exists())
        self.sources['cage'].write_text(original)
        lengths=self.root/'ref.fa.gz.fai';lengths.write_text('1\t35\t0\t60\t61\n')
        with self.assertRaisesRegex(ValueError,'reference length'):
            sites.prepare_sites(self.sources,self.scopes,self.root/'bad',lengths)
        with self.assertRaisesRegex(ValueError,'Duplicate'):
            sites.prepare_sites(self.sources,self.scopes+self.scopes,self.root/'bad')
        with self.assertRaisesRegex(ValueError,'new directory'):
            sites.prepare_sites(self.sources,self.scopes,self.root)
        link=self.root/'linked';link.symlink_to(self.root,target_is_directory=True)
        with self.assertRaisesRegex(ValueError,'new directory'):
            sites.prepare_sites(self.sources,self.scopes,link)

    def test_sort_failures_leave_no_published_directory(self):
        with patch.object(sites.subprocess,'run',side_effect=subprocess.CalledProcessError(1,'sort')):
            with self.assertRaises(subprocess.SubprocessError):
                sites.prepare_sites(self.sources,self.scopes,self.root/'bad')
        self.assertFalse((self.root/'bad').exists())
        self.assertFalse(list(self.root.glob('.bad-*')))

    def test_standalone_bundle_registers_public_offline_and_empty_views(self):
        run=fixture(self.root/'leap')
        output=self.root/'review'
        result=subprocess.run([sys.executable,str(ROOT/'bin/prepare_manual_review.py'),'run',str(run),str(output),
            '--chromosomes','1','--genome-id','hg38','--cage',str(self.sources['cage']),'--polya',str(self.sources['polya'])],
            capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stderr)
        for filename in ('samplesheet.csv','samplesheet_offline.csv'):
            with (output/filename).open() as handle:
                rows=[r for r in csv.DictReader(handle) if r['type']=='regions']
            self.assertEqual(len(rows),6)
            self.assertTrue(all(r['species']=='hg38' and r['data_category']=='regions' for r in rows))
            self.assertTrue(all(Path(r['file']).is_file() for r in rows))
            if filename.endswith('_offline.csv'):
                self.assertTrue(all(r['publication']=='offline' for r in rows))
        summary=json.loads((output/'review-summary.json').read_text())
        self.assertEqual(summary['site_evidence']['status'],'available')
        self.assertEqual(summary['site_evidence']['sources']['cage']['counts']['unstranded'],1)
        self.assertTrue((output/'region_presets/cageBed9.tsv').is_file())
        self.assertTrue((output/'install_review_region_presets.py').is_file())
        for script in ('run_jbrowser.sbatch','run_jbrowser_offline.sbatch'):
            self.assertIn('install_review_region_presets.py',(output/script).read_text())

    def test_original_scores_are_lossless_for_each_strand_and_unsplit(self):
        self.sources['cage'].write_text('chr1\t0\t5\ta\t2000\t+\nchr1\t10\t15\tb\t20000\t-\n'
                                       'chr1\t20\t25\tc\t0.1250\t+\nchr1\t30\t35\td\t.\t.\n')
        summary=sites.prepare_sites({'cage':self.sources['cage']},self.scopes,self.root/'scores')
        self.assertEqual(summary['sources']['cage']['region_format'],'cageBed6')
        for view,wanted in [('positive',['2000','0.1250']),('negative',['20000']),('unsplit',['2000','20000','0.1250','.'])]:
            rows=[line.split('\t') for line in (self.root/'scores/cage'/(view+'.bed')).read_text().splitlines()]
            self.assertEqual([r[-1] for r in rows],wanted)
            self.assertTrue(all(0<=int(r[4])<=1000 for r in rows))

    def test_preset_installation_is_repeatable_and_refuses_conflicts(self):
        spec=importlib.util.spec_from_file_location('install_presets',ROOT/'bin/install_review_region_presets.py')
        installer=importlib.util.module_from_spec(spec);spec.loader.exec_module(installer)
        bundle=self.root/'bundle';bundle.mkdir()
        sites.prepare_sites({'cage':self.sources['cage']},self.scopes,bundle/'evidence')
        shutil.copytree(bundle/'evidence/region_presets',bundle/'region_presets')
        (bundle/'samplesheet.csv').write_text('sample,region_format\ncage,cageBed9\n')
        jbrowser=self.root/'jbrowser';jbrowser.mkdir();(jbrowser/'main.nf').write_text('// fixture\n')
        for _ in range(2):
            self.assertEqual(installer.install_presets(bundle,jbrowser),['cageBed9'])
        destination=jbrowser/'assets/region_presets/cageBed9.tsv'
        self.assertEqual(destination.read_bytes(),(bundle/'region_presets/cageBed9.tsv').read_bytes())
        destination.write_text('different preset\n')
        with self.assertRaisesRegex(ValueError,'Conflicting'):
            installer.install_presets(bundle,jbrowser)
        self.assertEqual(destination.read_text(),'different preset\n')

    @unittest.skipUnless(os.environ.get('LEAP_JBROWSER_PIPELINE'),'set LEAP_JBROWSER_PIPELINE for region validation')
    def test_jbrowser_accepts_each_nonempty_region_view(self):
        summary=sites.prepare_sites(self.sources,self.scopes,self.root/'sites')
        # Test against an isolated JBrowser checkout with our bundled definitions.
        browser=self.root/'jbrowser';(browser/'bin').mkdir(parents=True)
        source=Path(os.environ['LEAP_JBROWSER_PIPELINE'])
        helper=browser/'bin/regions_bigbed_lib.py'
        shutil.copyfile(source/'bin/regions_bigbed_lib.py',helper)
        shutil.copytree(source/'assets/region_presets',browser/'assets/region_presets')
        for preset in (self.root/'sites/region_presets').glob('*.tsv'):
            shutil.copyfile(preset,browser/'assets/region_presets'/preset.name)
        for kind,source in summary['sources'].items():
            for view in sites.VIEWS:
                prefix=self.root/(kind+'_'+view)
                result=subprocess.run([sys.executable,str(helper),'validate','--input',str(self.root/'sites'/kind/(view+'.bed')),
                    '--sample',kind+'_'+view,'--region-format',source['region_format'],
                    '--validation-tsv',str(prefix)+'.tsv','--mqc-json',str(prefix)+'.json','--status-file',str(prefix)+'.status'],
                    capture_output=True,text=True)
                self.assertEqual(result.returncode,0,result.stdout+result.stderr)
                self.assertIn(Path(str(prefix)+'.status').read_text().strip(),('PASS','WARN'))
        spec=importlib.util.spec_from_file_location('jbrowser_regions',helper)
        regions=importlib.util.module_from_spec(spec);spec.loader.exec_module(regions)
        _,records,_,issues,schema=regions.validate_region_records(str(self.root/'sites/cage/positive.bed'),'cageBed9')
        self.assertEqual(regions.issue_counts(issues)[0],0)
        self.assertEqual(schema['bigbed_type'],'bed9+1')
        self.assertEqual(schema['fields'][-1][1],'originalCageScore')
        self.assertEqual(records[0][3][-1],'2972')


if __name__=='__main__':unittest.main()
