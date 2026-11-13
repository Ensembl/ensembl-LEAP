import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from test_check_duplicate_transcripts import annotation, SCRIPT


def model(tid, tags='', leap='', gene='g1', strand='+'):
    text = annotation(tid, gene=gene, strand=strand)
    return ''.join(line.rstrip('\n') + (';tag=' + tags + ';leap_tag=' + leap if '\ttranscript\t' in line else '') + '\n'
                   for line in text.splitlines(True))


class DuplicateRecommendationTests(unittest.TestCase):
    def check(self, content, recommend=True):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            source,output=root/'input.gff3',root/'report.csv'
            source.write_text(content)
            original=source.read_bytes()
            command=[sys.executable,str(SCRIPT),str(source),str(output)]
            if recommend: command.append('--recommend-removals')
            result=subprocess.run(command,capture_output=True,text=True)
            self.assertEqual(result.returncode,1,result.stderr)
            self.assertEqual(source.read_bytes(),original)
            with output.open() as handle: rows={r['transcript_id']:r for r in csv.DictReader(handle)}
            return result,rows

    def test_mane_all_protected_and_override_conflicting_labels(self):
        for strand in ('+','-'):
            content=(model('ENST100.1','MANE_Select,Ensembl_canonical_extended',strand=strand)
                     +model('ENST200.1','MANE_Plus_Clinical',strand=strand)
                     +model('ENST001.1',leap='LEAP_blacklisted',strand=strand))
            result,rows=self.check(content)
            self.assertEqual({r['transcript_id'] for r in rows.values() if r['recommendation']=='keep'}, {'ENST100.1','ENST200.1'})
            self.assertEqual(rows['ENST001.1']['recommendation'],'remove_from_gff')
            self.assertEqual(rows['ENST001.1']['retained_member_ids'],'ENST100.1;ENST200.1')
            self.assertIn('REVIEW:',result.stdout)

    def test_noncanonical_before_blacklist_then_numeric_id(self):
        _,rows=self.check(model('ENST001.1','Ensembl_canonical_extended',leap='LEAP_blacklisted')+
                          model('ENST200.1')+model('ENST100.2'))
        self.assertEqual(rows['ENST100.2']['recommendation'],'keep')
        self.assertEqual(rows['ENST001.1']['reason'],'noncanonical_member_preferred')

    def test_all_blacklist_sources_and_aliases_are_preferred(self):
        for tag,leap in [('upstream_ATG',''),('CAGE_supported_TSS',''),('RP_supported_TIS',''),
                         ('upstream%20ATG',''),('cage-supported-tss',''),('','LEAP_blacklisted')]:
            with self.subTest(tag=tag,leap=leap):
                _,rows=self.check(model('ENST1.1')+model('ENST20.1',tag,leap))
                self.assertEqual(rows['ENST20.1']['recommendation'],'keep')
                self.assertEqual(rows['ENST20.1']['is_blacklisted'],'True')

    def test_canonical_only_retains_one_then_applies_blacklist(self):
        for tags,expected in [('', 'ENST1.1'), (';LEAP_blacklisted', 'ENST20.1')]:
            _,rows=self.check(model('ENST1.1','Ensembl_canonical_extended')+
                              model('ENST20.1','ensembl_canonical_extended',tags.strip(';')))
            self.assertEqual({r['transcript_id'] for r in rows.values() if r['recommendation']=='keep'},{expected})
            self.assertTrue(all(r['reason'].startswith('canonical_only_group:') for r in rows.values()))

    def test_rank_is_numeric_and_versions_do_not_determine_preference(self):
        _,rows=self.check(model('ENST10.1')+model('ENST2.99'))
        self.assertEqual(rows['ENST2.99']['recommendation'],'keep')
        for ids in [('ENST2.1','ENST2.99'),('arbitrary','ENST1.1'),('ENST2.1','ENSMUST2.1')]:
            _,rows=self.check(model(ids[0])+model(ids[1]))
            self.assertEqual({r['recommendation'] for r in rows.values()},{'review'})
            self.assertEqual({r['retained_member_ids'] for r in rows.values()},{''})

    def test_cross_gene_same_rules_and_reordering_stable(self):
        first=model('ENST1.1',gene='other')
        second=model('ENST20.1','MANE_select',gene='g1')
        _,rows=self.check(first+second)
        _,reordered=self.check(second+first)
        self.assertEqual(rows['ENST1.1']['cross_gene_group'],'True')
        self.assertEqual(rows['ENST1.1']['recommendation'],'remove_from_gff')
        for tid in rows:
            for field in ('recommendation','reason','group_id','retained_member_ids'):
                self.assertEqual(rows[tid][field],reordered[tid][field])

    def test_detection_only_keeps_recommendation_fields_empty(self):
        result,rows=self.check(model('ENST1.1')+model('ENST20.1','MANE_Select'),recommend=False)
        self.assertEqual({r['recommendation'] for r in rows.values()},{''})
        self.assertEqual({r['retained_member_ids'] for r in rows.values()},{''})
        self.assertIn('FAIL:',result.stdout)


if __name__=='__main__': unittest.main()
