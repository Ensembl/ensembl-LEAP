import csv
import gzip
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / 'assets/check_gff_structural_integrity.py'


def row(kind, start, end, identifier=None, parent=None, tid=None, gene=None,
        strand='+', chromosome='1', phase=None):
    attrs = {}
    for name, value in [('ID', identifier), ('Parent', parent), ('transcript_id', tid), ('gene_id', gene)]:
        if value is not None:
            attrs[name] = value
    return '\t'.join([chromosome, 'test', kind, str(start), str(end), '.', strand,
                      phase if phase is not None else '0' if kind == 'CDS' else '.',
                      ';'.join(k+'='+v for k,v in attrs.items()) or '.']) + '\n'


def annotation(strand='+', transcript_type='transcript', tid='t1', gene='g1', with_gene=True):
    def feature(kind, start, end, identifier, parent, transcript=tid, phase=None):
        return row(kind,start,end,identifier,parent,transcript,gene,strand,phase=phase)
    result = row('gene',1,400,'gene:'+gene,gene=gene,strand=strand) if with_gene else ''
    result += feature(transcript_type,100,350,'transcript:'+tid,'gene:'+gene)
    result += feature('exon',100,150,'exon:'+tid+':1','transcript:'+tid)
    result += feature('exon',300,350,'exon:'+tid+':2','transcript:'+tid)
    result += feature('CDS',125,150,'CDS:'+tid,'transcript:'+tid,phase='0')
    result += feature('CDS',300,325,'CDS:'+tid,'transcript:'+tid,phase='1')
    result += feature('five_prime_UTR' if strand=='+' else 'three_prime_UTR',100,124,'UTR-low:'+tid,'transcript:'+tid)
    result += feature('three_prime_UTR' if strand=='+' else 'five_prime_UTR',326,350,'UTR-high:'+tid,'transcript:'+tid)
    return result


def edit(text, predicate, column=None, value=None, attrs=None):
    result = []
    for line in text.splitlines():
        fields = line.split('\t')
        if predicate(fields):
            if column is not None:
                fields[column] = str(value)
            if attrs is not None:
                values = dict(v.split('=',1) for v in fields[8].split(';') if '=' in v)
                values.update(attrs)
                fields[8] = ';'.join(k+'='+v for k,v in values.items())
        result.append('\t'.join(fields)+'\n')
    return ''.join(result)


class GffIntegrityTests(unittest.TestCase):
    def run_check(self, text, compressed=False, allow_missing_genes=False):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);source=root/('input.gff3.gz' if compressed else 'input.gff3');report=root/'issues.csv'
            with (gzip.open if compressed else open)(source,'wt') as handle:handle.write(text)
            command=[sys.executable,str(SCRIPT),str(source),str(report)]
            if allow_missing_genes:command.append('--allow-missing-genes')
            result=subprocess.run(command,capture_output=True,text=True)
            issues=None
            if report.exists():
                with report.open(newline='') as handle:issues=list(csv.DictReader(handle))
            return result,issues

    def assert_issue(self, text, code, **options):
        result,issues=self.run_check(text,**options)
        self.assertEqual(result.returncode,1,result.stderr)
        self.assertIn(code,{r['issue'] for r in issues},issues)
        return issues

    def test_valid_both_strands_order_gzip_multipart_and_mrna(self):
        for strand in ('+','-'):
            for kind in ('transcript','mRNA'):
                text=annotation(strand,kind)
                text='# comment\n\n'+''.join(reversed(text.splitlines(keepends=True)))+'##FASTA\n>1\nACGT\n'
                result,issues=self.run_check(text,compressed=True)
                self.assertEqual(result.returncode,0,result.stderr)
                self.assertEqual(issues,[])
                self.assertIn('PASS: 0 structural issue(s)',result.stdout)
                self.assertIn('gene_boundaries_checked: True',result.stdout)
                self.assertIn('transcripts: 1',result.stdout)

    def test_gene_and_transcript_contain_all_descendants(self):
        for strand in ('+','-'):
            text=annotation(strand)
            bad=edit(text,lambda f:f[2]=='gene',4,349)
            issues=self.assert_issue(bad,'outside_parent')
            self.assertTrue(any(r['issue']=='outside_ancestor' and r['feature_type']=='exon' and r['parent_id']=='gene:g1' for r in issues))
            self.assertTrue(all(r['line'].isdigit() for r in issues))
            bad=edit(text,lambda f:f[2]=='transcript',4,349)
            issues=self.assert_issue(bad,'outside_parent')
            self.assertTrue(any(r['feature_type']=='exon' for r in issues))
            # Containment allows margins; equality to the outer exon bounds is
            # already enforced separately by the stats/GFF consistency check.
            large=edit(text,lambda f:f[2]=='transcript',3,90)
            large=edit(large,lambda f:f[2]=='transcript',4,360)
            self.assertEqual(self.run_check(large)[0].returncode,0)

    def test_parent_ids_types_and_custom_identity_attributes(self):
        text=annotation()
        self.assert_issue(edit(text,lambda f:f[2]=='exon',attrs={'Parent':'absent'}),'unresolved_parent')
        self.assert_issue(edit(text,lambda f:f[2]=='transcript',attrs={'Parent':'exon:t1:1'}),'invalid_parent_type')
        self.assert_issue(edit(text,lambda f:f[2]=='exon',attrs={'Parent':'gene:g1'}),'invalid_parent_type')
        self.assert_issue(edit(text,lambda f:f[2]=='CDS',attrs={'transcript_id':'wrong'}),'transcript_id_mismatch')
        self.assert_issue(edit(text,lambda f:f[2]=='transcript',attrs={'gene_id':'wrong'}),'gene_id_mismatch')
        self.assert_issue(edit(text,lambda f:f[2]=='exon',attrs={'gene_id':'wrong'}),'gene_id_mismatch')
        missing=text.replace('ID=transcript:t1;','')
        self.assert_issue(missing,'missing_id')
        missing=text.replace('Parent=transcript:t1;','')
        self.assert_issue(missing,'missing_parent')
        self.assert_issue(edit(text,lambda f:f[2]=='transcript',attrs={'Parent':'gene:g1,other'}),'invalid_transcript_parent_count')
        self.assert_issue(edit(text,lambda f:f[2]=='gene',attrs={'Parent':'transcript:t1'}),'invalid_parent_type')
        # Y PAR GFF IDs carry a suffix even when stable gene_id attributes do not.
        par=text.replace('ID=gene:g1;', 'ID=gene:g1_PAR_Y;').replace('Parent=gene:g1;', 'Parent=gene:g1_PAR_Y;')
        self.assertEqual(self.run_check(par)[0].returncode,0)
        broken=par.replace('ID=gene:g1_PAR_Y;', 'ID=gene:g1;')
        issues=self.assert_issue(broken,'unresolved_parent')
        self.assertTrue(any(r['parent_id']=='gene:g1_PAR_Y' for r in issues))

    def test_shared_exons_multiple_parents_and_optional_leaf_ids(self):
        text=row('gene',1,400,'g',gene='g')
        for tid in ['t1','t2']:text+=row('transcript',100,350,tid,'g',tid,'g')
        for start,end in [(100,150),(300,350)]:text+=row('exon',start,end,parent='t1,t2',gene='g')
        text+=row('CDS',125,150,'cds','t1,t2',gene='g')+row('CDS',300,325,'cds','t1,t2',gene='g',phase='1')
        result,issues=self.run_check(text)
        self.assertEqual(result.returncode,0,result.stderr)
        self.assertEqual(issues,[])
        self.assertIn('transcripts: 2',result.stdout)
        self.assert_issue(edit(text,lambda f:f[2]=='exon',attrs={'transcript_id':'t1'}),'transcript_id_mismatch')
        self.assert_issue(edit(text,lambda f:f[2]=='exon',attrs={'Parent':'t1,t1'}),'duplicate_parent_reference')

    def test_gene_free_mode_is_explicit_and_does_not_hide_partial_missing_genes(self):
        text=annotation(with_gene=False)
        self.assert_issue(text,'unresolved_parent')
        result,issues=self.run_check(text,allow_missing_genes=True)
        self.assertEqual(result.returncode,0,result.stderr)
        self.assertEqual(issues,[])
        self.assertIn('gene_boundaries_checked: False',result.stdout)
        self.assertIn('omitted_gene_parent_checks: 1',result.stdout)
        self.assert_issue(edit(text,lambda f:f[2]=='exon',attrs={'Parent':'missing_tx'}),'unresolved_parent',allow_missing_genes=True)
        partial=annotation()+annotation(tid='t2',gene='g2',with_gene=False)
        self.assert_issue(partial,'unresolved_parent',allow_missing_genes=True)
        bad=text+annotation(tid='t2',gene='other',with_gene=False).replace('Parent=gene:other','Parent=gene:g1')
        self.assert_issue(bad,'gene_id_mismatch',allow_missing_genes=True)

    def test_duplicate_records_conflicting_ids_and_legitimate_split_ids(self):
        text=annotation();lines=text.splitlines(keepends=True)
        issues=self.assert_issue(text+lines[2],'duplicate_record')
        detail=json.loads(next(r['details'] for r in issues if r['issue']=='duplicate_record'))
        self.assertEqual(detail['first_line'],3)
        self.assert_issue(text+lines[1],'conflicting_id')
        self.assert_issue(edit(text,lambda f:f[2]=='exon',attrs={'ID':'transcript:t1'}),'conflicting_id')
        self.assert_issue(edit(text,lambda f:f[2]=='CDS' and f[3]=='300',0,'2'),'conflicting_id')
        self.assertEqual(self.run_check(text)[0].returncode,0)
        # Ordinary exon_id values are not global GFF IDs and may be reused.
        reused=edit(text,lambda f:f[2]=='exon',attrs={'exon_id':'same_stable_exon'})
        self.assertEqual(self.run_check(reused)[0].returncode,0)
        duplicate_tid=text+annotation(tid='t2',gene='g2').replace('transcript_id=t2','transcript_id=t1')
        self.assert_issue(duplicate_tid,'duplicate_transcript_id')
        for strand in ('+','-'):
            # GFF3 represents programmed frameshifts with overlapping CDS parts.
            frameshift=row('gene',1,400,'g',strand=strand)
            frameshift+=row('transcript',100,350,'t','g',strand=strand)
            frameshift+=row('exon',100,350,parent='t',strand=strand)
            frameshift+=row('CDS',125,250,'c','t',strand=strand)
            frameshift+=row('CDS',249,325,'c','t',strand=strand)
            self.assertEqual(self.run_check(frameshift)[0].returncode,0)

    def test_exon_integrity_and_containing_exons_for_cds_utrs_and_codons(self):
        text=annotation()
        no_exons=''.join(line+'\n' for line in text.splitlines() if line.split('\t')[2]!='exon')
        self.assert_issue(no_exons,'transcript_without_exons')
        self.assert_issue(edit(text,lambda f:f[2]=='exon' and f[3]=='100',4,300),'overlapping_exons')
        for kind in ['CDS','five_prime_UTR','start_codon','stop_codon','three_prime_UTR']:
            bad=text+row(kind,145,305,'crosses_'+kind,'transcript:t1','t1','g1')
            self.assert_issue(bad,'feature_outside_exon')
        self.assert_issue(edit(text,lambda f:f[2]=='CDS',7,'.'),'missing_cds_phase')
        self.assert_issue(text+row('gene',500,600,'empty_gene',gene='g2'),'gene_without_transcripts')

    def test_nested_unknown_features_cycles_and_discontinuous_parent_coverage(self):
        text=annotation()+row('custom_feature',125,127,'custom','CDS:t1','t1','g1')
        self.assertEqual(self.run_check(text)[0].returncode,0)
        bad=edit(text,lambda f:f[2]=='custom_feature',4,305)
        self.assert_issue(bad,'outside_parent')
        root=annotation()+row('region',1,500,'region')
        self.assertEqual(self.run_check(root)[0].returncode,0)
        cycle=annotation()+row('custom',200,210,'a','b')+row('custom',200,210,'b','a')
        self.assert_issue(cycle,'cyclic_parent_graph')
        self.assert_issue(annotation()+row('custom',200,210,'self','self'),'self_parent')
        encoded=annotation().replace('ID=CDS:t1','ID=CDS%2Ct1')
        encoded+=row('custom',125,127,'custom','CDS%2Ct1','t1','g1')
        self.assertEqual(self.run_check(encoded)[0].returncode,0)

    def test_chromosome_and_strand_mismatches(self):
        text=annotation()
        self.assert_issue(edit(text,lambda f:f[2]=='exon',0,'2'),'chromosome_mismatch')
        self.assert_issue(edit(text,lambda f:f[2]=='exon',6,'-'),'strand_mismatch')
        other=annotation().replace('1\ttest','2\ttest')
        self.assert_issue(text+other,'conflicting_id')

    def test_parse_errors_and_output_protection(self):
        text=annotation()
        cases=[('', 'no feature'),('# comment\n','no feature'),('1\ttest\tgene\n','9 GFF3 columns'),
               (row('gene',1,400,'g'),'no transcript'),
               (row('region',1,400,'r'),'no transcript'),
               (text.replace('100\t150','x\t150'),'integers'),
               (text.replace('100\t150','0\t150'),'invalid coordinates'),
               (text.replace('100\t150','151\t150'),'invalid coordinates'),
               (text.replace('\t+\t','\tx\t'),'invalid strand'),
               (text.replace('\t+\t0\t','\t+\t3\t'),'invalid strand or phase'),
               (edit(text,lambda f:f[2]=='gene',5,'nan'),'invalid score'),
               (text.replace('gene_id=g1','bad'),'malformed'),
               (text.replace('gene_id=g1','gene_id=g1;gene_id=g1'),'duplicate attribute'),
               (text.replace('ID=gene:g1','ID=a,b'),'ID must have one'),
               (text.replace('ID=gene:g1','ID=bad%ZZ'),'malformed identifier'),
               (text.replace('Parent=gene:g1','Parent='),'malformed identifier')]
        for invalid,message in cases:
            with self.subTest(message=message):
                result,issues=self.run_check(invalid)
                self.assertEqual(result.returncode,2,result.stderr)
                self.assertIn(message,result.stderr)
                self.assertIsNone(issues)
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);source=root/'input.gff3';source.write_text(text)
            alias=root/'alias.gff3';alias.symlink_to(source)
            hard=root/'hard.gff3';hard.hardlink_to(source)
            for output in [source,alias,hard,root/'missing/report.csv']:
                result=subprocess.run([sys.executable,str(SCRIPT),str(source),str(output)],capture_output=True,text=True)
                self.assertEqual(result.returncode,2,result.stderr);self.assertEqual(source.read_text(),text)
            report=root/'report.csv';report.write_text('previous report')
            for path in [root/'absent',root/'bad.gff3.gz']:
                if path.suffix=='.gz':path.write_bytes(b'not gzip')
                result=subprocess.run([sys.executable,str(SCRIPT),str(path),str(report)],capture_output=True,text=True)
                self.assertEqual(result.returncode,2);self.assertEqual(report.read_text(),'previous report')


if __name__=='__main__':
    unittest.main()
