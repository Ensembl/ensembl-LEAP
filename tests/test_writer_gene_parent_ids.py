"""Exercise the real writer without requiring pandas or read-processing tools."""
import csv
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
WRITER=ROOT/'bin/build_chromosome_gff_and_stats.py'
STRUCTURE=ROOT/'assets/check_gff_structural_integrity.py'


def annotation(tid='T.1', strand='+', parent='G.1_PAR_Y', start=200, end=450):
    tx_id=tid+'_PAR_Y' if parent and parent.endswith('_PAR_Y') else tid
    def feature(kind,a,b,number=None):
        attrs=dict(gene_id='G.1',transcript_id=tid,gene_type='protein_coding',
                   gene_name='G',transcript_type='protein_coding')
        if kind=='transcript':
            attrs['ID']=tx_id
            if parent is not None:attrs['Parent']=parent
        else:
            attrs['ID']=(kind+':'+tx_id+(':'+str(number) if kind=='exon' else ''))
            attrs['Parent']=tx_id
            if number is not None:attrs['exon_number']=str(number)
        return '\t'.join(['Y','test',kind,str(a),str(b),'.',strand,
                          '0' if kind=='CDS' else '.',
                          ';'.join(k+'='+v for k,v in attrs.items())])+'\n'
    rows=[feature('transcript',start,end)]
    for number,(a,b) in enumerate([(start,250),(400,end)] if strand=='+' else [(400,end),(start,250)],1):
        rows.append(feature('exon',a,b,number))
        ca,cb=(220,240) if a<300 else (400,420)
        rows.append(feature('CDS',ca,cb,number))
        if a<300:rows.append(feature('five_prime_UTR' if strand=='+' else 'three_prime_UTR',start,219,number))
        else:rows.append(feature('three_prime_UTR' if strand=='+' else 'five_prime_UTR',421,end,number))
    return ''.join(rows)


def records(path):
    result=[]
    for line in path.read_text().splitlines():
        f=line.split('\t');a=dict(v.split('=',1) for v in f[8].split(';') if '=' in v)
        result.append((f,a))
    return result


class WriterGeneParentTests(unittest.TestCase):
    def build(self,folder,text,strand='+'):
        sign='positive' if strand=='+' else 'negative'
        raw=folder/'raw.gff';raw.write_text(text)
        partial=folder/'partial.gff';partial.touch()
        blacklist=folder/'blacklist';blacklist.touch()
        fe=folder/(sign+'_FE.bed');le=folder/(sign+'_LE.bed')
        def bed(start,end,options,count):return f'Y\t{start}\t{end}\t{options}\t1\t{strand}\t.\t.\t{count}\tG.1\n'
        fe.write_text(bed(150,250,'185/200',3) if strand=='+' else bed(400,500,'470/450',3))
        le.write_text(bed(400,500,'470/450',4) if strand=='+' else bed(150,250,'185/200',4))
        gff,stats=folder/'extended.gff',folder/'stats.tsv'
        result=subprocess.run([sys.executable,str(WRITER),str(fe),str(le),str(raw),str(blacklist),str(gff),str(stats),'--partial-annotation',str(partial)],
                              cwd=folder,capture_output=True,text=True,env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',PYTHONHASHSEED='0'))
        return result,gff,stats

    def test_par_and_ordinary_parent_ids_both_strands(self):
        for strand in ('+','-'):
            for parent in ('G.1','G.1_PAR_Y','gene:G.1_PAR_Y'):
                with self.subTest(strand=strand,parent=parent),tempfile.TemporaryDirectory() as directory:
                    folder=Path(directory);text=annotation(strand=strand,parent=parent)
                    result,gff,stats=self.build(folder,text,strand)
                    self.assertEqual(result.returncode,0,result.stderr)
                    features=records(gff);gene=next((f,a) for f,a in features if f[2]=='gene')
                    self.assertEqual(gene[1]['ID'],parent)
                    self.assertEqual(gene[1]['gene_id'],'G.1')
                    self.assertEqual(gene[0][3:5],['150','500'])
                    tx=next((f,a) for f,a in features if f[2]=='transcript')
                    self.assertEqual(tx[1]['Parent'],parent)
                    self.assertEqual(tx[1]['transcript_id'],'T.1')
                    self.assertEqual(tx[1]['ID'],'T.1_PAR_Y' if parent.endswith('_PAR_Y') else 'T.1')
                    for f,a in features:
                        if f[2] not in ('gene','transcript'):self.assertEqual(a['Parent'],tx[1]['ID'])
                    with stats.open() as h:rows=list(csv.DictReader(h,delimiter='\t'))
                    self.assertEqual(rows[0]['gene_id'],'G.1')
                    self.assertEqual(rows[0]['five_prime_extension'],'50')
                    self.assertEqual(rows[0]['three_prime_extension'],'50')
                    self.assertEqual(rows[0]['five_prime_read_count'],'3')
                    self.assertEqual(rows[0]['three_prime_read_count'],'4')
                    check=subprocess.run([sys.executable,str(STRUCTURE),str(gff),str(folder/'issues.csv')],capture_output=True,text=True)
                    self.assertEqual(check.returncode,0,check.stdout+check.stderr)

    def test_gene_grouping_uses_parent_identity_and_expands_bounds(self):
        with tempfile.TemporaryDirectory() as directory:
            folder=Path(directory)
            # Distinct GFF parent identities must not be merged via stable gene_id.
            text=annotation('T.1',parent='G.1_PAR_Y')+annotation('T.2',parent='gene:G.1_PAR_Y')
            result,gff,stats=self.build(folder,text)
            self.assertEqual(result.returncode,0,result.stderr)
            genes=[(f,a) for f,a in records(gff) if f[2]=='gene']
            self.assertEqual({a['ID'] for f,a in genes},{'G.1_PAR_Y','gene:G.1_PAR_Y'})
            self.assertEqual({a['gene_id'] for f,a in genes},{'G.1'})
            self.assertTrue(all(f[3:5]==['150','500'] for f,a in genes))
        with tempfile.TemporaryDirectory() as directory:
            folder=Path(directory)
            text=annotation('T.1')+annotation('T.2',start=190,end=460)
            result,gff,stats=self.build(folder,text)
            self.assertEqual(result.returncode,0,result.stderr)
            genes=[(f,a) for f,a in records(gff) if f[2]=='gene']
            self.assertEqual(len(genes),1)
            self.assertEqual(genes[0][1]['ID'],'G.1_PAR_Y')

    def test_legacy_fallback_and_invalid_gene_parents(self):
        with tempfile.TemporaryDirectory() as directory:
            result,gff,stats=self.build(Path(directory),annotation(parent=None))
            self.assertEqual(result.returncode,0,result.stderr)
            self.assertEqual(next(a['ID'] for f,a in records(gff) if f[2]=='gene'),'G.1')
        for parent in ('','.','g1,g2'):
            with self.subTest(parent=parent),tempfile.TemporaryDirectory() as directory:
                result,gff,stats=self.build(Path(directory),annotation(parent=parent))
                self.assertNotEqual(result.returncode,0)
                self.assertIn('Expected one gene Parent',result.stderr)
                self.assertFalse(gff.exists())
                self.assertFalse(stats.exists())


if __name__=='__main__':unittest.main()
