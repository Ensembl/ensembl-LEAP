import argparse
from contextlib import contextmanager
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'bin'))
from gxf_dump.annotation import Feature
from gxf_dump.exporter import write_features
from gxf_dump.par import add_y_par_features
try:
    from sqlalchemy import create_engine
except ImportError:
    create_engine=None


def records(gid,tid,chrom,start,end,strand='+'):
    common=dict(gene_id=gid,gene_type='protein_coding')
    return [Feature(chrom,'test','gene',start,end,'.',strand,'.',dict(common,ID=gid)),
            Feature(chrom,'test','transcript',start,end,'.',strand,'.',dict(common,ID=tid,Parent=gid,transcript_id=tid,tag='MANE_Select')),
            Feature(chrom,'test','exon',start,end,'.',strand,'.',dict(common,ID='exon:'+tid+':1',Parent=tid,transcript_id=tid,exon_id='ENSE1.1',exon_number='1')),
            Feature(chrom,'test','CDS',start+3,end-3,'.',strand,'0',dict(common,ID='CDS:'+tid,Parent=tid,transcript_id=tid,protein_id='ENSP1.1'))]


@unittest.skipUnless(create_engine,'requires SQLAlchemy')
class PARExportTests(unittest.TestCase):
    def setUp(self):
        self.engine=create_engine('sqlite://');self.addCleanup(self.engine.dispose)
        self.db=self.engine.connect();self.addCleanup(self.db.close)
        for sql in '''
            CREATE TABLE coord_system (coord_system_id INTEGER,species_id INTEGER,name TEXT,version TEXT,attrib TEXT);
            INSERT INTO coord_system VALUES (1,1,'chromosome','GRCh38','default_version');
            CREATE TABLE meta (species_id INTEGER,meta_key TEXT,meta_value TEXT);
            INSERT INTO meta VALUES (1,'assembly.default','GRCh38');
            CREATE TABLE seq_region (seq_region_id INTEGER,name TEXT,length INTEGER,coord_system_id INTEGER);
            INSERT INTO seq_region VALUES (10,'X',2000,1),(20,'Y',2000,1);
            CREATE TABLE attrib_type (attrib_type_id INTEGER,code TEXT);
            INSERT INTO attrib_type VALUES (1,'toplevel');
            CREATE TABLE seq_region_attrib (seq_region_id INTEGER,attrib_type_id INTEGER,value TEXT);
            INSERT INTO seq_region_attrib VALUES (10,1,'1'),(20,1,'1');
            CREATE TABLE assembly_exception (seq_region_id INTEGER,seq_region_start INTEGER,seq_region_end INTEGER,
                exc_seq_region_id INTEGER,exc_seq_region_start INTEGER,exc_seq_region_end INTEGER,ori INTEGER,exc_type TEXT);
            INSERT INTO assembly_exception VALUES (20,1000,1099,10,100,199,1,'PAR'),(20,1500,1599,10,300,399,-1,'PAR');
            CREATE TABLE gene (stable_id TEXT,version INTEGER,gene_id INTEGER,is_current INTEGER,seq_region_id INTEGER);
            INSERT INTO gene VALUES ('G1',1,1,1,10),('G2',1,2,1,10),('G3',1,3,1,20);
        '''.split(';'):
            if sql.strip(): self.db.exec_driver_sql(sql)
        @contextmanager
        def connect(_): yield self.db
        self.patcher=patch('gxf_dump.exporter.sql_connection',connect);self.patcher.start();self.addCleanup(self.patcher.stop)
        self.args=argparse.Namespace(db_secrets='private-path',no_version=False,add_chr=False,format='gff3',
            start=None,end=None,vega=True,artifacts=False,include_hidden_remarks=False,source=None,analysis=None,
            gene_list=None,transcript_list=None,mane=None,biotype=None)
        self.regions=[dict(name='Y',seq_region_id=20)]

    def replay(self,args):
        self.assertIn('--vega',args)
        start=int(args[args.index('--start')+1])
        row=records('G1.1','T1.1','X',110,140) if start==100 else records('G2.1','T2.1','X',310,340,'-')
        write_features(row,args[args.index('--output')+1],args[args.index('--format')+1])

    def test_projection_ids_metadata_phases_and_gene_traversal(self):
        base=records('G3.1','T3.1','Y',200,240)
        result=add_y_par_features({},self.args,self.regions,base,self.replay)
        genes=[f for f in result if f.feature_type=='gene']
        self.assertEqual([f.attrs['ID'] for f in genes],['G1.1_PAR_Y','G2.1_PAR_Y','G3.1'])
        self.assertEqual([(f.start,f.end) for f in genes],[(1010,1040),(1559,1589),(200,240)])
        self.assertEqual([f.strand for f in genes],['+','+','+'])
        transcripts=[f for f in result if f.feature_type=='transcript']
        self.assertEqual(transcripts[0].attrs['Parent'],'G1.1_PAR_Y')
        self.assertEqual(transcripts[0].attrs['transcript_id'],'T1.1')
        exon=next(f for f in result if f.feature_type=='exon' and f.attrs['gene_id']=='G1.1')
        self.assertEqual(exon.attrs['ID'],'exon:T1.1_PAR_Y:1')
        self.assertEqual(exon.attrs['exon_id'],'ENSE1.1')
        self.assertEqual([f.phase for f in result if f.feature_type=='CDS'],['0','0','0'])
        self.assertEqual([f.attrs.get('tag') for f in genes],['PAR','PAR',None])
        self.assertEqual(base[0].attrs['ID'],'G3.1')

    def test_other_chromosomes_unchanged_and_no_duplicate_stored_copy(self):
        original=records('G1.1','T1.1','X',110,140)
        self.assertIs(add_y_par_features({},self.args,[dict(name='X',seq_region_id=10)],original,self.replay),original)
        stored=records('G1.1','T1.1','Y',1010,1040)
        result=add_y_par_features({},self.args,self.regions,stored,self.replay)
        self.assertEqual(sum(f.feature_type=='gene' for f in result),2)
        self.assertEqual(sum(f.attrs.get('ID')=='G1.1_PAR_Y' for f in result),1)

    def test_slice_prefix_and_invalid_mapping(self):
        self.args.start=1000;self.args.end=1100;self.args.add_chr=True
        result=add_y_par_features({},self.args,self.regions,[],self.replay)
        self.assertEqual({f.chrom for f in result},{'chrY'})
        self.assertEqual({f.attrs['gene_id'] for f in result},{'G1.1'})
        self.db.exec_driver_sql('UPDATE assembly_exception SET ori=0')
        with self.assertRaisesRegex(ValueError,'Invalid Y PAR'):
            add_y_par_features({},self.args,self.regions,[],self.replay)

    def test_gene_crossing_projection_boundary_is_not_partially_projected(self):
        def crossing(args): write_features(records('G1.1','T1.1','X',90,140),args[args.index('--output')+1],'gff3')
        from contextlib import redirect_stderr
        import io
        log=io.StringIO()
        with redirect_stderr(log):
            result=add_y_par_features({},self.args,self.regions,[],crossing)
        self.assertEqual(result,[])
        self.assertIn('Excluded non-projectable Y PAR gene',log.getvalue())


if __name__=='__main__': unittest.main()
