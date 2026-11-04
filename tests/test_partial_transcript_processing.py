"""Exercise the production scripts with partial models retained until output."""
import csv
import gzip
import importlib.util
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
BIN = ROOT / 'bin'
PANDAS = importlib.util.find_spec('pandas') is not None


def model(tid, gene='A', strand='+', exons=((200, 250), (400, 450)),
          cds=((220, 240), (400, 420)), five=True, three=True,
          biotype='protein_coding', chromosome='1'):
    attrs = f'gene_id={gene};transcript_id={tid};transcript_type={biotype};gene_type=protein_coding;gene_name={gene}'
    def record(kind, start, end, extra=''):
        return f'{chromosome}\ttest\t{kind}\t{start}\t{end}\t.\t{strand}\t.\t{attrs}{extra}\n'
    rows = [record('transcript', min(a for a, b in exons), max(b for a, b in exons))]
    ordered = sorted(exons, reverse=strand == '-')
    for number, (start, end) in enumerate(ordered, 1):
        extra = f';exon_id={tid}_E{number};exon_number={number}'
        rows.append(record('exon', start, end, extra))
        for a, b in cds:
            if start <= a <= b <= end:
                rows.append(record('CDS', a, b, extra))
        for kind, enabled, interval in [
            ('five_prime_UTR', five, (200, 219) if strand == '+' else (421, 450)),
            ('three_prime_UTR', three, (421, 450) if strand == '+' else (200, 219)),
        ]:
            a, b = interval
            if enabled and start <= a <= b <= end:
                rows.append(record(kind, a, b, extra))
    return ''.join(rows)


def parse(path):
    rows = []
    for line in path.read_text().splitlines():
        if not line or line.startswith('#'):
            continue
        fields = line.split('\t')
        attrs = dict(x.split('=', 1) for x in fields[8].split(';') if '=' in x)
        rows.append((fields, attrs))
    return rows


@unittest.skipUnless(PANDAS, 'production annotation classification requires pandas')
class PartialProcessingTests(unittest.TestCase):
    def run_script(self, folder, script, *args, ok=True):
        env = dict(os.environ, PYTHONHASHSEED='0', PYTHONDONTWRITEBYTECODE='1')
        result = subprocess.run([sys.executable, str(BIN / script), *map(str, args)],
                                cwd=folder, capture_output=True, text=True, env=env)
        if ok:
            self.assertEqual(result.returncode, 0, result.stderr)
        return result

    def classify(self, folder, text):
        source = folder / 'raw.gff3'; source.write_text(text)
        self.run_script(folder, 'filter_out_partials.py', source, folder / 'annotation')
        return {suffix: folder / f'annotation_{suffix}.gff' for suffix in ['FT', 'PT', 'PC', 'NPC']}

    def candidates(self, folder, strand='+', five_start=150, five_end=250, five_options='185/200'):
        sign = 'positive' if strand == '+' else 'negative'
        fe = folder / f'1_{sign}_FE.bed'; le = folder / f'1_{sign}_LE.bed'
        if strand == '+':
            fe.write_text(f'1\t{five_start}\t{five_end}\t{five_options}\t0\t+\t.\t.\t3\tA\n')
            le.write_text('1\t400\t500\t450\t0\t+\t.\t.\t4\tA\n')
        else:
            fe.write_text('1\t400\t500\t470/450\t0\t-\t.\t.\t3\tA\n')
            le.write_text('1\t150\t250\t200\t0\t-\t.\t.\t4\tA\n')
        return fe, le

    def build(self, folder, source, partial, fe=None, le=None, output='extended'):
        if fe is None:
            fe, le = self.candidates(folder)
        blacklist = folder / 'blacklist.txt'; blacklist.touch()
        gff, stats = folder / f'{output}.gff', folder / f'{output}.tsv'
        self.run_script(folder, 'build_chromosome_gff_and_stats.py', fe, le, source,
                        blacklist, gff, stats, '--partial-annotation', partial)
        return gff, stats

    def test_classification_keeps_all_coding_features_and_existing_outputs(self):
        for five, three in [(False, True), (True, False), (False, False)]:
            with self.subTest(five=five, three=three), tempfile.TemporaryDirectory() as tmp:
                p = Path(tmp)
                full = model('full.1'); partial = model('partial.1', five=five, three=three)
                npc = model('noncoding.1', biotype='lncRNA')
                files = self.classify(p, full + partial + npc)
                self.assertEqual(files['FT'].read_text(), full)
                self.assertEqual(files['PT'].read_text(), partial)
                self.assertEqual(files['NPC'].read_text(), npc)
                self.assertEqual(files['PC'].read_text(), full + partial)

    def test_empty_coding_or_partial_input_is_valid(self):
        for text in ['', model('npc.1', biotype='lncRNA'), model('full.1')]:
            with self.subTest(text=text[:40]), tempfile.TemporaryDirectory() as tmp:
                files = self.classify(Path(tmp), text)
                self.assertEqual(files['PT'].read_text(), '')
                self.assertEqual(files['PC'].read_text(), text if 'protein_coding;' in text.split('transcript_type=')[-1].split('gene_type=')[0] else '')

    def test_existing_cds_filter_shortens_candidates_with_partial_context(self):
        for strand in ['+', '-']:
            for target_gene in ['A', 'B']:
                with self.subTest(strand=strand, gene=target_gene), tempfile.TemporaryDirectory() as tmp:
                    p = Path(tmp)
                    blocker = model('partial.1', target_gene, strand,
                                    exons=((100, 170), (300, 350)) if strand == '+' else ((470, 500), (550, 600)),
                                    cds=((160, 170), (300, 330)) if strand == '+' else ((480, 490), (560, 570)),
                                    five=False, three=False)
                    files = self.classify(p, model('source.1', strand=strand) + blocker)
                    fe, le = self.candidates(p, strand)
                    single, npc, natural = p/'single.gff', p/'npc.gff', p/'natural.txt'
                    single.touch(); npc.touch(); natural.touch()
                    outputs = {}
                    for label in ['FT', 'PC']:
                        prefix = p / f'{label}_1_{"positive" if strand == "+" else "negative"}_FE'
                        self.run_script(p, 'stop_cds_overlap.py', files[label], fe,
                                        prefix, single, npc, natural)
                        outputs[label] = Path(str(prefix) + 'C.bed')
                    before = outputs['FT'].read_text().split('\t')
                    after = outputs['PC'].read_text().split('\t')
                    boundary = 1 if strand == '+' else 2
                    self.assertEqual(int(before[boundary]), 150 if strand == '+' else 500)
                    self.assertEqual(int(after[boundary]), 185 if strand == '+' else 470)
                    gff, stats = self.build(p, files['PC'], files['PT'], outputs['PC'], le)
                    self.assertNotIn('partial.1', gff.read_text())
                    rows = list(csv.DictReader(stats.read_text().splitlines(), delimiter='\t'))
                    self.assertEqual([r['transcript_id'] for r in rows], ['source.1'])
                    self.assertEqual(rows[0]['five_prime_extension'], '15' if strand == '+' else '20')
                    self.assertEqual(rows[0]['five_prime_support'], 'reads')

    def test_partial_unique_exons_enter_existing_flank_logic(self):
        spec = importlib.util.spec_from_file_location('production_cds', BIN/'stop_cds_overlap.py')
        module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
        for strand in ['+', '-']:
            with self.subTest(strand=strand), tempfile.TemporaryDirectory() as tmp:
                p = Path(tmp)
                # One unique CDS segment per gene; partial contributes distinct exons on both sides.
                full = model('full.1', strand=strand, exons=((200,250),(300,350),(400,450)), cds=((310,320),))
                partial = model('partial.1', strand=strand, exons=((100,150),(300,350),(500,550)), cds=((310,320),), five=False, three=False)
                files = self.classify(p, full + partial)
                _, special, exons = module.split_gff(str(files['PC']), None, {'A'})
                self.assertEqual(special, ['A'])
                left = module.filter_exons_by_side(exons, 'left', 'positive' if strand=='+' else 'negative')
                right = module.filter_exons_by_side(exons, 'right', 'positive' if strand=='+' else 'negative')
                for frame in [left, right]:
                    self.assertIn('partial.1', set(frame['transcript_id']))
                coords = {(int(r.start), int(r.end)) for frame in [left,right] for r in frame.itertuples()}
                self.assertTrue({(100,150),(500,550)} <= coords)

    def test_distinct_partial_cds_updates_existing_single_cds_gene_classification(self):
        spec = importlib.util.spec_from_file_location('production_cds', BIN/'stop_cds_overlap.py')
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        for strand in ['+', '-']:
            with self.subTest(strand=strand), tempfile.TemporaryDirectory() as tmp:
                folder = Path(tmp)
                full = model('full.1', strand=strand, cds=((220, 240),))
                partial = model('partial.1', strand=strand, cds=((225, 240),),
                                five=False, three=False)
                files = self.classify(folder, full + partial)
                original_cds, original_special, _ = module.split_gff(
                    str(files['FT']), None, {'A'})
                combined_cds, combined_special, _ = module.split_gff(
                    str(files['PC']), None, {'A'})
                self.assertEqual(len(original_cds), 1)
                self.assertEqual(original_special, ['A'])
                self.assertEqual(len(combined_cds), 2)
                self.assertEqual(combined_special, [])

    def test_partial_features_stats_and_gene_extents_do_not_leak(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)
            full = model('source.1')
            shared = model('source.2', exons=((1,100),(900,1000)), cds=((20,80),), five=False, three=False)
            only_partial = model('partial_only.1', 'P', five=False)
            files = self.classify(p, full + shared + only_partial)
            gff, stats = self.build(p, files['PC'], files['PT'])
            records = parse(gff)
            self.assertEqual({a['transcript_id'] for f,a in records if f[2]!='gene'}, {'source.1'})
            genes = [(f,a) for f,a in records if f[2]=='gene']
            self.assertEqual([(a['gene_id'], f[3], f[4]) for f,a in genes], [('A','150','500')])
            self.assertEqual([r['transcript_id'] for r in csv.DictReader(stats.read_text().splitlines(),delimiter='\t')], ['source.1'])
            all_partial = p/'all_partial.gff'; all_partial.write_text(shared + only_partial)
            empty_gff, empty_stats = self.build(p, all_partial, files['PT'], output='empty')
            self.assertEqual(empty_gff.read_text(), '')
            self.assertEqual(list(csv.DictReader(empty_stats.read_text().splitlines(),delimiter='\t')), [])

    def test_empty_partial_builder_retains_full_models_and_both_end_stats(self):
        for strand in ['+', '-']:
            with self.subTest(strand=strand), tempfile.TemporaryDirectory() as tmp:
                p=Path(tmp);source=p/'source.gff';source.write_text(model('source.1',strand=strand))
                partial=p/'empty.gff';partial.touch();fe,le=self.candidates(p,strand)
                gff,stats=self.build(p,source,partial,fe,le)
                transcript=next(f for f,a in parse(gff) if f[2]=='transcript')
                self.assertEqual(transcript[3:5],['150','500'])
                original_cds={tuple(f[3:5]) for f,a in parse(source) if f[2]=='CDS'}
                self.assertEqual({tuple(f[3:5]) for f,a in parse(gff) if f[2]=='CDS'},original_cds)
                row=list(csv.DictReader(stats.read_text().splitlines(),delimiter='\t'))[0]
                self.assertEqual((row['five_prime_extension'],row['three_prime_extension']),('50','50'))
                self.assertEqual((row['five_prime_support'],row['three_prime_support']),('reads','reads'))

    def test_exclusion_uses_exact_version_chromosome_and_strand_and_accepts_gzip(self):
        for other in [model('source.2',five=False), model('source.1',strand='-',five=False),
                      model('source.1',chromosome='2',five=False)]:
            with self.subTest(other=other[:100]), tempfile.TemporaryDirectory() as tmp:
                p = Path(tmp); source = p/'source.gff'; source.write_text(model('source.1'))
                partial = p/'partial.gff.gz'
                with gzip.open(partial,'wt') as handle:handle.write(other)
                gff, _ = self.build(p, source, partial)
                self.assertIn('transcript_id=source.1;', gff.read_text())

    def test_missing_or_malformed_partial_input_stops_before_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp); source = p/'source.gff'; source.write_text(model('source.1'))
            fe,le = self.candidates(p); blacklist=p/'blacklist';blacklist.touch()
            for text in ['bad\n','1\ttest\texon\t1\t2\t.\t+\t.\tgene_id=A\n']:
                partial=p/'bad.gff';partial.write_text(text)
                result=self.run_script(p,'build_chromosome_gff_and_stats.py',fe,le,source,blacklist,p/'out.gff',p/'out.tsv','--partial-annotation',partial,ok=False)
                self.assertNotEqual(result.returncode,0)
                self.assertIn('bad.gff:1:',result.stderr)
                self.assertFalse((p/'out.gff').exists())
            result=self.run_script(p,'build_chromosome_gff_and_stats.py',fe,le,source,blacklist,p/'out.gff',p/'out.tsv',ok=False)
            self.assertEqual(result.returncode,2)

    @unittest.skipUnless(shutil.which('bedtools'), 'requires production bedtools executable')
    def test_partial_terminals_and_single_exons_reach_existing_overlap_processing(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)
            partial=model('partial.1','B',exons=((140,180),(300,350)),cds=((150,170),),five=False,three=False)
            single=model('single.1','A',exons=((205,225),),cds=((210,220),),five=False,three=False)
            files=self.classify(p,model('full.1')+partial+single)
            self.run_script(p,'extract_single_exons.py',files['PC'],'1',p/'annotation_PC')
            multi=p/'annotation_PC_M.gff3';self.assertIn('partial.1',multi.read_text())
            self.assertIn('single.1',(p/'annotation_PC_S.gff3').read_text())
            self.assertIn('single_exon',(p/'1_positive_SE.bed').read_text())
            self.run_script(p,'extract_first_and_last_exon.py',multi,p/'annotation_1_positive')
            self.assertIn('partial.1',(p/'annotation_1_positive_FE.gff').read_text())
            self.assertIn('partial.1',(p/'annotation_1_positive_LE.gff').read_text())

    @unittest.skipUnless(shutil.which('bedtools'), 'requires production bedtools executable')
    def test_partial_single_exon_without_cds_blocks_an_unknown_terminal_contact(self):
        for strand in ['+', '-']:
            with self.subTest(strand=strand), tempfile.TemporaryDirectory() as tmp:
                p=Path(tmp)
                partial=model('partial.1','B',strand,exons=((160,180),),cds=(),five=False,three=False)
                files=self.classify(p,model('full.1',strand=strand)+partial)
                self.run_script(p,'extract_single_exons.py',files['PC'],'1',p/'annotation_PC')
                fe,le=self.candidates(p,strand);sign='positive' if strand=='+' else 'negative'
                single=p/f'1_{sign}_SE.bed';overlap=p/'contacts.bed'
                cmd=['bedtools','intersect','-a',str(fe if strand=='+' else single),'-b',str(single if strand=='+' else le),'-wa','-wb']
                result=subprocess.run(cmd,capture_output=True,text=True,check=True)
                self.assertTrue(result.stdout)
                overlap.write_text(result.stdout)
                natural=p/'natural.txt';natural.touch();(p/'safe_FE.bed').touch();(p/'safe_LE.bed').touch()
                self.run_script(p,'overlap_safety_checker.py',natural,overlap,p/'safe_FE',p/'safe_LE',p/'unsafe.bed')
                self.assertTrue((p/'unsafe.bed').read_text())
                self.assertEqual((p/'safe_FEC.bed').read_text(),'')
                self.assertEqual((p/'safe_LEC.bed').read_text(),'')

    @unittest.skipUnless(shutil.which('bedtools'), 'requires production bedtools executable')
    def test_partial_terminal_features_contribute_recorded_natural_pairs(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)
            partial=model('partial.1','B',exons=((50,100),(180,220)),cds=((60,90),),five=False,three=False)
            files=self.classify(p,model('full.1')+partial)
            self.run_script(p,'extract_first_and_last_exon.py',files['PC'],p/'annotation_1_positive')
            self.run_script(p,'find_natural_overlaps.py',p/'annotation_1_positive_FE.gff',p/'annotation_1_positive_LE.gff',p/'natural')
            natural=p/'natural_NO.txt';self.assertIn('A\tB\n',natural.read_text())
            fe,_=self.candidates(p)
            target=p/'partial_LE.bed';target.write_text('1\t180\t220\t.\t0\t+\t.\t.\t0\tB\n')
            contact=subprocess.check_output(['bedtools','intersect','-a',str(fe),'-b',str(target),'-wa','-wb'],text=True)
            overlap=p/'contacts.bed';overlap.write_text(contact)
            (p/'safe_FE.bed').touch();(p/'safe_LE.bed').touch()
            self.run_script(p,'overlap_safety_checker.py',natural,overlap,p/'safe_FE',p/'safe_LE',p/'unsafe.bed')
            self.assertEqual((p/'unsafe.bed').read_text(),'')
            self.assertTrue((p/'safe_FEC.bed').read_text())

    def test_partial_annotation_participates_in_shared_splice_site_unification(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)
            partial=model('partial.1',exons=((100,250),(400,450)),five=False)
            files=self.classify(p,model('full.1')+partial)
            self.run_script(p,'extract_first_and_last_exon.py',files['PC'],p/'annotation_1_positive')
            self.run_script(p,'extract_unique_splice_sites_from_current_annotation.py',p/'annotation_1_positive_FE.gff',p/'annotation_1_positive_FE')
            rows=list(csv.DictReader((p/'annotation_1_positive_FE_USS.csv').read_text().splitlines()))
            self.assertEqual([(r['splice_id'],r['non_splice_site']) for r in rows],[('1_250_+','100')])
            fe,le=self.candidates(p)
            fe.write_text('1\t100\t250\t.\t0\t+\t.\t.\t0\tA\n')
            gff,stats=self.build(p,files['PC'],files['PT'],fe,le)
            self.assertNotIn('partial.1',gff.read_text())
            row=list(csv.DictReader(stats.read_text().splitlines(),delimiter='\t'))[0]
            self.assertEqual(row['five_prime_extension'],'100')
            self.assertEqual(row['five_prime_support'],'current_annotation')


if __name__ == '__main__':
    unittest.main()
