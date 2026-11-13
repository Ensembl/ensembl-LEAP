import csv
import gzip
import io
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / 'assets/check_stats_gff_consistency.py'
HEADER = ['gene_id', 'transcript_id', 'transcript_location', 'five_prime_extension',
          'five_prime_support', 'five_prime_read_count', 'three_prime_extension',
          'three_prime_support', 'three_prime_read_count']


def annotation(start=100, end=250, strand='+', chromosome='1', tid='t1.1', gene='g1', tags=None):
    def row(kind, a, b):
        attrs = 'gene_id={};transcript_id={}'.format(gene, tid)
        if kind == 'transcript' and tags is not None:
            attrs += ';leap_tag=' + tags
        return '{}\ttest\t{}\t{}\t{}\t.\t{}\t.\t{}\n'.format(chromosome, kind, a, b, strand, attrs)
    # Exons before transcript deliberately match pipeline output ordering.
    return row('exon', start, 150) + row('exon', 200, end) + row('transcript', start, end)


def fixture(strand='+', start=90, end=260, chromosome='1', tid='t1.1', gene='g1', blacklisted=False):
    five, three = (100 - start, end - 250) if strand == '+' else (end - 250, 100 - start)
    tags = [tag for length, tag in [(five, 'LEAP_5'), (three, 'LEAP_3')] if length > 0]
    if blacklisted:
        tags.append('LEAP_blacklisted')
    entry = dict(gene_id=gene, transcript_id=tid, transcript_location=f'{chromosome}:{start}-{end}',
                 five_prime_extension=str(max(0, five)), three_prime_extension=str(max(0, three)),
                 five_prime_support='reads' if five > 0 else 'blacklisted' if blacklisted else 'none',
                 five_prime_read_count='7' if five > 0 else 'blacklisted' if blacklisted else '0',
                 three_prime_support='current_annotation' if three > 0 else 'none',
                 three_prime_read_count='2' if three > 0 else '0')
    return (annotation(strand=strand, chromosome=chromosome, tid=tid, gene=gene),
            annotation(start, end, strand, chromosome, tid, gene, ','.join(tags) or 'LEAP_none'), entry)


def tsv(rows):
    handle = io.StringIO();writer = csv.DictWriter(handle, fieldnames=HEADER, delimiter='\t')
    writer.writeheader();writer.writerows(rows)
    return handle.getvalue()


def context(folder, strand='+', start=90, end=260, five_support='reads', three_support='current_annotation', five_count=7, three_count=2):
    fixed = folder / 'fixed_extensions';fixed.mkdir()
    for kind, side, support, count in [('FE', 'five', five_support, five_count), ('LE', 'three', three_support, three_count)]:
        left = (kind == 'FE') == (strand == '+')
        a, b = (start, 150) if left else (200, end)
        options = '.' if support == 'current_annotation' else str(100 if left else 250)
        (fixed / ('1_' + ('positive' if strand == '+' else 'negative') + '_' + kind + '_fixed.bed')).write_text(
            f'1\t{a}\t{b}\t{options}\t0\t{strand}\t.\t.\t{count}\tg1\n')


class StatsGffTests(unittest.TestCase):
    def run_check(self, raw, extended, entries, compressed=False, candidate_options=None, directory=False):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp);paths = []
            content = [('raw.gff3', raw), ('extended.gff3', extended), ('stats.tsv', tsv(entries) if isinstance(entries, list) else entries)]
            for name, data in content:
                path = root / (name + ('.gz' if compressed else ''))
                with (gzip.open if compressed else open)(path, 'wt') as handle:handle.write(data)
                paths.append(path)
            output = root / 'report.csv';stats = paths[2]
            if directory:
                stats = root / 'stats';stats.mkdir();paths[2].rename(stats / paths[2].name)
                (stats / 'notes.txt').write_text('ignored')
            command = [sys.executable, str(SCRIPT), str(paths[0]), str(paths[1]), str(output), '--extension-stats', str(stats)]
            if candidate_options is not None:
                options = dict(candidate_options);mutator = options.pop('mutate', None)
                context(root, **options)
                if mutator:mutator(root / 'fixed_extensions')
                command += ['--leap-output-dir', str(root)]
            result = subprocess.run(command, capture_output=True, text=True)
            rows = None
            if output.exists():
                with output.open(newline='') as handle:rows = list(csv.DictReader(handle))
            return result, rows

    def test_strand_mapping_unchanged_and_either_end_extensions(self):
        for strand in ('+', '-'):
            for start, end in [(100, 250), (90, 250), (100, 260), (90, 260)]:
                with self.subTest(strand=strand, start=start, end=end):
                    raw, extended, entry = fixture(strand, start, end)
                    result, rows = self.run_check(raw, extended, [entry], compressed=True, directory=True)
                    self.assertEqual(result.returncode, 0, result.stderr)
                    self.assertEqual(rows[0]['status'], 'consistent')
                    self.assertEqual(rows[0]['five_prime_expected_extension'], entry['five_prime_extension'])
                    self.assertEqual(rows[0]['three_prime_expected_extension'], entry['three_prime_extension'])
                    self.assertIn('PASS: 0 failing transcript key(s)', result.stdout)
                    self.assertIn('candidate_provenance: disabled', result.stdout)

    def test_annotation_counts_and_paired_mixed_lists(self):
        raw, extended, entry = fixture()
        for reads in ['0', '12']:
            entry['five_prime_support'] = 'current_annotation';entry['five_prime_read_count'] = reads
            self.assertEqual(self.run_check(raw, extended, [entry])[0].returncode, 0)
        entry['five_prime_support'] = 'reads,current_annotation';entry['five_prime_read_count'] = '3,0'
        self.assertEqual(self.run_check(raw, extended, [entry])[0].returncode, 0)
        entry['five_prime_support'] = 'reads,reads';entry['five_prime_read_count'] = '3,7'
        self.assertEqual(self.run_check(raw, extended, [entry])[0].returncode, 0)

    def test_blacklist_rules_both_strands_and_three_prime_extension(self):
        for strand in ('+', '-'):
            for other_extended in (False, True):
                start, end = (100, 260 if other_extended else 250) if strand == '+' else (90 if other_extended else 100, 250)
                raw, extended, entry = fixture(strand, start, end, blacklisted=True)
                result, rows = self.run_check(raw, extended, [entry])
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertIn('LEAP_blacklisted', rows[0]['expected_leap_tags'])
        raw, extended, entry = fixture(blacklisted=True)
        entry.update(five_prime_support='blacklisted', five_prime_read_count='blacklisted')
        result, rows = self.run_check(raw, extended, [entry])
        self.assertEqual(result.returncode, 1)
        self.assertEqual(rows[0]['support_mismatch'], 'True')
        raw, extended, entry = fixture(start=100, end=250)
        entry.update(three_prime_support='blacklisted', three_prime_read_count='blacklisted')
        self.assertEqual(self.run_check(raw, extended, [entry])[0].returncode, 1)

    def test_lengths_support_counts_and_tags_mismatches_independently(self):
        raw, extended, original = fixture()
        cases = [('five_prime_extension', '11', 'extension_length_mismatch'),
                 ('three_prime_extension', '0', 'extension_length_mismatch'),
                 ('five_prime_support', 'none', 'support_mismatch'),
                 ('five_prime_read_count', '0', 'read_count_mismatch'),
                 ('five_prime_read_count', 'blacklisted', 'read_count_mismatch'),
                 ('gene_id', 'g2', 'gene_mismatch'),
                 ('transcript_location', '1:91-260', 'location_mismatch')]
        for column, value, flag in cases:
            with self.subTest(column=column, value=value):
                entry = dict(original, **{column:value});result, rows = self.run_check(raw, extended, [entry])
                self.assertEqual(result.returncode, 1, result.stderr)
                self.assertEqual(rows[0][flag], 'True')
                self.assertTrue(any(d['issue']==flag for d in json.loads(rows[0]['details'])))
        for tag in ['', 'LEAP_5', 'LEAP_none,LEAP_5,LEAP_3', 'LEAP_5,LEAP_5,LEAP_3', 'LEAP_unknown', 'LEAP_blacklisted,LEAP_3']:
            bad = extended.replace('leap_tag=LEAP_5,LEAP_3', 'leap_tag=' + tag)
            result, rows = self.run_check(raw, bad, [original])
            self.assertEqual(result.returncode, 1, result.stderr)
            self.assertEqual(rows[0]['leap_tag_mismatch'], 'True')
        reversed_tags = extended.replace('LEAP_5,LEAP_3', 'LEAP_3,LEAP_5')
        self.assertEqual(self.run_check(raw, reversed_tags, [original])[0].returncode, 0)
        self.assertEqual(self.run_check(raw, extended.replace(';leap_tag=LEAP_5,LEAP_3', ''), [original])[0].returncode, 1)

    def test_zero_extensions_and_restored_ends_clear_support_and_counts(self):
        raw, extended, entry = fixture(start=100, end=250)
        for changes, flag in [({'five_prime_support':'reads','five_prime_read_count':'2'}, 'support_mismatch'),
                              ({'three_prime_read_count':'2'}, 'read_count_mismatch'),
                              ({'five_prime_support':'blacklisted','five_prime_read_count':'0'}, 'read_count_mismatch')]:
            result, rows = self.run_check(raw, extended, [dict(entry, **changes)])
            self.assertEqual(result.returncode, 1)
            self.assertEqual(rows[0][flag], 'True')
        # Processed output uses corrected stats; original saved read endpoints need not be present.
        result, rows = self.run_check(raw, extended, [entry], candidate_options={})
        self.assertEqual(result.returncode, 0, result.stderr)
        bad = extended.replace('LEAP_none', 'LEAP_blacklisted')
        result, rows = self.run_check(raw, bad, [entry])
        self.assertEqual(result.returncode, 1)
        self.assertEqual(rows[0]['leap_tag_mismatch'], 'True')

    def test_candidate_endpoint_support_counts_gene_and_last_row_selection(self):
        for strand in ('+', '-'):
            raw, extended, entry = fixture(strand)
            result, rows = self.run_check(raw, extended, [entry], candidate_options={'strand':strand})
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn('candidate_provenance: enabled', result.stdout)
            for changes in [{'five_count':8}, {'five_support':'current_annotation'}, {'start':89}]:
                options = dict(changes, strand=strand)
                result, rows = self.run_check(raw, extended, [entry], candidate_options=options)
                self.assertEqual(result.returncode, 1, result.stderr)
                self.assertEqual(rows[0]['candidate_mismatch'], 'True')
        def prefix_annotation(files):
            f=next(files.glob('*FE*'));f.write_text('1\t90\t150\t.\t0\t+\t.\t.\t9\tg1\n'+f.read_text())
        raw, extended, entry = fixture()
        self.assertEqual(self.run_check(raw, extended, [entry], candidate_options={'mutate':prefix_annotation})[0].returncode, 0)
        self.assertEqual(self.run_check(raw, extended, [dict(entry, five_prime_read_count='0007')], candidate_options={})[0].returncode, 0)
        def wrong_gene(files):
            f=next(files.glob('*FE*'));f.write_text(f.read_text().replace('g1', 'g2'))
        result, rows = self.run_check(raw, extended, [entry], candidate_options={'mutate':wrong_gene})
        self.assertEqual(result.returncode, 1)
        self.assertEqual(rows[0]['candidate_mismatch'], 'True')
        def single_exon_context(files):
            for f in files.glob('*.bed'):
                with f.open('a') as handle:
                    handle.write('1\t500\t600\tsingle_exon\t0\t+\tg2\t.\t.\t0\n')
        self.assertEqual(self.run_check(raw, extended, [entry], candidate_options={'mutate':single_exon_context})[0].returncode, 0)
        raw, extended, entry = fixture(start=100, end=250)
        raw = raw.replace('exon\t100\t150', 'exon\t100\t250')
        extended = extended.replace('exon\t100\t150', 'exon\t100\t250')
        raw = ''.join(line+'\n' for line in raw.splitlines() if 'exon\t200\t250' not in line)
        extended = ''.join(line+'\n' for line in extended.splitlines() if 'exon\t200\t250' not in line)
        self.assertEqual(self.run_check(raw, extended, [entry], candidate_options={'mutate':single_exon_context})[0].returncode, 0)

    def test_membership_exact_versions_and_par_copies(self):
        raw, extended, entry = fixture();raw2, extended2, entry2 = fixture(tid='t2.1')
        result, rows = self.run_check(raw+raw2, extended, [entry])
        self.assertEqual(result.returncode, 0)
        self.assertEqual(rows[1]['status'], 'missing_from_extended')
        result, rows = self.run_check(raw+raw2, extended+extended2, [entry])
        self.assertEqual(result.returncode, 1)
        self.assertEqual(rows[1]['missing_stats'], 'True')
        result, rows = self.run_check(raw+raw2, extended, [entry,entry2])
        self.assertEqual(result.returncode, 1)
        self.assertEqual(rows[1]['extra_stats'], 'True')
        raw2, extended2, entry2 = fixture(tid='t1.2')
        result, rows = self.run_check(raw, extended2, [entry2])
        self.assertEqual(result.returncode, 1)
        self.assertEqual(rows[1]['missing_from_raw'], 'True')
        x, xe, xs = fixture(chromosome='X');y, ye, ys = fixture(chromosome='Y', strand='-')
        result, rows = self.run_check(x+y, xe+ye, [xs,ys])
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual([r['chromosome'] for r in rows], ['X','Y'])

    def test_geometry_strand_and_shortening(self):
        raw, extended, entry = fixture()
        bad = extended.replace('transcript\t90\t260', 'transcript\t91\t260')
        result, rows = self.run_check(raw, bad, [entry])
        self.assertEqual(result.returncode, 1)
        self.assertEqual(rows[0]['exon_span_mismatch'], 'True')
        self.assertEqual(rows[0]['location_mismatch'], 'True')
        result, rows = self.run_check(raw, extended.replace('\t+\t', '\t-\t'), [entry])
        self.assertEqual(result.returncode, 1)
        self.assertEqual(rows[0]['location_mismatch'], 'True')
        raw, extended, entry = fixture(start=101, end=260)
        result, rows = self.run_check(raw, extended, [entry])
        self.assertEqual(result.returncode, 1)
        self.assertEqual(rows[0]['shortened'], 'True')
        # Unmatched exonless raw models are irrelevant to a stats/coordinate comparison.
        exonless = '1\ttest\ttranscript\t500\t600\t.\t+\t.\tgene_id=g2;transcript_id=filtered\n'
        raw, extended, entry = fixture()
        self.assertEqual(self.run_check(raw+exonless, extended, [entry])[0].returncode, 0)

    def test_invalid_stats_and_annotation_inputs(self):
        raw, extended, entry = fixture();valid=tsv([entry])
        invalid_stats = [('', 'columns'), (tsv([]), 'no stats rows'),
                         (tsv([entry,entry]), 'duplicate stats'),
                         (valid.replace('gene_id\t', 'transcript_id\t',1), 'columns'),
                         (tsv([dict(entry, five_prime_extension='-1')]), 'nonnegative integer'),
                         (tsv([dict(entry, five_prime_read_count='-1')]), 'nonnegative integer'),
                         (tsv([dict(entry, five_prime_read_count='x')]), 'nonnegative integer'),
                         (tsv([dict(entry, five_prime_support='other')]), 'unknown support'),
                         (tsv([dict(entry, five_prime_support='reads,current_annotation')]), 'different lengths'),
                         (tsv([dict(entry, transcript_location='bad')]), 'transcript_location'),
                         (tsv([dict(entry, transcript_id='')]), 'incomplete'),
                         (valid.rstrip()+'\textra\n', 'extra stats fields')]
        for invalid, message in invalid_stats:
            result, rows = self.run_check(raw, extended, invalid)
            self.assertEqual(result.returncode, 2, result.stderr)
            self.assertIn(message, result.stderr)
            self.assertIsNone(rows)
        lines=extended.splitlines(keepends=True)
        for invalid, message in [('', 'no transcript'), (extended+lines[-1], 'duplicate transcript'),
                                 (''.join(lines[:2]), 'no transcript'), (lines[-1], 'no exons'),
                                 (extended+lines[0], 'duplicate exons'),
                                 (extended.replace('gene_id=g1', 'broken'), 'malformed'),
                                 (extended.replace('transcript_id=t1.1','missing=t1.1'), 'missing transcript_id'),
                                 (extended.replace('\t90\t150\t','\t0\t150\t'), 'invalid coordinates')]:
            result, rows=self.run_check(raw, invalid, [entry]);self.assertEqual(result.returncode,2,result.stderr)
            self.assertIn(message,result.stderr);self.assertIsNone(rows)

    def test_context_errors_and_output_input_protection(self):
        raw, extended, entry = fixture()
        for mutator in [lambda d: next(d.glob('*FE*')).unlink(),
                        lambda d: (d/'copy_1_positive_FE_fixed.bed').write_text(next(d.glob('*FE*')).read_text()),
                        lambda d: next(d.glob('*FE*')).write_text('invalid\n')]:
            result, rows = self.run_check(raw, extended, [entry], candidate_options={'mutate':mutator})
            self.assertEqual(result.returncode, 2, result.stderr);self.assertIsNone(rows)
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);r=root/'raw.gff3';e=root/'extended.gff3';s=root/'stats.tsv'
            r.write_text(raw);e.write_text(extended);s.write_text(tsv([entry]));context(root)
            alias=root/'alias.tsv';alias.symlink_to(s)
            hard=root/'hard.tsv';hard.hardlink_to(s)
            for output in [r,e,s,alias,hard,next((root/'fixed_extensions').glob('*FE*')),root/'absent/report.csv']:
                command=[sys.executable,str(SCRIPT),str(r),str(e),str(output),'--extension-stats',str(s),'--leap-output-dir',str(root)]
                result=subprocess.run(command,capture_output=True,text=True)
                self.assertEqual(result.returncode,2,result.stderr)
                self.assertEqual(r.read_text(),raw);self.assertEqual(e.read_text(),extended);self.assertEqual(s.read_bytes(),tsv([entry]).encode())
            report=root/'report.csv';report.write_text('previous report')
            result=subprocess.run([sys.executable,str(SCRIPT),str(r),str(e),str(report),'--extension-stats',str(root/'missing')],capture_output=True,text=True)
            self.assertEqual(result.returncode,2);self.assertEqual(report.read_text(),'previous report')


if __name__ == '__main__':
    unittest.main()
