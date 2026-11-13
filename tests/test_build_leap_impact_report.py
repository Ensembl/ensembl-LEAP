import csv
import gzip
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

from datacheck_fixtures import write_run

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / 'assets/build_leap_impact_report.py'
sys.path.insert(0, str(ROOT / 'assets'))
spec = importlib.util.spec_from_file_location('leap_impact', SCRIPT)
impact = importlib.util.module_from_spec(spec); spec.loader.exec_module(impact)


class ImpactTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(); self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.run = write_run(self.root / 'run')

    def invoke(self, expected=0, extra=()):
        output = self.root / ('out' + str(len(list(self.root.glob('out*')))))
        sources = {p: p.read_bytes() for p in self.run.rglob('*') if p.is_file()}
        result = subprocess.run([sys.executable, str(SCRIPT), '--leap-output-dir', str(self.run),
            '--output', str(output), '--chromosomes', '1', *extra], capture_output=True, text=True)
        self.assertEqual(result.returncode, expected, result.stdout + result.stderr)
        self.assertNotIn('Traceback', result.stderr)
        self.assertEqual({p: p.read_bytes() for p in sources}, sources)
        return output, json.loads((output / 'summary.json').read_text())

    def edit_stats(self, change):
        path = self.run / 'extended_transcripts_stats/1_positive_stats.tsv'
        with path.open() as handle:
            reader = csv.DictReader(handle, delimiter='\t'); fields = reader.fieldnames; rows = list(reader)
        change(rows)
        with path.open('w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t'); writer.writeheader(); writer.writerows(rows)

    def test_both_strands_outputs_counts_lengths_and_offline_html(self):
        output, summary = self.invoke()
        totals = summary['totals']['leap']
        self.assertEqual((totals['transcripts'], totals['extended_genes'], totals['added_bases']), (2, 2, 40))
        self.assertEqual(totals['novel_read_supported_bases'], 40)
        self.assertEqual(totals['newly_covered_bases'], 40)
        self.assertEqual(totals['newly_exonic_bases'], 40)
        self.assertEqual(summary['totals']['external']['added_bases'], 40)
        self.assertEqual(summary['datachecks']['status'], 'NOT_RUN')
        self.assertEqual(summary['length_statistics'][0]['median'], 10)
        with (output / 'transcripts.csv').open() as handle:
            rows = list(csv.DictReader(handle))
        negative = next(r for r in rows if r['strand'] == 'negative')
        self.assertEqual(negative['jbrowse_label'], '1:490-760')
        self.assertEqual(negative['five_prime_external_extension'], '10')
        page = (output / 'index.html').read_text()
        self.assertNotIn('__LEAP_IMPACT_DATA__', page)
        self.assertNotIn('src="http', page)
        self.assertIn('New read-supported genomic bases', page)
        self.assertIn('validationStatus', page)
        self.assertIn('2 / 2 models extended', (output / 'summary.txt').read_text())
        with (output / 'overall_statistics.csv').open() as handle:
            self.assertEqual(len(list(csv.DictReader(handle))), 3)

    def test_duplicate_models_do_not_inflate_unique_sequence_or_gene_counts(self):
        for path in (self.run / 'raw_annotation/current_annotation_1.gff3', self.run / 'extended_transcripts/1_positive_extended.gff'):
            text = path.read_text()
            copied = ''.join(line.replace('t_plus.1', 't_copy.1') for line in text.splitlines(True)
                             if 'transcript_id=t_plus.1;' in line)
            path.write_text(text + copied)
        self.edit_stats(lambda rows: rows.append(dict(rows[0], transcript_id='t_copy.1')))
        _, summary = self.invoke()
        totals = summary['totals']['leap']
        self.assertEqual(totals['added_bases'], 60)
        self.assertEqual(totals['unique_added_bases'], 40)
        self.assertEqual(totals['novel_read_supported_bases'], 40)
        self.assertEqual(totals['extended_genes'], 2)
        self.assertEqual(totals['newly_covered_bases'], 40)

    def test_all_feature_coverage_includes_gene_and_transcript_spans_and_unknown_types(self):
        raw = self.run / 'raw_annotation/current_annotation_1.gff3'
        extra = '1\ttest\tgene\t90\t95\t.\t-\t.\tID=opposite_gene;gene_id=opposite_gene\n'
        extra += '1\ttest\tcustom_feature\t98\t99\t.\t.\t.\tID=unknown\n'
        extra += '1\ttest\tregulatory_region\t352\t353\t.\t+\t.\tID=stranded_unknown\n'
        raw.write_text(raw.read_text() + extra)
        _, summary = self.invoke()
        totals = summary['totals']['leap']
        self.assertEqual(totals['newly_covered_bases'], 30)
        self.assertEqual(totals['newly_exonic_bases'], 40)
        self.assertEqual(totals['new_same_strand_feature_bases'], 36)

    def test_genomic_bases_merge_opposite_strands_and_stages_independently(self):
        coverage = {'1': dict(all_features=[(100, 200)], all_exons=[(100, 110)],
                             positive_features=[(100, 200)], negative_features=[])}
        def row(strand, side, start, end):
            result = dict(chromosome='1', strand=strand, raw_start=start, raw_end=end,
                          dropped_no_cds=False, external_start=start, external_end=end)
            for s in impact.SIDES:
                result[s+'_extension'] = 10 if s == side else 0
                result[s+'_external_extension'] = 10 if s == side else 0
                result[s+'_internal_extension'] = 0
                result[s+'_support'] = 'reads' if s == side else 'none'
                result[s+'_read_count'] = '1' if s == side else '0'
            return result
        a = row('positive', 'five_prime', 100, 200)
        b = row('negative', 'three_prime', 100, 200)
        records, intervals = impact.genomic_coverage([a, b, dict(a)], coverage)
        record = next(r for r in records if r['stage']=='leap' and r['strand']=='both')
        self.assertEqual(record['unique_genomic_extension_bases'], 10)
        self.assertEqual(record['newly_covered_bases'], 10)
        self.assertEqual(record['new_read_supported_bases'], 10)
        self.assertEqual(record['unique_genomic_read_supported_bases'], 10)
        self.assertEqual(sum(r['bases'] for r in intervals if r['stage']=='leap'), 10)
        b['three_prime_external_extension'] = 0
        b['three_prime_internal_extension'] = 10
        records, _ = impact.genomic_coverage([a, b], coverage)
        self.assertEqual([r['newly_covered_bases'] for r in records if r['strand']=='both'], [10, 10, 10])
        coverage['1']['all_features'] = [(80, 210)]
        records, _ = impact.genomic_coverage([a, b], coverage)
        record = next(r for r in records if r['stage']=='leap' and r['strand']=='both')
        self.assertEqual(record['newly_covered_bases'], 0)
        self.assertEqual(record['newly_exonic_bases'], 10)
        coverage['1']['all_features'] = [(99, 200)]
        records, _ = impact.genomic_coverage([a, b], coverage)
        record = next(r for r in records if r['stage']=='leap' and r['strand']=='both')
        self.assertEqual(record['newly_covered_bases'], 9)
        coverage['2'] = coverage['1']
        records, _ = impact.genomic_coverage([a, dict(a, chromosome='2')], coverage)
        self.assertEqual(sum(r['newly_covered_bases'] for r in records if r['stage']=='leap' and r['strand']=='both'), 18)

    def test_novel_sequence_subtracts_partial_and_noncoding_raw_exons(self):
        path = self.run / 'raw_annotation/current_annotation_1.gff3'
        extra = ''
        for kind in ('transcript', 'exon'):
            extra += '1\ttest\t{}\t95\t100\t.\t+\t.\tID={};Parent={};transcript_id=partial.1;gene_id=g_plus;transcript_type=lncRNA\n'.format(
                kind, 'partial.1' if kind == 'transcript' else 'partial_exon', 'g_plus' if kind == 'transcript' else 'partial.1')
        path.write_text(path.read_text() + extra)
        _, summary = self.invoke()
        self.assertEqual(summary['totals']['leap']['unique_read_supported_bases'], 40)
        self.assertEqual(summary['totals']['leap']['novel_read_supported_bases'], 35)
        self.assertEqual(summary['totals']['leap']['eligible_raw_transcripts'], 2)

    def test_opposite_ends_can_use_different_sources_and_repeated_read_counts_are_not_summed(self):
        self.edit_stats(lambda rows: rows[0].update(five_prime_support='current_annotation', five_prime_read_count='0',
            three_prime_support='reads,reads', three_prime_read_count='7,7'))
        _, summary = self.invoke()
        totals = summary['totals']['leap']
        self.assertEqual(totals['annotation_added_bases'], 10)
        self.assertNotIn('mixed_added_bases', totals)
        self.assertEqual(totals['reads_added_bases'], 30)
        self.assertEqual(totals['unique_read_supported_bases'], 30)

    def test_conflicting_sources_on_one_end_are_an_input_error(self):
        self.edit_stats(lambda rows: rows[0].update(five_prime_support='reads,current_annotation', five_prime_read_count='7,0'))
        _, summary = self.invoke(expected=2)
        self.assertIn('one support source', summary['diagnostics'][0]['message'])

    def test_internal_external_and_no_cds_policy_match_existing_processor(self):
        fixture = ROOT / 'tests/data/internal_extensions/internal_left_positive'
        shutil.rmtree(self.run); shutil.copytree(fixture, self.run)
        _, summary = self.invoke(extra=('--strand', 'positive'))
        self.assertEqual(summary['totals']['leap']['added_bases'], 100)
        self.assertEqual(summary['totals']['external']['added_bases'], 50)
        self.assertEqual(summary['totals']['internal']['added_bases'], 50)
        self.assertEqual(summary['totals']['leap']['dropped_no_cds'], 1)
        self.assertEqual(summary['totals']['external']['transcripts'], 2)
        self.assertEqual(summary['totals']['internal']['annotation_added_bases'], 50)

    def test_zero_growth_blacklist_empty_partition_and_gzip(self):
        shutil.rmtree(self.run); self.run = write_run(self.root / 'run', empty_negative=True)
        raw = self.run / 'raw_annotation/current_annotation_1.gff3'
        final = self.run / 'extended_transcripts/1_positive_extended.gff'
        final.write_text(''.join(line.replace(';transcript_type=protein_coding',
            ';transcript_type=protein_coding;leap_tag=LEAP_blacklisted') if '\ttranscript\t' in line else line
            for line in raw.read_text().splitlines(True)))
        self.edit_stats(lambda rows: rows[0].update(transcript_location='1:100-350', five_prime_extension='0',
            three_prime_extension='0', five_prime_support='blacklisted', five_prime_read_count='blacklisted',
            three_prime_support='none', three_prime_read_count='0'))
        for path in (raw, final, self.run / 'extended_transcripts_stats/1_positive_stats.tsv'):
            with gzip.open(str(path) + '.gz', 'wt') as handle: handle.write(path.read_text())
            path.unlink()
        _, summary = self.invoke()
        self.assertEqual(summary['totals']['external']['extended_transcripts'], 0)
        self.assertEqual(summary['totals']['external']['blacklisted_transcripts'], 1)
        self.assertEqual(summary['completed_partitions'], 2)

    def test_missing_ambiguous_malformed_or_stale_inputs_are_visible_incomplete_reports(self):
        stats = self.run / 'extended_transcripts_stats/1_positive_stats.tsv'
        original = stats.read_text()
        self.edit_stats(lambda rows: rows[0].update(five_prime_extension='11'))
        output, summary = self.invoke(expected=2)
        self.assertEqual(summary['status'], 'INPUT_ERROR')
        self.assertEqual(summary['completed_partitions'], 1)
        self.assertIn('extension_length_mismatch', (output / 'diagnostics.csv').read_text())
        stats.write_text(original)
        with gzip.open(str(stats) + '.gz', 'wt') as handle: handle.write(original)
        _, summary = self.invoke(expected=2)
        self.assertIn('exactly one stats', summary['diagnostics'][0]['message'])
        Path(str(stats) + '.gz').unlink(); stats.unlink()
        self.invoke(expected=2)

    def test_datachecks_and_hostile_labels_are_reported_without_html_injection(self):
        checks = self.root / 'checks.json'
        checks.write_text(json.dumps(dict(status='FAIL', counts={'FAIL': 1}, results=[dict(
            chromosome='1', strand='positive', check='test', status='FAIL', message='<unsafe>')])) )
        output, summary = self.invoke(extra=('--datacheck-summary', str(checks), '--label', '</script><script>bad()</script>'))
        self.assertEqual(summary['datachecks']['status'], 'FAIL')
        self.assertEqual(summary['status'], 'COMPLETE')
        self.assertNotIn('</script><script>bad()', (output / 'index.html').read_text())
        self.assertIn('\\u003c/script', (output / 'index.html').read_text())

    def test_existing_outputs_and_linked_outputs_are_protected(self):
        for output in (self.run, self.root / 'linked'):
            if not output.exists(): output.symlink_to(self.run, target_is_directory=True)
            result = subprocess.run([sys.executable, str(SCRIPT), '--leap-output-dir', str(self.run),
                '--output', str(output), '--chromosomes', '1'], capture_output=True, text=True)
            self.assertEqual(result.returncode, 2)
            self.assertIn('new directory', result.stderr)

    def test_interval_arithmetic_and_exact_quantiles(self):
        self.assertEqual(impact.outside([(1, 20), (30, 50)], [(5, 7), (10, 35), (45, 60)]), [(1, 4), (8, 9), (36, 44)])
        self.assertEqual(impact.describe([1, 2, 3, 100])['median'], 2.5)
        self.assertEqual(impact.describe([])['count'], 0)


if __name__ == '__main__': unittest.main()
