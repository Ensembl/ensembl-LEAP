import csv
import gzip
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / 'assets/check_output_completeness.py'
HEADER = 'gene_id\ttranscript_id\ttranscript_location\tfive_prime_support\n'


def model(tid='t.1', chrom='1', strand='+', kind='protein_coding', utrs=True, cds=True,
          single=False, par=False, tag=''):
    def row(feature, start, end):
        return '{}\ttest\t{}\t{}\t{}\t.\t{}\t.\tID={}:{}{};Parent=g{};gene_id=g;transcript_id={};transcript_type={};{}\n'.format(
            chrom, feature, start, end, strand, feature, tid, '_PAR_Y' if par else '',
            '_PAR_Y' if par else '', tid, kind, tag)
    rows = [row('transcript', 10, 100), row('exon', 10, 30)]
    rows += [row('exon', 70, 100)] if not single else []
    if single:
        rows[0] = row('transcript', 10, 30)
    if utrs:
        rows += [row('five_prime_UTR', 10, 15), row('three_prime_UTR', 25 if single else 90, 30 if single else 100)]
    if cds:
        rows += [row('CDS', 16, 24)]
    return ''.join(rows)


def stats(tid='t.1', chrom='1', support='none'):
    return '{}\t{}\t{}:10-100\t{}\n'.format('g', tid, chrom, support)


class CompletenessTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.raw = self.root / 'raw.gff3'
        self.gff = self.root / 'extended.gff'
        self.tsv = self.root / 'stats.tsv'
        self.output = self.root / 'report.csv'

    def put(self, path, content):
        with (gzip.open if str(path).endswith('.gz') else open)(path, 'wt') as f:
            f.write(content)

    def run_check(self, raw=None, gff=None, tsv=None, extra=(), paths=None):
        if paths is None:
            for path, value in ((self.raw, raw), (self.gff, gff), (self.tsv, tsv)):
                if value is not None:
                    self.put(path, value)
            paths = self.raw, self.gff, self.tsv
        result = subprocess.run([sys.executable, str(SCRIPT), str(paths[0]), str(paths[1]),
                                 str(self.output), '--extension-stats', str(paths[2]),
                                 '--chromosomes', '1', '--strand', 'positive'] + list(extra),
                                capture_output=True, text=True)
        rows = []
        if self.output.exists() and result.returncode != 2:
            with self.output.open() as handle:
                rows = list(csv.DictReader(handle))
        return result, rows

    def directory_paths(self, raw=None, gff=None, tsv=None):
        dirs = [self.root / x for x in ('raws', 'gffs', 'statss')]
        for d in dirs:
            d.mkdir(exist_ok=True)
        for d, name, value in zip(dirs, ('current_annotation_1.gff3', '1_positive_extended.gff', '1_positive_stats.tsv'), (raw, gff, tsv)):
            if value is not None:
                self.put(d / name, value)
        return dirs

    def test_retained_zero_extension_blacklist_and_both_strands(self):
        for strand, label in (('+', 'positive'), ('-', 'negative')):
            result, rows = self.run_check(model(strand=strand), model(strand=strand, tag='leap_tag=LEAP_blacklisted;'),
                HEADER + stats(support='blacklisted'), extra=['--strand', label])
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(rows[0]['status'], 'retained')
            self.assertIn('expected_retained: 1', result.stdout)

    def test_exclusions_and_no_cds_stages(self):
        for kwargs, reason in (({'kind': 'lncRNA'}, 'non_protein_coding'), ({'utrs': False}, 'missing_five_prime_UTR'),
                               ({'single': True}, 'single_exon')):
            result, rows = self.run_check(model(**kwargs), '', HEADER)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn(reason, rows[0]['exclusion_reasons'])
            result, rows = self.run_check(gff=model(**kwargs), tsv=HEADER + stats())
            self.assertEqual(result.returncode, 1)
            self.assertEqual(rows[0]['status'], 'unexpected_retained')
        result, rows = self.run_check(model(cds=False), model(cds=False), HEADER + stats())
        self.assertEqual(result.returncode, 0, result.stderr)
        result, rows = self.run_check(gff='', tsv=HEADER, extra=['--stage', 'processed'])
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(rows[0]['exclusion_reasons'], 'processed_no_CDS')
        result, rows = self.run_check(model(), model(tag='transcript_class=internal;'), HEADER + stats(), extra=['--stage', 'processed'])
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_missing_independently_and_from_both(self):
        for gff, tsv, detail in (('', HEADER, 'missing_from_gff;missing_from_stats'),
                                 ('', HEADER + stats(), 'missing_from_gff'),
                                 (model(), HEADER, 'missing_from_stats')):
            result, rows = self.run_check(model(), gff, tsv)
            self.assertEqual(result.returncode, 1, result.stderr)
            self.assertEqual(rows[0]['details'], detail)
            self.assertIn('missing_expected: 1', result.stdout)

    def test_extra_and_version_changes(self):
        for gff, tsv in ((model() + model('extra'), HEADER + stats()), (model(), HEADER + stats() + stats('extra')),
                          (model('t.2'), HEADER + stats('t.2'))):
            result, rows = self.run_check(model(), gff, tsv)
            self.assertEqual(result.returncode, 1, result.stderr)
            self.assertTrue(any(r['status'] == 'unexpected_retained' for r in rows))
            self.assertTrue(any(r['exclusion_reasons'] == 'absent_from_raw_exact_id' for r in rows))

    def test_raw_scope_ignored_output_scope_fails(self):
        raw = model() + model('other', strand='-') + model('foreign', chrom='2')
        result, rows = self.run_check(raw, model(), HEADER + stats())
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(len(rows), 1)
        for gff, tsv in ((model(strand='-'), HEADER + stats()), (model(), HEADER + stats(chrom='2'))):
            result, rows = self.run_check(raw, gff, tsv)
            self.assertEqual(result.returncode, 1, result.stderr)
            self.assertTrue(any(r['status'].endswith('outside_file_scope') for r in rows))
        result, _ = self.run_check(model(strand='-'), '', HEADER)
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_directory_both_strands_par_and_gzip(self):
        paths = self.directory_paths()
        for chrom in ('X', 'Y'):
            self.put(paths[0] / ('current_annotation_' + chrom + '.gff3.gz'), model(chrom=chrom, par=chrom == 'Y') + model('negative', chrom=chrom, strand='-'))
            for label, strand, tid in (('positive', '+', 't.1'), ('negative', '-', 'negative')):
                self.put(paths[1] / (chrom + '_' + label + '_extended.gff.gz'), model(tid, chrom, strand, par=chrom == 'Y'))
                self.put(paths[2] / (chrom + '_' + label + '_stats.tsv.gz'), HEADER + stats(tid, chrom))
        result, rows = self.run_check(paths=paths, extra=['--chromosomes', 'X', 'Y', '--strand', 'both'])
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(len(rows), 4)
        self.assertTrue(any(r['raw_gff_id'].endswith('_PAR_Y') for r in rows))
        self.assertIn('retained: 4', result.stdout)

    def test_directory_missing_files_and_ambiguous_duplicates(self):
        paths = self.directory_paths(model(), model(), HEADER + stats())
        result, rows = self.run_check(paths=paths, extra=['--strand', 'both'])
        self.assertEqual(result.returncode, 1, result.stderr)
        self.assertEqual({r['status'] for r in rows if r['row_type'] == 'file'}, {'missing_gff_file', 'missing_stats_file'})
        self.put(paths[1] / '1_positive_extended.gff.gz', model())
        result, _ = self.run_check(paths=paths)
        self.assertEqual(result.returncode, 2)
        self.assertIn('ambiguous', result.stderr)
        (paths[1] / '1_positive_extended.gff.gz').unlink()
        (paths[0] / 'current_annotation_1.gff3').unlink()
        result, rows = self.run_check(paths=paths)
        self.assertEqual(result.returncode, 1, result.stderr)
        self.assertTrue(any(r['status'] == 'raw_unavailable' for r in rows))

    def test_zero_eligible_still_requires_files(self):
        paths = self.directory_paths(model(kind='lncRNA'), '', HEADER)
        result, _ = self.run_check(paths=paths)
        self.assertEqual(result.returncode, 0, result.stderr)
        (paths[1] / '1_positive_extended.gff').unlink()
        result, rows = self.run_check(paths=paths)
        self.assertEqual(result.returncode, 1)
        self.assertTrue(any(r['status'] == 'missing_gff_file' for r in rows))

    def test_invalid_raw_models(self):
        raw = model()
        for invalid, message in ((raw.replace('transcript_type=protein_coding;', ''), 'transcript_type'),
                (raw.replace('transcript_id=t.1;', ''), 'transcript_id'),
                (raw + raw.splitlines(True)[0], 'duplicate'),
                (raw + raw.splitlines(True)[1], 'duplicate'),
                (raw.replace('\texon\t', '\tother\t'), 'no exons'),
                (''.join(raw.splitlines(True)[1:]), 'orphan'),
                (raw.replace('\t70\t100\t', '\t20\t100\t'), 'overlapping'),
                (raw.replace('\t10\t30\t', '\t0\t30\t'), 'invalid'),
                (raw.replace('gene_id=g;', 'gene_id=other;', 1), 'inconsistent'),
                ('# empty\n', 'no transcript'), ('1\tbad\n', '9 GFF3')):
            with self.subTest(message=message):
                self.output.unlink(missing_ok=True)
                result, _ = self.run_check(invalid, model(), HEADER + stats())
                self.assertEqual(result.returncode, 2, result.stderr)
                self.assertIn(message, result.stderr)
                self.assertFalse(self.output.exists())

    def test_invalid_stats_and_final(self):
        for gff, tsv in ((model() + model(), HEADER + stats()), (model(), HEADER + stats() + stats()),
                         (model(), 'wrong\n'), (model(), HEADER + stats().replace('10-100', '100-10')),
                         (''.join(model().splitlines(True)[1:]), HEADER)):
            result, _ = self.run_check(model(), gff, tsv)
            self.assertEqual(result.returncode, 2, result.stderr)

    def test_partial_single_missing_files_and_raw_scope_errors(self):
        raw = ''.join(line for line in model().splitlines(True) if '\tthree_prime_UTR\t' not in line)
        result, rows = self.run_check(raw, '', HEADER)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(rows[0]['exclusion_reasons'], 'missing_three_prime_UTR')
        self.gff.unlink()
        result, _ = self.run_check()
        self.assertEqual(result.returncode, 1, result.stderr)
        result, _ = self.run_check(model(chrom='2'), model(), HEADER + stats())
        self.assertEqual(result.returncode, 2)
        self.assertIn('selected raw chromosome', result.stderr)
        result, _ = self.run_check(model().replace('\t+\t', '\t.\t'), model(), HEADER + stats())
        self.assertEqual(result.returncode, 2)

    def test_ambiguous_raw_stats_and_missing_directory(self):
        paths = self.directory_paths(model(), model(), HEADER + stats())
        for path, content in ((paths[0] / 'current_annotation_1.gff', model()),
                              (paths[2] / '1_positive_stats.tsv.gz', HEADER + stats())):
            self.put(path, content)
            result, _ = self.run_check(paths=paths)
            self.assertEqual(result.returncode, 2)
            self.assertIn('ambiguous', result.stderr)
            path.unlink()
        paths[1] = self.root / 'absent-gff-directory'
        result, rows = self.run_check(paths=paths)
        self.assertEqual(result.returncode, 1, result.stderr)
        self.assertTrue(any(r['status'] == 'missing_gff_file' for r in rows))

    def test_reordered_records_and_invalid_syntax(self):
        reordered = ''.join(reversed(model().splitlines(True)))
        result, _ = self.run_check(reordered, reordered, HEADER + stats())
        self.assertEqual(result.returncode, 0, result.stderr)
        for invalid, message in ((model().replace('\t.\t+\t.', '\tbad\t+\t.'), 'score'),
                                 (model().replace('\t+\t.\t', '\t+\t9\t'), 'phase')):
            result, _ = self.run_check(invalid, model(), HEADER + stats())
            self.assertEqual(result.returncode, 2)
            self.assertIn(message, result.stderr)

    def test_directory_overwrite_unselected_and_missing_slots(self):
        paths = self.directory_paths(model(), model(), HEADER + stats())
        unselected = paths[0] / 'current_annotation_2.gff3'
        self.put(unselected, model(chrom='2'))
        self.output = unselected
        result, _ = self.run_check(paths=paths)
        self.assertEqual(result.returncode, 2)
        self.assertEqual(unselected.read_text(), model(chrom='2'))
        missing = paths[1] / '1_positive_extended.gff'
        missing.unlink()
        self.output = missing
        result, _ = self.run_check(paths=paths)
        self.assertEqual(result.returncode, 2)
        self.assertFalse(missing.exists())

    def test_unreadable_gzip_returns_input_error(self):
        self.put(self.raw, model())
        self.put(self.tsv, HEADER + stats())
        compressed = self.root / 'broken.gff.gz'
        # A gzip header followed by an invalid DEFLATE block.
        compressed.write_bytes(bytes.fromhex('1f8b080000000000000007'))
        result, _ = self.run_check(paths=(self.raw, compressed, self.tsv))
        self.assertEqual(result.returncode, 2, result.stderr)
        self.assertIn('ERROR:', result.stderr)
        self.assertNotIn('Traceback', result.stderr)

    def test_single_scope_and_overwrite_protection(self):
        result, _ = self.run_check(model(), model(), HEADER + stats(), extra=['--strand', 'both'])
        self.assertEqual(result.returncode, 2)
        for target in (self.raw, self.gff, self.tsv):
            for link in ('direct', 'symlink', 'hardlink'):
                self.output = self.root / ('report-' + link)
                self.output.unlink(missing_ok=True)
                if link == 'direct':
                    self.output = target
                elif link == 'symlink':
                    self.output.symlink_to(target)
                else:
                    os.link(target, self.output)
                result, _ = self.run_check(model(), model(), HEADER + stats())
                self.assertEqual(result.returncode, 2, result.stderr)
                self.assertIn('overwrite', result.stderr)


if __name__ == '__main__':
    unittest.main()
