import csv
import gzip
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / 'assets/check_transcript_preservation.py'


def annotation(identifier='t1.1', strand='+', single=False, noncoding=False):
    def record(kind, start, end, number=None):
        attrs = dict(ID=(identifier if kind == 'transcript' else kind + ':' + identifier + ':' + str(number or start)),
                     Parent='g1' if kind == 'transcript' else identifier,
                     gene_id='g1', transcript_id=identifier, transcript_type='protein_coding')
        if number is not None:
            attrs.update(exon_number=str(number), exon_id='e' + str(number))
        return '\t'.join(['1', 'test', kind, str(start), str(end), '.', strand,
                          '0' if kind in ('CDS', 'start_codon', 'stop_codon') else '.',
                          ';'.join(k + '=' + v for k, v in attrs.items())]) + '\n'
    exons = [(100, 350)] if single else [(100, 150), (200, 250), (300, 350)]
    rows = [record('transcript', 100, 350)]
    rows += [record('exon', start, end, index + 1 if strand == '+' else len(exons) - index)
             for index, (start, end) in enumerate(exons)]
    if not noncoding:
        low, high = ('five_prime_UTR', 'three_prime_UTR') if strand == '+' else ('three_prime_UTR', 'five_prime_UTR')
        if single:
            rows += [record(low, 100, 124), record('CDS', 125, 325), record(high, 326, 350)]
            codons = [(125, 127), (323, 325)]
        else:
            rows += [record(low, 100, 150), record(low, 200, 204), record('CDS', 205, 245),
                     record('CDS', 300, 320), record(high, 321, 350)]
            codons = [(205, 207), (318, 320)]
        if strand == '-':
            codons.reverse()
        rows += [record('start_codon', *codons[0]), record('stop_codon', *codons[1])]
    return ''.join(rows)


def edit(text, predicate, column=None, value=None, attributes=None):
    result = []
    for line in text.splitlines():
        fields = line.split('\t')
        if predicate(fields):
            if column is not None:
                fields[column] = str(value)
            if attributes is not None:
                attrs = dict(v.split('=', 1) for v in fields[8].split(';') if v)
                attrs.update(attributes)
                fields[8] = ';'.join(k + '=' + v for k, v in attrs.items())
        result.append('\t'.join(fields) + '\n')
    return ''.join(result)


def extend(text, left=90, right=360):
    strand = text.splitlines()[0].split('\t')[6]
    low, high = ('five_prime_UTR', 'three_prime_UTR') if strand == '+' else ('three_prime_UTR', 'five_prime_UTR')
    text = edit(text, lambda f: f[2] in ('transcript', 'exon', low) and f[3] == '100', 3, left)
    return edit(text, lambda f: f[2] in ('transcript', 'exon', high) and f[4] == '350', 4, right)


class TranscriptPreservationTests(unittest.TestCase):
    def run_check(self, raw, extended, compressed=False):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            inputs = []
            for label, content in [('raw', raw), ('extended', extended)]:
                path = root / (label + ('.gff3.gz' if compressed else '.gff3'))
                with (gzip.open if compressed else open)(path, 'wt') as handle:
                    handle.write(content)
                inputs.append(str(path))
            output = root / 'report.csv'
            result = subprocess.run([sys.executable, str(SCRIPT)] + inputs + [str(output)], capture_output=True, text=True)
            rows = None
            if output.exists():
                with output.open(newline='') as handle:
                    rows = list(csv.DictReader(handle))
            return result, rows

    def test_valid_extensions_both_strands_and_single_exons(self):
        for strand in ('+', '-'):
            for single in (False, True):
                for left, right in [(100, 350), (90, 350), (100, 360), (90, 360)]:
                    with self.subTest(strand=strand, single=single, left=left, right=right):
                        raw = annotation(strand=strand, single=single)
                        result, rows = self.run_check(raw, extend(raw, left, right), compressed=True)
                        self.assertEqual(result.returncode, 0, result.stderr)
                        self.assertEqual(rows[0]['status'], 'unchanged' if (left, right) == (100, 350) else 'extended')
                        self.assertEqual(rows[0]['extended_start'], str(left))
                        self.assertEqual(rows[0]['extended_end'], str(right))
                        self.assertEqual(rows[0]['feature_differences'], '{}')
                        self.assertIn('PASS: 0 failing transcript(s)', result.stdout)

    def test_record_attribute_order_and_only_generated_metadata_ignored(self):
        raw = annotation()
        new = extend(raw)
        new = edit(new, lambda f: True, attributes={'leap_tag': 'LEAP_3,LEAP_5', 'transcript_class': 'external'})
        rows = []
        for line in new.splitlines():
            fields = line.split('\t');fields[8] = ';'.join(reversed(fields[8].split(';')))
            rows.append('\t'.join(fields) + '\n')
        gene = '1\ttest\tgene\t1\t1000\t.\t+\t.\tID=g1;gene_id=g1\n'
        result, report = self.run_check('##gff-version 3\n' + gene + raw,
                                        gene.replace('\t1000\t', '\t2000\t') + ''.join(reversed(rows)))
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(report[0]['status'], 'extended')

    def test_every_protected_exon_boundary(self):
        for strand in ('+', '-'):
            raw = annotation(strand=strand)
            for start, end, col, value in [(100, 150, 4, 151), (200, 250, 3, 199),
                                            (200, 250, 4, 251), (300, 350, 3, 299)]:
                with self.subTest(strand=strand, start=start, column=col):
                    new = edit(raw, lambda f: f[2] == 'exon' and f[3:5] == [str(start), str(end)], col, value)
                    result, rows = self.run_check(raw, new)
                    self.assertEqual(result.returncode, 1, result.stderr)
                    self.assertEqual(rows[0]['feature_records_differ'], 'True')
                    self.assertIn('exon', rows[0]['feature_differences'])

    def test_cds_codons_internal_utrs_and_unknown_features(self):
        raw = annotation()
        for kind, start, col, value in [('CDS', 205, 3, 206), ('CDS', 300, 4, 319),
                                       ('start_codon', 205, 4, 206), ('stop_codon', 318, 3, 319),
                                       ('five_prime_UTR', 200, 3, 201), ('five_prime_UTR', 100, 4, 149),
                                       ('three_prime_UTR', 321, 3, 322)]:
            with self.subTest(kind=kind, start=start, column=col):
                result, rows = self.run_check(raw, edit(raw, lambda f: f[2] == kind and f[3] == str(start), col, value))
                self.assertEqual(result.returncode, 1, result.stderr)
                self.assertIn(kind, rows[0]['feature_differences'])
        unknown = raw + '1\ttest\tcustom_feature\t210\t220\t.\t+\t.\tID=u;Parent=t1.1;transcript_id=t1.1;gene_id=g1\n'
        self.assertEqual(self.run_check(unknown, extend(unknown))[0].returncode, 0)
        self.assertEqual(self.run_check(unknown, edit(unknown, lambda f: f[2] == 'custom_feature', 4, 221))[0].returncode, 1)
        # A new feature wholly in the added part must not disappear during clipping.
        added = extend(raw) + '1\ttest\tcustom_feature\t90\t99\t.\t+\t.\tID=u;Parent=t1.1;transcript_id=t1.1;gene_id=g1\n'
        result, rows = self.run_check(raw, added)
        self.assertEqual(result.returncode, 1)
        self.assertEqual(rows[0]['feature_count_mismatch'], 'True')

    def test_source_score_phase_and_attributes_preserved(self):
        raw = annotation()
        changes = [(1, 'other_source', None), (5, '2.5', None), (7, '1', None),
                   (None, None, {'protein_id': 'p2'}), (None, None, {'ID': 'changed'}),
                   (None, None, {'exon_number': '99'}), (None, None, {'exon_id': 'changed'}),
                   (None, None, {'tag': 'changed'})]
        for column, value, attrs in changes:
            with self.subTest(column=column, attributes=attrs):
                result, rows = self.run_check(raw, edit(raw, lambda f: f[2] == 'CDS' and f[3] == '205', column, value, attrs))
                self.assertEqual(result.returncode, 1, result.stderr)
                details = json.loads(rows[0]['feature_differences'])
                self.assertTrue(details['current_only'])
                self.assertTrue(details['normalized_extended_only'])
        new = edit(raw, lambda f: f[2] == 'transcript', attributes={'ID': 'transcript:t1.1'})
        new = edit(new, lambda f: f[2] != 'transcript', attributes={'Parent': 'transcript:t1.1'})
        self.assertEqual(self.run_check(raw, new)[0].returncode, 1)

    def test_feature_removal_and_exon_addition(self):
        raw = annotation()
        for kind in ('CDS', 'start_codon', 'stop_codon', 'five_prime_UTR'):
            new = ''.join(line + '\n' for line in raw.splitlines() if line.split('\t')[2] != kind)
            result, rows = self.run_check(raw, new)
            self.assertEqual(result.returncode, 1, result.stderr)
            self.assertEqual(rows[0]['feature_count_mismatch'], 'True')
        extra = '1\ttest\texon\t260\t270\t.\t+\t.\tID=new;Parent=t1.1;gene_id=g1;transcript_id=t1.1;exon_number=4\n'
        result, rows = self.run_check(raw, raw + extra)
        self.assertEqual(result.returncode, 1, result.stderr)
        self.assertEqual(rows[0]['feature_count_mismatch'], 'True')
        two = ''.join(line + '\n' for line in raw.splitlines() if not 200 <= int(line.split('\t')[3]) <= 250)
        self.assertEqual(self.run_check(two, extend(two))[0].returncode, 0)
        self.assertEqual(self.run_check(raw, two)[0].returncode, 1)

    def test_inconsistent_terminal_extensions(self):
        for strand in ('+', '-'):
            raw = annotation(strand=strand)
            # Transcript only; exon only; UTR not extended; UTR extended too little.
            cases = [(edit(raw, lambda f: f[2] == 'transcript', 3, 90), 1),
                     (edit(raw, lambda f: f[2] == 'exon' and f[3] == '100', 3, 90), 2),
                     (edit(extend(raw), lambda f: f[2].endswith('_UTR') and f[3] == '90', 3, 100), 1),
                     (edit(extend(raw), lambda f: f[2].endswith('_UTR') and f[3] == '90', 3, 91), 1),
                     (edit(extend(raw), lambda f: f[2].endswith('_UTR') and f[4] == '360', 4, 350), 1)]
            for new, expected_code in cases:
                result, rows = self.run_check(raw, new)
                self.assertEqual(result.returncode, expected_code, result.stderr)
                if rows:
                    self.assertEqual(rows[0]['status'], 'changed_features')
                    self.assertEqual(rows[0]['inconsistent_extension'], 'True')
        raw = annotation(noncoding=True)
        self.assertEqual(self.run_check(raw, extend(raw))[0].returncode, 0)
        for left, right in [(101, 360), (90, 349)]:
            raw = annotation(single=True)
            result, rows = self.run_check(raw, extend(raw, left, right))
            self.assertEqual(result.returncode, 1, result.stderr)
            self.assertEqual(rows[0]['shortened'], 'True')

    def test_location_and_gene_changes(self):
        raw = annotation()
        for column, value in [(0, '2'), (6, '-')]:
            result, rows = self.run_check(raw, edit(raw, lambda f: True, column, value))
            self.assertEqual(result.returncode, 1, result.stderr)
            self.assertEqual(rows[0]['location_mismatch'], 'True')
        result, rows = self.run_check(raw, edit(raw, lambda f: True, attributes={'gene_id': 'g2'}))
        self.assertEqual(result.returncode, 1)
        self.assertIn('g2', rows[0]['feature_differences'])

    def test_unmatched_exact_versions_and_summary(self):
        result, rows = self.run_check(annotation() + annotation('filtered'), annotation())
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(rows[0]['status'], 'missing_from_extended')
        self.assertEqual(rows[0]['extended_start'], '')
        result, rows = self.run_check(annotation(), annotation('t1.2'))
        self.assertEqual(result.returncode, 1, result.stderr)
        self.assertEqual([r['status'] for r in rows], ['missing_from_extended', 'missing_from_raw'])
        self.assertEqual([r['failed'] for r in rows], ['False', 'True'])
        self.assertIn('FAIL: 1 failing transcript(s)', result.stdout)
        self.assertIn('missing_from_extended: 1', result.stdout)
        self.assertIn('missing_from_raw: 1', result.stdout)

    def test_invalid_inputs(self):
        raw = annotation();lines = raw.splitlines(keepends=True)
        cases = [('', 'no transcript'), ('# comment\n', 'no transcript'),
                 (raw + lines[0], 'duplicate transcript'), (raw + lines[1], 'duplicate feature'),
                 (raw.replace('transcript_id=t1.1', 'wrong=t1.1'), 'missing transcript_id'),
                 (raw.replace('gene_id=g1', 'wrong=g1'), 'missing gene_id'),
                 (raw.replace('gene_id=g1', 'broken'), 'malformed'),
                 (raw.replace('gene_id=g1', 'gene_id=g1;gene_id=g1'), 'duplicate attribute'),
                 (raw.replace('\t100\t150\t', '\t0\t150\t'), 'invalid coordinates'),
                 (raw.replace('\t100\t150\t', '\t151\t150\t'), 'invalid coordinates'),
                 (raw.replace('\t100\t150\t', '\tx\t150\t'), 'integers'),
                 (raw.replace('\t+\t', '\tx\t'), 'invalid chromosome'),
                 (raw.replace('\t+\t0\t', '\t+\t3\t'), 'invalid phase'),
                 (edit(raw, lambda f: f[2] == 'CDS', 5, 'nan'), 'invalid score'),
                 (lines[0], 'no exons'), (''.join(lines[1:]), 'orphan child'),
                 (edit(raw, lambda f: f[2] == 'exon', attributes={'Parent': 'absent'}), 'orphan Parent'),
                 (edit(raw, lambda f: f[2] == 'exon', 0, '2'), 'inconsistent feature location'),
                 (edit(raw, lambda f: f[2] == 'exon' and f[3] == '100', 4, 200), 'overlapping'),
                 (edit(raw, lambda f: f[2] == 'CDS' and f[3] == '205', 4, 251), 'containing exon'),
                 ('1\ttest\ttranscript\n', '9 GFF3 columns')]
        for invalid, message in cases:
            for before, after in [(invalid, raw), (raw, invalid)]:
                with self.subTest(message=message, raw=before is invalid):
                    result, rows = self.run_check(before, after)
                    self.assertEqual(result.returncode, 2, result.stderr)
                    self.assertIn(message, result.stderr)
                    self.assertIsNone(rows)

    def test_files_input_protection_and_existing_report(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory);raw = root / 'raw.gff3';raw.write_text(annotation())
            alias = root / 'alias.gff3';alias.symlink_to(raw)
            hardlink = root / 'hardlink.gff3';hardlink.hardlink_to(raw)
            for output in [raw, alias, hardlink, root / 'missing/report.csv']:
                result = subprocess.run([sys.executable, str(SCRIPT), str(raw), str(raw), str(output)], capture_output=True, text=True)
                self.assertEqual(result.returncode, 2)
                self.assertEqual(raw.read_text(), annotation())
            output = root / 'report.csv';output.write_text('previous report')
            for missing in [root / 'absent', root / 'bad.gff3.gz']:
                if missing.suffix == '.gz':missing.write_bytes(b'not gzip')
                result = subprocess.run([sys.executable, str(SCRIPT), str(raw), str(missing), str(output)], capture_output=True, text=True)
                self.assertEqual(result.returncode, 2)
                self.assertEqual(output.read_text(), 'previous report')


if __name__ == '__main__':
    unittest.main()
