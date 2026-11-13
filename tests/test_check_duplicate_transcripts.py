import csv
import gzip
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / 'assets/check_duplicate_transcripts.py'
EXONS = ((100, 150), (200, 250), (300, 350))
CDS = ((120, 150, '0'), (200, 250, '2'), (300, 320, '2'))


def annotation(identifier='t1.1', chromosome='1', strand='+', exons=EXONS,
               cds=CDS, gene='g1', source='test', reverse=False, suffix=''):
    def feature(kind, start, end, phase='.', extra=''):
        return '{}\t{}\t{}\t{}\t{}\t.\t{}\t{}\ttranscript_id={};gene_id={};{}\n'.format(
            chromosome, source, kind, start, end, strand, phase, identifier, gene, extra)
    rows = [feature('transcript', exons[0][0], exons[-1][1],
                    extra='ID={};gene_name=GENE{};transcript_name={};'.format(identifier, suffix, identifier))]
    rows += [feature('exon', start, end, extra='exon_id=e{}{};exon_number={};'.format(index, suffix, index))
             for index, (start, end) in enumerate(exons, 1)]
    rows += [feature('CDS', start, end, phase, extra='ID=CDS:{};'.format(identifier))
             for start, end, phase in cds]
    return ''.join(reversed(rows) if reverse else rows)


class DuplicateTranscriptTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.annotation = self.root / 'annotation.gff3'
        self.output = self.root / 'report.csv'

    def put(self, path, content):
        with (gzip.open if str(path).endswith('.gz') else open)(path, 'wt') as handle:
            handle.write(content)

    def check(self, content=None):
        if content is not None:
            self.put(self.annotation, content)
        result = subprocess.run([sys.executable, str(SCRIPT), str(self.annotation), str(self.output)],
                                capture_output=True, text=True)
        rows = []
        if result.returncode != 2:
            with self.output.open(newline='') as handle:
                rows = list(csv.DictReader(handle))
        return result, rows

    def test_exact_group_members_details_and_summary(self):
        result, rows = self.check(annotation() + annotation('t2.7') + annotation('t3.2') +
                                  annotation('unique', exons=((90, 150),) + EXONS[1:]))
        self.assertEqual(result.returncode, 1, result.stderr)
        self.assertEqual(len(rows), 3)
        self.assertEqual({row['transcript_id'] for row in rows}, {'t1.1', 't2.7', 't3.2'})
        self.assertEqual(len({row['group_id'] for row in rows}), 1)
        for row in rows:
            self.assertEqual(row['group_size'], '3')
            self.assertEqual(row['jbrowse_label'], '1:100-350')
            self.assertEqual(row['exon_coordinates'], '100-150;200-250;300-350')
            self.assertEqual(row['cds_coordinates'], '120-150;200-250;300-320')
            self.assertEqual(row['cds_phases'], '0;2;2')
            self.assertEqual(row['spliced_length'], '153')
            self.assertEqual(row['cds_length'], '103')
            self.assertEqual(row['source_path'], str(self.annotation.resolve()))
            self.assertGreater(int(row['transcript_line']), 0)
        for text in ('transcripts: 4', 'duplicate_groups: 1', 'duplicate_transcripts: 3',
                     'redundant_transcripts: 2', 'unique_transcripts: 1'):
            self.assertIn(text, result.stdout)

    def test_both_strands_and_multiple_groups(self):
        result, rows = self.check(''.join(annotation(identifier, strand=strand)
                                  for strand in ('+', '-') for identifier in ('a' + strand, 'b' + strand)))
        self.assertEqual(result.returncode, 1, result.stderr)
        self.assertEqual(len({row['group_id'] for row in rows}), 2)
        self.assertEqual({row['strand'] for row in rows}, {'+', '-'})

    def test_location_exon_and_cds_differences_are_not_duplicates(self):
        for changes in (dict(chromosome='2'), dict(strand='-'),
                        dict(exons=((99, 150),) + EXONS[1:]),
                        dict(exons=EXONS[:-1] + ((300, 351),)),
                        dict(exons=((100, 151),) + EXONS[1:]),
                        dict(exons=EXONS[:1] + ((199, 250),) + EXONS[2:]),
                        dict(exons=EXONS + ((400, 450),)),
                        dict(cds=((121, 150, '0'),) + CDS[1:]),
                        dict(cds=((120, 150, '1'),) + CDS[1:]),
                        dict(cds=CDS[:2]), dict(cds=())):
            with self.subTest(changes=changes):
                result, rows = self.check(annotation() + annotation('t2.1', **changes))
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual(rows, [])
                self.assertIn('duplicate_groups: 0', result.stdout)

    def test_record_order_metadata_ids_and_versions_do_not_prevent_matches(self):
        first = annotation('t.1')
        second = annotation('t.2', gene='different_gene', source='other', suffix='new', reverse=True)
        second = '\n'.join('\t'.join(line.split('\t')[:8] +
                    [';'.join(reversed(line.split('\t')[8].strip(';').split(';')))])
                    for line in second.strip().splitlines()) + '\n'
        result, rows = self.check(first + second)
        self.assertEqual(result.returncode, 1, result.stderr)
        group = rows[0]['group_id']
        self.assertEqual({r['transcript_id'] for r in rows}, {'t.1', 't.2'})
        self.assertEqual({r['gene_id'] for r in rows}, {'g1', 'different_gene'})
        result, rows = self.check(second + first)
        self.assertEqual(result.returncode, 1, result.stderr)
        self.assertEqual(rows[0]['group_id'], group)

    def test_single_exon_and_noncoding_models(self):
        for cds in ((), ((120, 150, '0'),)):
            result, rows = self.check(annotation(exons=((100, 200),), cds=cds) +
                                      annotation('t2', exons=((100, 200),), cds=cds))
            self.assertEqual(result.returncode, 1, result.stderr)
            self.assertEqual(rows[0]['exon_count'], '1')
            self.assertEqual(rows[0]['cds_count'], str(len(cds)))

    def test_directory_gzip_and_cross_file_groups(self):
        self.annotation = self.root / 'gffs'
        self.annotation.mkdir()
        self.put(self.annotation / 'one.gff.gz', annotation())
        self.put(self.annotation / 'two.gff3', annotation('t2'))
        self.put(self.annotation / 'empty.gff', '##gff-version 3\n')
        (self.annotation / 'ignored.txt').write_text('not a GFF')
        result, rows = self.check()
        self.assertEqual(result.returncode, 1, result.stderr)
        self.assertEqual(len(rows), 2)
        self.assertIn('files: 3', result.stdout)
        self.assertEqual({Path(row['source_path']).name for row in rows}, {'one.gff.gz', 'two.gff3'})
        self.put(self.annotation / 'one.gff3', annotation())
        result, _ = self.check()
        self.assertEqual(result.returncode, 2)
        self.assertIn('ambiguous annotation copies', result.stderr)

    def test_par_copies_are_separate_models(self):
        result, rows = self.check(annotation('t1', chromosome='X') + annotation('t1', chromosome='Y') +
                                  annotation('t2', chromosome='Y'))
        self.assertEqual(result.returncode, 1, result.stderr)
        self.assertEqual({row['chromosome'] for row in rows}, {'Y'})
        self.assertEqual(len(rows), 2)
        self.assertIn('transcripts: 3', result.stdout)

    def test_repeated_id_is_input_error_not_structural_duplicate(self):
        result, _ = self.check(annotation() * 2)
        self.assertEqual(result.returncode, 2)
        self.assertIn('duplicate transcript record', result.stderr)
        self.annotation = self.root / 'gffs'
        self.annotation.mkdir()
        for name in ('one.gff', 'two.gff'):
            self.put(self.annotation / name, annotation())
        result, _ = self.check()
        self.assertEqual(result.returncode, 2)
        self.assertIn('also occurs in', result.stderr)

    def test_invalid_inputs_preserve_existing_report(self):
        content = annotation()
        invalid = (
            'malformed\n', content.replace('100\t150', 'bad\t150'),
            content.replace('100\t150', '0\t150'),
            content.replace('100\t150', '151\t150'),
            content.replace('transcript_id=t1.1;', ''),
            content.replace('gene_id=g1;', 'gene_id=g1;gene_id=other;'),
            content.replace('\t+\t', '\t?\t'),
            content.replace('\t0\t', '\t3\t'),
            content.replace('\t.\t+\t', '\tNaN\t+\t'),
            '\n'.join(content.splitlines()[1:]) + '\n',
            content.splitlines()[0] + '\n',
            content + content.splitlines()[1] + '\n',
            content + content.splitlines()[-1] + '\n',
            content.replace('100\t150', '99\t150'),
            content.replace('120\t150', '90\t150'),
            content.replace('200\t250\t.\t+\t2', '151\t250\t.\t+\t2'),
            content.replace('transcript_id=t1.1;gene_id=g1;exon_id=e1',
                            'transcript_id=t1.1;gene_id=other;exon_id=e1'),
            content.replace('200\t250\t.\t+\t.', '200\t250\t.\t-\t.'),
        )
        for value in invalid:
            with self.subTest(value=value):
                self.output.write_text('existing report\n')
                result, _ = self.check(value)
                self.assertEqual(result.returncode, 2, result.stderr)
                self.assertTrue(result.stderr.startswith('ERROR:'), result.stderr)
                self.assertEqual(self.output.read_text(), 'existing report\n')

    def test_empty_inputs_and_fasta(self):
        for content in ('', '##gff-version 3\n', annotation() + '##FASTA\n>1\nACGT\n'):
            result, rows = self.check(content)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(rows, [])

    def test_missing_directory_and_unreadable_gzip(self):
        result, _ = self.check()
        self.assertEqual(result.returncode, 2)
        self.annotation.mkdir()
        result, _ = self.check()
        self.assertEqual(result.returncode, 2)
        self.assertIn('no GFF/GFF3', result.stderr)
        self.annotation = self.root / 'bad.gff.gz'
        self.annotation.write_bytes(b'not gzip')
        result, _ = self.check()
        self.assertEqual(result.returncode, 2)

    def test_overwrite_protection_including_links(self):
        self.put(self.annotation, annotation())
        before = self.annotation.read_bytes()
        for kind in ('same', 'symlink', 'hardlink'):
            self.output = self.root / (kind + '.csv')
            if kind == 'same':
                self.output = self.annotation
            elif kind == 'symlink':
                self.output.symlink_to(self.annotation)
            else:
                os.link(self.annotation, self.output)
            result, _ = self.check()
            self.assertEqual(result.returncode, 2)
            self.assertIn('overwrite input', result.stderr)
            self.assertEqual(self.annotation.read_bytes(), before)

    def test_real_cox14_annotation_convergence(self):
        fixture = SCRIPT.parent.parent / 'tests/data/duplicate_transcripts/cox14'
        self.annotation = fixture / 'raw.gff3'
        result, rows = self.check()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(rows, [])
        self.annotation = fixture / 'extended.gff3'
        result, rows = self.check()
        self.assertEqual(result.returncode, 1, result.stderr)
        self.assertEqual({row['transcript_id'] for row in rows},
                         {'ENST00000550487.6', 'ENST00001138219.1'})
        self.assertEqual({row['jbrowse_label'] for row in rows}, {'12:50112029-50120463'})
        self.assertEqual({row['cds_coordinates'] for row in rows}, {'50120044-50120217'})
        self.assertEqual({row['cds_phases'] for row in rows}, {'0'})
        self.assertEqual(len({row['group_id'] for row in rows}), 1)


if __name__ == '__main__':
    unittest.main()
