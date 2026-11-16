import csv
import gzip
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / 'bin/prepare_manual_review.py'
DATA = ROOT / 'tests/data/internal_extensions'
spec = importlib.util.spec_from_file_location('manual_review', SCRIPT)
review = importlib.util.module_from_spec(spec)
spec.loader.exec_module(review)


def fixture(root, empty_negative=False):
    for strand in ('positive', 'negative'):
        source = DATA / ('internal_left_' + strand)
        for folder in ('extended_transcripts', 'extended_transcripts_stats'):
            target = root / folder
            target.mkdir(parents=True, exist_ok=True)
            for path in (source / folder).iterdir():
                text = path.read_text()
                if strand == 'negative':
                    for old, new in (('anchor.1', 'anchor_negative.1'), ('source.1', 'source_negative.1'),
                                     ('utr_only.1', 'utr_negative.1'), ('geneA', 'geneAnegative'), ('geneB', 'geneBnegative')):
                        text = text.replace(old, new)
                if empty_negative and strand == 'negative':
                    text = '##gff-version 3\n' if folder == 'extended_transcripts' else text.splitlines()[0] + '\n'
                (target / path.name).write_text(text)
    raw = root / 'raw_annotation'
    raw.mkdir()
    text = (DATA / 'internal_left_positive/raw_annotation/current_annotation_1.gff3').read_text()
    if not empty_negative:
        negative = (DATA / 'internal_left_negative/raw_annotation/current_annotation_1.gff3').read_text()
        for old, new in (('anchor.1', 'anchor_negative.1'), ('source.1', 'source_negative.1'),
                         ('utr_only.1', 'utr_negative.1'), ('geneA', 'geneAnegative'), ('geneB', 'geneBnegative')):
            negative = negative.replace(old, new)
        text += negative
    (raw / 'current_annotation_1.gff3').write_text(text)
    return root


@unittest.skipUnless(shutil.which('sort'), 'requires local sort')
class ManualReviewTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(); self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def run_bundle(self, source, extra=()):
        output = self.root / 'review'
        result = subprocess.run([sys.executable, str(SCRIPT), 'run', str(source), str(output),
            '--chromosomes', '1', '--genome-id', 'hg38', '--species', 'homo_sapiens',
            *extra], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        return output

    def test_all_twelve_tracks_preserve_originals_and_show_counterfactual_internal_models(self):
        source = fixture(self.root / 'input')
        before = {p: p.read_bytes() for p in source.rglob('*') if p.is_file()}
        output = self.run_bundle(source)
        annotations = list((output / 'annotations').glob('*/*.gff3'))
        self.assertEqual(len(annotations), 12)
        for group in review.GROUPS:
            pos = list(review.records(output / 'annotations' / group / 'positive.gff3'))
            neg = list(review.records(output / 'annotations' / group / 'negative.gff3'))
            both = list(review.records(output / 'annotations' / group / 'unsplit.gff3'))
            self.assertTrue(pos and neg, group)
            self.assertTrue(all(f[6] == '+' for f in pos))
            self.assertTrue(all(f[6] == '-' for f in neg))
            self.assertCountEqual(both, pos + neg)
        original = list(review.records(source / 'extended_transcripts/1_positive_extended.gff'))
        self.assertCountEqual(list(review.records(output / 'annotations/leaped/positive.gff3')), original)
        def transcript(group):
            return next(f for f in review.records(output / 'annotations' / group / 'positive.gff3')
                        if f[2] == 'transcript' and 'transcript_id=source.1;' in f[8])
        self.assertEqual(transcript('internal')[3:5], ['250', '950'])
        self.assertEqual(transcript('external')[3:5], ['300', '950'])
        self.assertIn('transcript_class=internal_left', transcript('internal')[8])
        self.assertNotIn('utr_only.1', (output / 'annotations/external/positive.gff3').read_text())
        self.assertEqual({p: p.read_bytes() for p in source.rglob('*') if p.is_file()}, before)
        self.assertEqual(json.loads((output / 'review-summary.json').read_text())['read_evidence']['status'], 'not_supplied')
        with (output / 'samplesheet.csv').open() as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(len(rows), 14)
        sheet = next(r for r in rows if r['type']=='spreadsheet')
        self.assertEqual(sheet['sample'], 'LEAP_transcripts')
        self.assertEqual(sheet['data_category'], 'leap_review')
        self.assertEqual(Path(sheet['file']), output / 'spreadsheets/transcripts.csv')
        with Path(sheet['file']).open() as handle:
            transcripts = list(csv.DictReader(handle))
        self.assertEqual(len(transcripts), 6)
        positive = next(r for r in transcripts if r['transcript_id']=='source.1')
        self.assertEqual(positive['jbrowse_label'], '1:250-950')
        self.assertEqual(positive['five_prime_internal_extension'], '50')
        self.assertEqual(positive['three_prime_external_extension'], '50')
        self.assertEqual(sum(r['dropped_no_cds']=='True' for r in transcripts), 2)
        self.assertEqual({r['publication'] for r in rows}, {'public'})
        self.assertTrue(all(r['species'] == 'homo_sapiens' for r in rows))
        with (output / 'sample_metadata.csv').open() as handle:
            self.assertEqual(len(list(csv.DictReader(handle))), 14)
        with (output / 'processing/1_positive/stats/1_positive_stats.tsv').open() as handle:
            stats = {r['transcript_id']: r for r in csv.DictReader(handle, delimiter='\t')}
        self.assertEqual(stats['source.1']['five_prime_extension'], '0')
        self.assertEqual(stats['source.1']['three_prime_extension'], '50')

    def test_empty_strand_gzip_and_fai_become_valid_sizes(self):
        source = fixture(self.root / 'input', empty_negative=True)
        path = source / 'raw_annotation/current_annotation_1.gff3'
        with gzip.open(str(path) + '.gz', 'wt') as handle:
            handle.write(path.read_text())
        path.unlink()
        lengths = self.root / 'reference.fa.gz.fai'
        lengths.write_text('1\t2000\t3\t60\t61\n')
        reference = self.root / "reference's genome.fa"
        reference.write_text('>1\nACGT\n')
        output = self.run_bundle(source, ['--reference-lengths', str(lengths), '--reference-genome', str(reference)])
        self.assertEqual((output / 'reference.chrom.sizes').read_text(), '1\t2000\n')
        with (output / 'track_manifest.tsv').open() as handle:
            manifest = list(csv.DictReader(handle, delimiter='\t'))
        self.assertEqual(len(manifest), 13)
        self.assertTrue(all(r['transcripts'] == '0' and r['included_in_samplesheet'] == 'False'
                            for r in manifest if r['strand'] == 'negative'))
        with (output / 'samplesheet.csv').open() as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(len(rows), 10)
        self.assertEqual(rows[0]['file'], str(reference.resolve()))
        self.assertTrue(all(Path(r['file']).is_file() for r in rows))

    def test_review_spreadsheet_omits_audit_coordinates_and_preserves_other_impact_fields(self):
        source = fixture(self.root / 'input')
        output = self.run_bundle(source)
        impact_output = self.root / 'impact'
        result = subprocess.run([sys.executable, str(ROOT / 'assets/build_leap_impact_report.py'),
            '--leap-output-dir', str(source), '--output', str(impact_output), '--chromosomes', '1'],
            capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        omitted = {'raw_start', 'raw_end', 'leap_start', 'leap_end', 'external_start', 'external_end'}
        with (output/'spreadsheets/transcripts.csv').open() as handle:
            reader = csv.DictReader(handle); fields = reader.fieldnames; review_rows = list(reader)
        with (impact_output/'transcripts.csv').open() as handle:
            reader = csv.DictReader(handle); impact_fields = reader.fieldnames; impact_rows = list(reader)
        self.assertTrue(omitted.isdisjoint(fields))
        self.assertTrue(omitted <= set(impact_fields))
        self.assertIn('jbrowse_label', fields)
        self.assertEqual(fields, [f for f in impact_fields if f not in omitted])
        self.assertEqual(review_rows, [{k:v for k,v in row.items() if k not in omitted} for row in impact_rows])
        with (output / 'samplesheet_offline.csv').open() as handle:
            sheet = next(r for r in csv.DictReader(handle) if r['type']=='spreadsheet')
        self.assertEqual(sheet['publication'], 'offline')
        self.assertEqual(sheet['file'], str(output / 'spreadsheets/transcripts.csv'))
        with (output / 'sample_metadata.csv').open() as handle:
            metadata = next(r for r in csv.DictReader(handle) if r['sample']=='LEAP_transcripts')
        self.assertIn('JBrowse search labels', metadata['description'])

    def test_empty_annotations_still_share_header_only_spreadsheet_without_reference_row(self):
        source = fixture(self.root / 'input', empty_negative=True)
        (source / 'raw_annotation/current_annotation_1.gff3').write_text('##gff-version 3\n')
        (source / 'extended_transcripts/1_positive_extended.gff').write_text('##gff-version 3\n')
        path = source / 'extended_transcripts_stats/1_positive_stats.tsv'
        path.write_text(path.read_text().splitlines()[0] + '\n')
        output = self.run_bundle(source)
        with (output / 'samplesheet.csv').open() as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual([r['type'] for r in rows], ['spreadsheet'])
        with (output / 'spreadsheets/transcripts.csv').open() as handle:
            reader = csv.DictReader(handle)
            self.assertIn('jbrowse_label', reader.fieldnames)
            self.assertEqual(list(reader), [])

    def test_missing_or_invalid_partition_spreadsheet_stops_bundle_publication(self):
        source = DATA / 'internal_left_positive'
        raw = source / 'raw_annotation/current_annotation_1.gff3'
        part = self.root / 'raw_1'; review.raw_partition('1', raw, part)
        extended = self.root / 'extended_1_positive'
        review.extension_partition('1', 'positive', raw,
            source / 'extended_transcripts/1_positive_extended.gff',
            source / 'extended_transcripts_stats/1_positive_stats.tsv', extended)
        sheet = extended / 'transcripts.csv'
        original = sheet.read_bytes(); sheet.unlink()
        with self.assertRaises(OSError):
            review.build_bundle([part, extended], [('1','positive')], dict(genome_id='hg38'), self.root/'missing')
        sheet.write_text('wrong,columns\n')
        with self.assertRaisesRegex(ValueError, 'spreadsheet columns'):
            review.build_bundle([part, extended], [('1','positive')], dict(genome_id='hg38'), self.root/'invalid')
        self.assertFalse((self.root/'missing').exists())
        self.assertFalse((self.root/'invalid').exists())
        sheet.write_bytes(original)

    @unittest.skipUnless(os.environ.get('LEAP_JBROWSER_PIPELINE'), 'set LEAP_JBROWSER_PIPELINE for spreadsheet sharing validation')
    def test_jbrowser_public_and_offline_sharing_preserve_transcript_spreadsheet(self):
        source = fixture(self.root / 'input')
        output = self.run_bundle(source, ['--species', 'hg38'])
        code = '''
import json,shutil,subprocess,sys
from pathlib import Path
pipeline,bundle,work=map(Path,sys.argv[1:])
sys.path.insert(0,str(pipeline/'bin'))
import sample_metadata
import ftp_copy_script_lib as ftp
normalized=work/'samples.json'
sample_metadata.samplesheet_to_json(bundle/'samplesheet.csv',normalized,bundle/'sample_metadata.csv')
data=json.loads(normalized.read_text())
sheet=next(r for r in data['samples'] if r['type']=='spreadsheet')
metadata=data['metadata'][sheet['sample']]
prepared=work/'prepared'
target=prepared/'spreadsheet'/sheet['sample']/Path(sheet['file']).name
target.parent.mkdir(parents=True);shutil.copyfile(sheet['file'],target)
record=dict(sample=sheet['sample'],input_type='spreadsheet',kind='spreadsheet',
    data_file=target.name,species=sheet['species'],data_category=sheet['data_category'],
    ftp_project=sheet['ftp_project'],publication=sheet['publication'],**metadata)
ftp.generate([record],prepared/'pipeline_info',target_dir=str(work/'ftp'),
    public_base_url='https://example.invalid',archive_date='20261116')
subprocess.run(['bash',str(prepared/'pipeline_info/copy_to_ftp.sh')],check=True,capture_output=True)
subprocess.run(['bash',str(prepared/'pipeline_info/build_offline_bundle.sh'),str(work/'offline')],check=True,capture_output=True)
'''
        result = subprocess.run([sys.executable, '-c', code, os.environ['LEAP_JBROWSER_PIPELINE'],
            str(output), str(self.root)], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        original = (output/'spreadsheets/transcripts.csv').read_bytes()
        public = self.root/'ftp/spreadsheet/hg38/leap_review/20261116_1'
        self.assertEqual((public/'transcripts.csv').read_bytes(), original)
        self.assertEqual((self.root/'offline/hg38/LEAP_review/spreadsheet/LEAP_transcripts/transcripts.csv').read_bytes(), original)
        metadata = json.loads((public/'metadata.json').read_text())
        self.assertIn('spreadsheet', metadata['samples'][0]['tags'])
        self.assertEqual(metadata['samples'][0]['source'], 'local')

    def test_missing_ambiguous_and_overwrite_errors(self):
        source = fixture(self.root / 'input')
        gff = source / 'extended_transcripts/1_positive_extended.gff'
        with gzip.open(str(gff) + '.gz', 'wt') as handle:
            handle.write(gff.read_text())
        result = subprocess.run([sys.executable, str(SCRIPT), 'run', str(source), str(self.root / 'bad'),
            '--chromosomes', '1'], capture_output=True, text=True)
        self.assertEqual(result.returncode, 2)
        self.assertIn('ambiguous', result.stderr)
        self.assertFalse((self.root / 'bad').exists())
        with self.assertRaisesRegex(ValueError, 'new directory'):
            review.build_bundle([], [], {}, source)
        Path(str(gff) + '.gz').unlink(); gff.unlink()
        result = subprocess.run([sys.executable, str(SCRIPT), 'run', str(source), str(self.root / 'bad'),
            '--chromosomes', '1'], capture_output=True, text=True)
        self.assertEqual(result.returncode, 2)
        self.assertIn('Missing', result.stderr)

    def test_checked_sort_failure_does_not_publish_bundle(self):
        part = self.root / 'raw_1'
        annotation = DATA / 'internal_left_positive/raw_annotation/current_annotation_1.gff3'
        review.raw_partition('1', annotation, part)
        extended = self.root / 'extended_1_positive'
        review.extension_partition('1', 'positive', annotation,
            DATA / 'internal_left_positive/extended_transcripts/1_positive_extended.gff',
            DATA / 'internal_left_positive/extended_transcripts_stats/1_positive_stats.tsv', extended)
        with patch.object(review.subprocess, 'run', side_effect=subprocess.CalledProcessError(1, 'sort')):
            with self.assertRaises(subprocess.CalledProcessError):
                review.build_bundle([part, extended], [('1', 'positive')], {}, self.root / 'bad')
        self.assertFalse((self.root / 'bad').exists())
        self.assertFalse(list(self.root.glob('.bad-*')))

    def test_sbatch_quotes_settings_and_uses_published_folder_under_slurm(self):
        script = review.launch_script(dict(published_root=str(self.root / "review's folder"),
            jbrowser_pipeline='/somewhere/$(touch SHOULD_NOT_EXIST)', public_base_url='https://example.org'))
        path = self.root / 'run.sbatch'; path.write_text(script)
        result = subprocess.run(['bash', '-n', str(path)], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('#SBATCH', script)
        self.assertIn('REVIEW_DIR=', script)
        self.assertNotIn('dirname -- "$0"', script)
        self.assertIn("JBROWSER_PIPELINE='/somewhere/$(touch SHOULD_NOT_EXIST)'", script)
        self.assertNotIn('eval ', script)

    def test_generated_sbatch_invokes_jbrowser_with_correct_files_and_ftp_roots(self):
        source = fixture(self.root / 'input')
        reference = self.root / "Homo_sapiens.GRCh38.dna primary.fa"
        reference.write_text('>1\nACGT\n')
        pipeline = self.root / "JBrowser's checkout"
        pipeline.mkdir(); (pipeline / 'main.nf').write_text('// test stub\n')
        output = self.run_bundle(source, ['--reference-genome', str(reference), '--jbrowser-pipeline', str(pipeline),
            '--public-base-url', 'https://ftp.example.org/data', '--ftp-target-dir', '/ftp root',
            '--hub-email', 'review@example.org'])
        executables = self.root / 'executables'; executables.mkdir()
        capture = self.root / 'arguments.json'
        fake = executables / 'nextflow'
        fake.write_text('#!' + sys.executable + '\nimport json,sys\nfrom pathlib import Path\nPath(' +
                        repr(str(capture)) + ').write_text(json.dumps(sys.argv[1:]))\n')
        fake.chmod(0o755)
        result = subprocess.run(['bash', str(output / 'run_jbrowser.sbatch'), '-resume'], cwd=self.root,
            env=dict(os.environ, PATH=str(executables) + os.pathsep + os.environ['PATH']), capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        args = json.loads(capture.read_text())
        for key, value in (('--input', str(output / 'samplesheet.csv')),
                           ('--sample_metadata', str(output / 'sample_metadata.csv')),
                           ('--ftp_target_dir', '/ftp root'), ('--public_base_url', 'https://ftp.example.org/data')):
            self.assertEqual(args[args.index(key) + 1], value)
        self.assertEqual(args[-1], '-resume')

    def test_offline_sbatch_refreshes_shared_edits_without_ftp_settings(self):
        source = fixture(self.root / 'input')
        pipeline = self.root / 'JBrowser'; pipeline.mkdir()
        (pipeline / 'main.nf').write_text('// test stub\n')
        output = self.run_bundle(source, ['--jbrowser-pipeline', str(pipeline)])
        with (output / 'samplesheet_offline.csv').open() as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual({r['publication'] for r in rows}, {'offline'})
        self.assertEqual(rows[0]['file'], 'EDIT_ME/reference.fa')
        reference = self.root / "reviewer's reference.fa"
        reference.write_text('>1\nACGT\n')
        with (output / 'samplesheet.csv').open() as handle:
            reader = csv.DictReader(handle); fields = reader.fieldnames; rows = list(reader)
        rows[0]['file'] = str(reference)
        review.write_csv(output / 'samplesheet.csv', fields, rows)
        original = (output / 'samplesheet.csv').read_bytes()
        executables = self.root / 'executables'; executables.mkdir()
        capture = self.root / 'offline-arguments.json'
        fake = executables / 'nextflow'
        fake.write_text('#!' + sys.executable + '\nimport json,sys\nfrom pathlib import Path\nPath(' +
                        repr(str(capture)) + ').write_text(json.dumps(sys.argv[1:]))\n')
        fake.chmod(0o755)
        script = output / 'run_jbrowser_offline.sbatch'
        self.assertEqual(subprocess.run(['bash', '-n', str(script)]).returncode, 0)
        result = subprocess.run(['bash', str(script), '-resume'], cwd=self.root,
            env=dict(os.environ, PATH=str(executables) + os.pathsep + os.environ['PATH']), capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        args = json.loads(capture.read_text())
        self.assertEqual(args[args.index('--input') + 1], str(output / 'samplesheet_offline.csv'))
        self.assertEqual(args[args.index('--outdir') + 1], str(output / 'jbrowser_output_offline'))
        self.assertEqual(args[args.index('-work-dir') + 1], str(output / 'jbrowser_work_offline'))
        self.assertFalse({'--public_base_url', '--ftp_target_dir', '--hub_email'} & set(args))
        self.assertEqual(args[-1], '-resume')
        with (output / 'samplesheet_offline.csv').open() as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual({r['publication'] for r in rows}, {'offline'})
        self.assertEqual(rows[0]['file'], str(reference))
        self.assertEqual((output / 'samplesheet.csv').read_bytes(), original)

    def test_ftp_directory_defaults_to_resolved_genome_and_preserves_explicit_overrides(self):
        part = self.root / 'raw_1'
        source = DATA / 'internal_left_positive'
        raw = source / 'raw_annotation/current_annotation_1.gff3'
        review.raw_partition('1', raw, part)
        extended = self.root / 'extended_1_positive'
        review.extension_partition('1', 'positive', raw,
            source / 'extended_transcripts/1_positive_extended.gff',
            source / 'extended_transcripts_stats/1_positive_stats.tsv', extended)
        cases = [
            (dict(reference_genome='/server/Mus_musculus.GRCm39.dna.fa'), 'mm39', 'mm39'),
            (dict(reference_genome='/server/Homo_sapiens.GRCh38.fa'), 'hg38', 'hg38'),
            (dict(genome_id='custom_assembly'), 'custom_assembly', 'custom_assembly'),
            (dict(reference_genome='/server/Homo_sapiens.GRCh38.fa', genome_id='explicit_id'), 'explicit_id', 'explicit_id'),
            (dict(genome_id='hg38', species='homo_sapiens'), 'hg38', 'homo_sapiens'),
        ]
        for index, (supplied, genome_id, ftp_namespace) in enumerate(cases):
            with self.subTest(settings=supplied):
                output = self.root / ('resolved_' + str(index))
                review.build_bundle([part, extended], [('1', 'positive')], supplied, output)
                settings = json.loads((output / 'review-summary.json').read_text())['settings']
                self.assertEqual(settings['genome_id'], genome_id)
                self.assertEqual(settings['species'], ftp_namespace)
                for sheet in ('samplesheet.csv', 'samplesheet_offline.csv'):
                    with (output / sheet).open() as handle:
                        rows = list(csv.DictReader(handle))
                    self.assertTrue(rows)
                    self.assertTrue(all(r['genome_id']==genome_id and r['species']==ftp_namespace for r in rows))

    def test_internal_track_requires_growth_at_an_internal_end_both_strands(self):
        source = fixture(self.root / 'input')
        for strand in ('positive', 'negative'):
            path = source / 'extended_transcripts' / ('1_' + strand + '_extended.gff')
            name = 'source.1' if strand == 'positive' else 'source_negative.1'
            rows = []
            for fields in review.records(path):
                if 'transcript_id=' + name + ';' in fields[8] and fields[3] == '250':
                    fields[3] = '300'
                rows.append('\t'.join(fields) + '\n')
            path.write_text(''.join(rows))
            path = source / 'extended_transcripts_stats' / ('1_' + strand + '_stats.tsv')
            with path.open() as handle:
                reader = csv.DictReader(handle, delimiter='\t'); fields = reader.fieldnames; stats = list(reader)
            row = next(r for r in stats if r['transcript_id'] == name)
            prime = 'five_prime' if strand == 'positive' else 'three_prime'
            row[prime + '_extension'] = '0'; row[prime + '_support'] = 'none'; row[prime + '_read_count'] = '0'
            row['transcript_location'] = '1:300-950'
            review.write_csv(path, fields, stats, '\t')
        output = self.run_bundle(source)
        self.assertEqual(list(review.records(output / 'annotations/internal/unsplit.gff3')), [])

    def test_multiple_chromosomes_are_merged_without_losing_models(self):
        source = fixture(self.root / 'input')
        for folder in ('raw_annotation', 'extended_transcripts', 'extended_transcripts_stats'):
            for path in list((source / folder).iterdir()):
                text = path.read_text()
                if folder == 'extended_transcripts_stats':
                    text = text.replace('1:', '2:')
                else:
                    text = '\n'.join('2\t' + line[2:] if line.startswith('1\t') else line for line in text.splitlines()) + '\n'
                (path.parent / path.name.replace('_1', '_2').replace('1_positive', '2_positive').replace('1_negative', '2_negative')).write_text(text)
        output = self.run_bundle(source, ['--chromosomes', '1', '2'])
        rows = list(review.records(output / 'annotations/leaped/unsplit.gff3'))
        self.assertEqual({f[0] for f in rows}, {'1', '2'})
        self.assertEqual(sum(f[2] == 'transcript' for f in rows), 12)


if __name__ == '__main__':
    unittest.main()
