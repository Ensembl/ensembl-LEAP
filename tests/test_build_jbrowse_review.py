"""Integration checks for the reusable review-bundle builder."""

import csv
import gzip
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
BUILDER = ROOT / 'assets/build_jbrowse_review.py'
CHECK = ROOT / 'assets/check_cds_overlaps.py'


from leap_test_context import write_context


def model(start=200, target=False):
    identifier, gene = ('target', 'geneA') if target else ('source', 'geneA')
    features = [('transcript', start, 450), ('exon', start, 250), ('exon', 400, 450)]
    features += [('CDS', 160, 180)] if target else [('five_prime_UTR', start, 220)]
    return ''.join(
        '1\ttest\t{}\t{}\t{}\t.\t+\t.\tID={};Parent=geneA;transcript_id={};gene_id={};gene_name=Gene%20A;transcript_type=protein_coding;\n'.format(
            kind,
            left,
            right,
            identifier if kind == 'transcript' else identifier + kind,
            identifier,
            gene,
        )
        for kind, left, right in features
    )


@unittest.skipUnless(
    all(shutil.which(tool) for tool in ('sort', 'bgzip', 'tabix')),
    'requires sort, bgzip and tabix',
)
class BuildReviewTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.directory = Path(self.temp.name)
        self.raw = self.directory / 'raw.gff3'
        self.extended = self.directory / 'extended.gff3'
        self.raw.write_text(model() + model(140, target=True))
        self.extended.write_text(model(150))
        self.report = self.directory / 'report.csv'
        context = write_context(self.directory, self.raw, self.extended)
        result = subprocess.run(
            [
                sys.executable,
                str(CHECK),
                str(self.raw),
                str(self.extended),
                str(self.report),
                '--leap-output-dir',
                str(context),
            ],
            capture_output=True,
            text=True,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.output = self.directory / 'review bundle'

    def run_builder(self, *arguments, raw=None, extended=None):
        return subprocess.run(
            [
                sys.executable,
                str(BUILDER),
                '--raw-annotation',
                str(raw or self.raw),
                '--extended-annotation',
                str(extended or self.extended),
                '--output',
                str(self.output),
            ]
            + list(arguments),
            capture_output=True,
            text=True,
        )

    def test_build_config_indexes_and_coordinate_preservation(self):
        self.extended.write_text(
            self.extended.read_text().replace(
                'gene_name=Gene%20A;',
                'gene_name=Gene%20A;transcript_class=internal_left;',
            )
        )
        result = self.run_builder('--reports', str(self.report))
        self.assertEqual(result.returncode, 0, result.stderr)
        config = json.loads((self.output / 'config.json').read_text())
        self.assertEqual(
            config,
            json.loads((self.output / 'LEAP-overlap-review.jbrowse').read_text()),
        )
        self.assertEqual(len(config['tracks']), 10)
        self.assertIn('shared internal splice site', config['defaultSession']['name'])
        self.assertEqual(len(config['preConfiguredSessions']), 1)
        for track in config['tracks']:
            adapter = track['adapter']
            path = Path(adapter['gffGzLocation']['localPath'])
            self.assertTrue(path.is_file())
            self.assertTrue(Path(adapter['index']['location']['localPath']).is_file())
            self.assertEqual(path.parent, self.output.resolve() / 'data')
            self.assertNotIn('-build-', str(path))
        records = subprocess.check_output(
            ['tabix', str(self.output / 'data/raw.gff3.gz'), '1:1-500'], text=True
        )
        transcript = next(
            line.split('\t')
            for line in records.splitlines()
            if '\ttranscript\t' in line and 'ID=source;' in line
        )
        self.assertEqual(transcript[3:5], ['200', '450'])
        self.assertNotIn('Parent=', transcript[8])
        self.assertIn('Name=Gene%20A%20source', transcript[8])
        processed = subprocess.check_output(
            ['tabix', str(self.output / 'data/extended.gff3.gz'), '1:1-500'], text=True
        )
        self.assertIn('transcript_class=internal_left', processed)
        marker = subprocess.check_output(
            ['tabix', str(self.output / 'data/shared_site_allowed.gff3.gz'), '1:1-500'],
            text=True,
        )
        self.assertIn('\t160\t180\t', marker)
        self.assertIn('shared_splice_sites=end:250', marker)
        self.assertFalse(list(self.output.rglob('*.unsorted.gff3')))
        self.assertIn('1:1..2180', (self.output / 'review-loci.tsv').read_text())
        summary = json.loads((self.output / 'build-summary.json').read_text())
        self.assertEqual(summary['feature_counts']['shared_site_allowed'], 1)

    def test_mouse_reference_and_chromosome_aliases(self):
        result = self.run_builder('--reports', str(self.report), '--assembly', 'GRCm39')
        self.assertEqual(result.returncode, 0, result.stderr)
        config = json.loads((self.output / 'config.json').read_text())
        assembly = config['assemblies'][0]
        self.assertEqual(assembly['name'], 'GRCm39')
        self.assertEqual(assembly['sequence']['adapter']['type'], 'TwoBitAdapter')
        self.assertTrue(
            assembly['sequence']['adapter']['twoBitLocation']['uri'].endswith(
                '/mm39.2bit'
            )
        )
        aliases = Path(assembly['refNameAliases']['adapter']['location']['localPath'])
        self.assertIn('chr1\t1\n', aliases.read_text())
        self.assertIn('chrX\tX\n', aliases.read_text())
        self.assertEqual(
            config['defaultSession']['views'][0]['displayedRegions'][0]['end'],
            195154279,
        )
        self.assertTrue(all(t['assemblyNames'] == ['GRCm39'] for t in config['tracks']))
        self.assertEqual(
            json.loads((self.output / 'build-summary.json').read_text())['assembly'],
            'GRCm39',
        )

    def test_gzip_directory_inputs_header_only_reports_and_fallback_session(self):
        rawdir, extdir, reportdir = (
            self.directory / name for name in ('rawdir', 'extdir', 'reports')
        )
        for path in (rawdir, extdir, reportdir):
            path.mkdir()
        for directory, source in [(rawdir, self.raw), (extdir, self.extended)]:
            with gzip.open(directory / 'annotation.gff3.gz', 'wt') as handle:
                handle.write(source.read_text())
        (reportdir / 'empty.csv').write_text(
            self.report.read_text().splitlines()[0] + '\n'
        )
        result = self.run_builder(
            '--reports', str(reportdir), raw=rawdir, extended=extdir
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        config = json.loads((self.output / 'config.json').read_text())
        self.assertIn('no overlap markers', config['defaultSession']['name'])
        self.assertEqual(len(config['tracks']), 10)

    def test_extra_datacheck_track_without_overlap_reports(self):
        marker = self.directory / 'other check.gff3'
        marker.write_text(
            '1\tcheck\tsequence_feature\t200\t205\t.\t+\t.\tID=issue1;Name=Changed%20exon;\n'
        )
        result = self.run_builder(
            '--extra-track', 'Internal exon changes=' + str(marker)
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        config = json.loads((self.output / 'config.json').read_text())
        self.assertEqual(len(config['tracks']), 5)
        self.assertEqual(config['tracks'][-1]['name'], 'Internal exon changes')
        self.assertIn('Internal exon changes', config['defaultSession']['name'])
        records = subprocess.check_output(
            ['tabix', str(self.output / 'data/extra_1.gff3.gz'), '1:1-500'], text=True
        )
        self.assertEqual(records, marker.read_text())

    def test_mixed_comparison_goes_to_both_tracks(self):
        with self.report.open(newline='') as handle:
            reader = csv.DictReader(handle)
            fields = reader.fieldnames
            row = next(reader)
        row.update(
            status='new_cds_overlap', failing_intervals='190-199', failing_bases='10'
        )
        with self.report.open('w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerow(row)
        result = self.run_builder('--reports', str(self.report))
        self.assertEqual(result.returncode, 0, result.stderr)
        for track, interval in [
            ('same_gene_failures', '\t190\t199\t'),
            ('shared_site_allowed', '\t160\t180\t'),
        ]:
            record = subprocess.check_output(
                ['tabix', str(self.output / ('data/' + track + '.gff3.gz')), '1:1-500'],
                text=True,
            )
            self.assertIn(interval, record)
            self.assertIn('comparison_status=new_cds_overlap', record)

    def test_current_annotation_allowance_has_its_own_track(self):
        with self.report.open(newline='') as handle:
            reader = csv.DictReader(handle)
            fields, row = reader.fieldnames, next(reader)
        row.update(
            status='allowed_current_annotation',
            annotation_allowed_intervals=row['allowed_intervals'],
            splice_site_allowed_intervals='',
            annotation_extension_ends='five_prime',
            shared_splice_sites='',
        )
        with self.report.open('w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerow(row)
        result = self.run_builder('--reports', str(self.report))
        self.assertEqual(result.returncode, 0, result.stderr)
        config = json.loads((self.output / 'config.json').read_text())
        self.assertIn('current-annotation extension', config['defaultSession']['name'])
        record = subprocess.check_output(
            ['tabix', str(self.output / 'data/annotation_allowed.gff3.gz'), '1:1-500'],
            text=True,
        )
        self.assertIn('\t160\t180\t', record)
        self.assertIn('annotation_extension_ends=five_prime', record)
        record = subprocess.check_output(
            ['tabix', str(self.output / 'data/shared_site_allowed.gff3.gz'), '1:1-500'],
            text=True,
        )
        self.assertEqual(record, '')

    def test_cross_gene_cds_routing_and_leap_audit_track(self):
        with self.report.open(newline='') as handle:
            reader = csv.DictReader(handle)
            fields = reader.fieldnames
            row = next(reader)
        row.update(
            target_gene='geneB',
            gene_relationship='cross_gene',
            status='new_cross_gene_cds_overlap',
            failing_intervals='160-170',
            splice_site_allowed_intervals='',
            allowed_intervals='',
            check_rule='full_length_cds',
            target_is_partial='false',
        )
        audit = dict(
            row,
            status='allowed_leap_overlap',
            failing_intervals='',
            allowed_intervals='171-180',
            leap_allowed_intervals='171-180',
            check_rule='outside_leap_protected_regions',
        )
        with self.report.open('w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerows([row, audit])
        result = self.run_builder('--reports', str(self.report))
        self.assertEqual(result.returncode, 0, result.stderr)
        cross = subprocess.check_output(
            ['tabix', str(self.output / 'data/cross_gene_failures.gff3.gz'), '1:1-500'],
            text=True,
        )
        self.assertIn('check_rule=full_length_cds', cross)
        self.assertIn('gene_relationship=cross_gene', cross)
        same = subprocess.check_output(
            ['tabix', str(self.output / 'data/same_gene_failures.gff3.gz'), '1:1-500'],
            text=True,
        )
        self.assertEqual(same, '')
        allowed = subprocess.check_output(
            ['tabix', str(self.output / 'data/leap_allowed.gff3.gz'), '1:1-500'],
            text=True,
        )
        self.assertIn('\t171\t180\t', allowed)
        self.assertIn('check_rule=outside_leap_protected_regions', allowed)

    def test_grouped_annotation_allowance_has_its_own_track(self):
        with self.report.open(newline='') as handle:
            row = next(csv.DictReader(handle))
        row.update(status='allowed_grouped_annotation',failing_intervals='',failing_bases='0',
                   allowed_intervals='160-180',allowed_bases='21',annotation_allowed_intervals='',
                   splice_site_allowed_intervals='',leap_allowed_intervals='',shared_splice_sites='',
                   grouped_annotation_allowed_intervals='160-180',
                   grouped_annotation_boundaries='five_prime:160',grouped_annotation_donors='donor.1')
        with self.report.open('w',newline='') as handle:
            writer = csv.DictWriter(handle,fieldnames=list(row));writer.writeheader();writer.writerow(row)
        result = self.run_builder('--reports',str(self.report))
        self.assertEqual(result.returncode,0,result.stderr)
        marker = subprocess.check_output(['tabix',str(self.output/'data/grouped_annotation_allowed.gff3.gz'),'1:1-500'],text=True)
        self.assertIn('\t160\t180\t',marker)
        self.assertIn('grouped_annotation_boundaries=five_prime:160',marker)
        self.assertIn('grouped_annotation_donors=donor.1',marker)
        for name in ('same_gene_failures','annotation_allowed','shared_site_allowed'):
            self.assertEqual(subprocess.check_output(['tabix',str(self.output/('data/'+name+'.gff3.gz')),'1:1-500'],text=True),'')

    def test_legacy_report_without_reason_columns(self):
        with self.report.open(newline='') as handle:
            reader = csv.DictReader(handle)
            row = next(reader)
        for field in (
            'annotation_allowed_intervals',
            'splice_site_allowed_intervals',
            'annotation_extension_ends',
        ):
            row.pop(field)
        with self.report.open('w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(row))
            writer.writeheader()
            writer.writerow(row)
        result = self.run_builder('--reports', str(self.report))
        self.assertEqual(result.returncode, 0, result.stderr)
        record = subprocess.check_output(
            ['tabix', str(self.output / 'data/shared_site_allowed.gff3.gz'), '1:1-500'],
            text=True,
        )
        self.assertIn('\t160\t180\t', record)

    def test_reject_existing_output_and_clean_up_failed_builds(self):
        self.output.mkdir()
        sentinel = self.output / 'keep.txt'
        sentinel.write_text('keep')
        result = self.run_builder()
        self.assertEqual(result.returncode, 2)
        self.assertEqual(sentinel.read_text(), 'keep')
        sentinel.unlink()
        self.output.rmdir()
        self.report.write_text('invalid,header\n')
        result = self.run_builder('--reports', str(self.report))
        self.assertEqual(result.returncode, 2)
        self.assertIn('missing report columns', result.stderr)
        self.assertFalse(self.output.exists())
        self.assertFalse(list(self.directory.glob('.review bundle-build-*')))

    def test_reject_wrong_assembly_and_missing_inputs(self):
        self.raw.write_text(self.raw.read_text().replace('1\ttest', 'chr1\ttest'))
        result = self.run_builder()
        self.assertEqual(result.returncode, 2)
        self.assertIn('chromosome names', result.stderr)
        self.assertFalse(self.output.exists())
        result = self.run_builder(raw=self.directory / 'absent')
        self.assertEqual(result.returncode, 2)
        self.assertIn('No matching input files', result.stderr)


if __name__ == '__main__':
    unittest.main()
