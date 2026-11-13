"""Exercise the live SQL path against a small synthetic Ensembl schema."""
import sys
import tempfile
import unittest
from contextlib import contextmanager
from pathlib import Path
from unittest.mock import patch
from urllib.parse import unquote

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'bin'))
try:
    from sqlalchemy import create_engine, text
except ImportError:
    create_engine = None

from gxf_dump.exporter import (
    export_main, fetch_annotation_rows, fetch_attribute_rows, fetch_chromosomes_from_db,
    gene_features_from_rows, records_to_features, rows_to_records,
)
from gxf_dump.regions import fetch_regions


@unittest.skipUnless(create_engine, 'Install runtime requirements to test database queries')
class AnnotationParityTests(unittest.TestCase):
    def setUp(self):
        self.engine = create_engine('sqlite://')
        self.connection = self.engine.connect()
        self.addCleanup(self.engine.dispose)
        self.addCleanup(self.connection.close)
        schema = '''
            CREATE TABLE coord_system (coord_system_id INTEGER, species_id INTEGER, name TEXT, version TEXT, attrib TEXT);
            CREATE TABLE meta (species_id INTEGER, meta_key TEXT, meta_value TEXT);
            CREATE TABLE seq_region (seq_region_id INTEGER, name TEXT, length INTEGER, coord_system_id INTEGER);
            CREATE TABLE attrib_type (attrib_type_id INTEGER, code TEXT);
            CREATE TABLE seq_region_attrib (seq_region_id INTEGER, attrib_type_id INTEGER, value TEXT);
            CREATE TABLE gene (gene_id INTEGER, stable_id TEXT, version INTEGER, source TEXT, biotype TEXT,
                is_current INTEGER, seq_region_id INTEGER, seq_region_start INTEGER, seq_region_end INTEGER,
                seq_region_strand INTEGER, display_xref_id INTEGER, analysis_id INTEGER);
            CREATE TABLE transcript (transcript_id INTEGER, gene_id INTEGER, stable_id TEXT, version INTEGER,
                source TEXT, biotype TEXT, seq_region_start INTEGER, seq_region_end INTEGER,
                seq_region_strand INTEGER, display_xref_id INTEGER);
            CREATE TABLE exon (exon_id INTEGER, stable_id TEXT, version INTEGER, seq_region_id INTEGER,
                seq_region_start INTEGER, seq_region_end INTEGER, phase INTEGER, end_phase INTEGER);
            CREATE TABLE exon_transcript (transcript_id INTEGER, exon_id INTEGER, rank INTEGER);
            CREATE TABLE translation (translation_id INTEGER, transcript_id INTEGER, stable_id TEXT, version INTEGER,
                start_exon_id INTEGER, end_exon_id INTEGER, seq_start INTEGER, seq_end INTEGER);
            CREATE TABLE xref (xref_id INTEGER, display_label TEXT);
            CREATE TABLE analysis (analysis_id INTEGER, logic_name TEXT);
            CREATE TABLE dna (seq_region_id INTEGER, sequence TEXT);
            CREATE TABLE gene_attrib (gene_id INTEGER, attrib_type_id INTEGER, value TEXT);
            CREATE TABLE transcript_attrib (transcript_id INTEGER, attrib_type_id INTEGER, value TEXT);
            INSERT INTO coord_system VALUES (77,1,'chromosome','GRCh38','default_version'),
                (88,1,'chromosome','GRCh37',''), (99,1,'scaffold','GRCh38','default_version');
            INSERT INTO meta VALUES (1,'assembly.default','GRCh38');
            INSERT INTO seq_region VALUES (10,'22',100,77), (20,'22',100,88),
                (30,'PATCH',100,99), (40,'SCAFFOLD',100,99), (50,'OLD_PATCH',100,88);
            INSERT INTO attrib_type VALUES (900,'toplevel'), (901,'non_ref'), (902,'remark'),
                (903,'artef_dupl'), (904,'hidden_remark'), (905,'havana_cv'), (906,'name'),
                (907,'cds_start_NF'), (908,'mRNA_end_NF');
            INSERT INTO seq_region_attrib VALUES (10,900,'1'), (20,900,'1'), (30,900,'1'),
                (40,900,'1'), (50,900,'1'), (30,901,'1'), (50,901,'1');
        '''
        for sql in schema.split(';'):
            if sql.strip():
                self.connection.exec_driver_sql(sql)

        @contextmanager
        def connect(_):
            yield self.connection
        self.patcher = patch('gxf_dump.exporter.sql_connection', connect)
        self.patcher.start()
        self.addCleanup(self.patcher.stop)

    def insert(self, table, **values):
        columns = ','.join(values)
        binds = ','.join(':' + key for key in values)
        self.connection.execute(text(f'INSERT INTO {table} ({columns}) VALUES ({binds})'), values)

    def model(self, n, region=10, gene_type='protein_coding', tx_type='lncRNA', with_tx=True):
        self.insert('gene', gene_id=n, stable_id=f'G{n}', version=1, source='havana', biotype=gene_type,
                    is_current=1, seq_region_id=region, seq_region_start=3, seq_region_end=12, seq_region_strand=1)
        if with_tx:
            self.insert('transcript', transcript_id=n, gene_id=n, stable_id=f'T{n}', version=1,
                        source='havana', biotype=tx_type, seq_region_start=3, seq_region_end=12, seq_region_strand=1)
            self.insert('exon', exon_id=n, stable_id=f'E{n}', version=1, seq_region_id=region,
                        seq_region_start=3, seq_region_end=12, phase=0, end_phase=0)
            self.insert('exon_transcript', transcript_id=n, exon_id=n, rank=1)

    def features(self, **options):
        rows = fetch_annotation_rows({}, **options)
        ga, ta = fetch_attribute_rows({}, {r['gene_db_id'] for r in rows},
                                      {r['transcript_db_id'] for r in rows if r['transcript_db_id'] is not None})
        records = rows_to_records(rows, gene_attrs=ga, transcript_attrs=ta)
        return records, records_to_features(records, 'gff3', genes=gene_features_from_rows(rows, ga, 'gff3'))

    def test_historical_current_genes_excluded_with_all(self):
        self.model(1)
        self.model(2, region=20)
        self.model(3, region=30)
        self.model(4, region=50)
        rows = fetch_annotation_rows({}, chromosome='22')
        self.assertEqual({r['gene_id'] for r in rows}, {'G1'})
        rows = fetch_annotation_rows({}, region_set='ALL')
        self.assertEqual({r['gene_id'] for r in rows}, {'G1', 'G3'})
        self.assertEqual({r['gene_id'] for r in fetch_annotation_rows({})}, {'G1'})

    def test_named_default_patch_matches_perl_slice_selection(self):
        self.model(3, region=30)
        self.assertEqual(fetch_annotation_rows({}, chromosome='PATCH')[0]['gene_id'], 'G3')
        with self.assertRaisesRegex(ValueError, 'default assembly'):
            fetch_annotation_rows({}, chromosome='OLD_PATCH')

    def test_discovery_switches_use_metadata_not_numeric_ids(self):
        self.assertEqual(fetch_chromosomes_from_db({}, True, False, False), ['22'])
        self.assertEqual(fetch_chromosomes_from_db({}, False, True, False), ['SCAFFOLD'])
        self.assertEqual(fetch_chromosomes_from_db({}, False, False, True), ['PATCH'])
        self.assertEqual(fetch_chromosomes_from_db({}, True, False, True), ['22', 'PATCH'])
        self.assertEqual(fetch_chromosomes_from_db({}, True, True, True), ['22', 'PATCH', 'SCAFFOLD'])

    def test_ambiguous_region_and_species_fail(self):
        self.insert('seq_region', seq_region_id=60, name='22', length=100, coord_system_id=77)
        self.insert('seq_region_attrib', seq_region_id=60, attrib_type_id=900, value='1')
        with self.assertRaisesRegex(ValueError, 'Ambiguous sequence-region'):
            fetch_regions(self.connection, '22')
        self.insert('coord_system', coord_system_id=100, species_id=2, name='chromosome',
                    version='other', attrib='default_version')
        with self.assertRaisesRegex(ValueError, 'species'):
            fetch_regions(self.connection, '22')

    def test_missing_and_ambiguous_default_fail(self):
        self.connection.exec_driver_sql('DELETE FROM meta')
        with self.assertRaisesRegex(ValueError, 'assembly.default'):
            fetch_regions(self.connection)
        self.insert('meta', species_id=1, meta_key='assembly.default', meta_value='GRCh38')
        self.insert('meta', species_id=1, meta_key='assembly.default', meta_value='GRCh37')
        with self.assertRaisesRegex(ValueError, 'assembly.default'):
            fetch_regions(self.connection)

    def test_reference_top_level_unversioned_system_and_non_top_exclusion(self):
        self.insert('coord_system', coord_system_id=101, species_id=1, name='contig',
                    version=None, attrib='default_version,sequence_level')
        self.insert('seq_region', seq_region_id=61, name='CONTIG', length=100, coord_system_id=101)
        self.insert('seq_region_attrib', seq_region_id=61, attrib_type_id=900, value='1')
        self.insert('seq_region', seq_region_id=62, name='COMPONENT', length=100, coord_system_id=101)
        regions = fetch_regions(self.connection, region_set='PRI')
        self.assertIn('CONTIG', {r['name'] for r in regions})
        self.assertNotIn('COMPONENT', {r['name'] for r in regions})

    def test_other_species_default_assembly(self):
        self.connection.exec_driver_sql("UPDATE coord_system SET version='GRCm39' WHERE version='GRCh38'")
        self.connection.exec_driver_sql("UPDATE meta SET meta_value='GRCm39'")
        self.model(1)
        self.assertEqual(fetch_annotation_rows({}, chromosome='22')[0]['gene_id'], 'G1')

    def test_region_presets_and_explicit_lists(self):
        self.insert('seq_region', seq_region_id=70, name='HAP', length=100, coord_system_id=77)
        self.insert('seq_region_attrib', seq_region_id=70, attrib_type_id=900, value='1')
        self.insert('seq_region_attrib', seq_region_id=70, attrib_type_id=901, value='1')
        for preset, expected in [(None, ['22']), ('CHR', ['22']),
                                 ('PRI', ['22', 'SCAFFOLD']),
                                 ('ALL', ['22', 'HAP', 'PATCH', 'SCAFFOLD'])]:
            with self.subTest(preset=preset):
                self.assertEqual([r['name'] for r in fetch_regions(self.connection, region_set=preset)], expected)
        self.assertEqual([r['name'] for r in fetch_regions(self.connection, ' PATCH,22, PATCH ')], ['22', 'PATCH'])
        with self.assertRaisesRegex(ValueError, 'MISSING'):
            fetch_regions(self.connection, '22,MISSING')
        with self.assertRaisesRegex(ValueError, 'not both'):
            fetch_regions(self.connection, '22', 'ALL')

    def test_nonhuman_chromosome_names_use_metadata(self):
        self.connection.exec_driver_sql("UPDATE seq_region SET name='LG_Alpha' WHERE seq_region_id=10")
        self.connection.exec_driver_sql("UPDATE coord_system SET version='FishAssembly' WHERE version='GRCh38'")
        self.connection.exec_driver_sql("UPDATE meta SET meta_value='FishAssembly'")
        self.assertEqual([r['name'] for r in fetch_regions(self.connection)], ['LG_Alpha'])
        self.assertEqual([r['name'] for r in fetch_regions(self.connection, 'LG_Alpha,PATCH')], ['LG_Alpha', 'PATCH'])

    def test_combined_live_export_equals_single_exports(self):
        from collections import Counter
        for n, region in [(1, 10), (2, 20), (3, 30), (4, 40)]:
            self.model(n, region=region)
        with tempfile.TemporaryDirectory() as tmp, patch('gxf_dump.exporter.load_secrets', return_value={}):
            output = Path(tmp) / 'result.gff3'
            base = ['--db-secrets', 'unused', '--output', str(output)]
            singles = []
            for name in ['22', 'PATCH', 'SCAFFOLD']:
                export_main(base + ['--chromosome', name])
                singles.extend(line for line in output.read_text().splitlines() if not line.startswith('#'))
            export_main(base + ['--chromosomes', 'SCAFFOLD,22,PATCH,22'])
            content = output.read_text()
            self.assertEqual(content.count('##gff-version 3'), 1)
            self.assertEqual(content.count('##sequence-region'), 3)
            self.assertEqual(content.count('\ttranscript\t'), 3)
            self.assertEqual(content.count('\texon\t'), 3)
            self.assertEqual(Counter(line for line in content.splitlines() if not line.startswith('#')), Counter(singles))
            export_main(base + ['--region-set', 'ALL'])
            self.assertEqual(output.read_text(), content)
            export_main(base)
            self.assertEqual(output.read_text().count('\tgene\t'), 1)
            self.assertNotIn('PATCH', output.read_text())
            self.assertNotIn('SCAFFOLD', output.read_text())
            export_main(base + ['--region-set', 'PRI'])
            self.assertEqual(output.read_text().count('\tgene\t'), 2)
            gtf_singles = []
            for name in ['22', 'PATCH']:
                export_main(base + ['--chromosome', name, '--format', 'gtf'])
                gtf_singles.extend(output.read_text().splitlines())
            export_main(base + ['--chromosome', '22,PATCH', '--format', 'gtf'])
            self.assertEqual(Counter(output.read_text().splitlines()), Counter(gtf_singles))
            self.assertNotIn('##', output.read_text())

    def test_multiple_region_bounds_and_unknown_names_preserve_output(self):
        with tempfile.TemporaryDirectory() as tmp, patch('gxf_dump.exporter.load_secrets', return_value={}):
            output = Path(tmp) / 'result.gff3'
            output.write_text('existing output')
            base = ['--db-secrets', 'unused', '--output', str(output)]
            for selection in [['--chromosome', '22,PATCH'], ['--region-set', 'PRI']]:
                with self.assertRaises(SystemExit):
                    export_main(base + selection + ['--start', '1'])
                self.assertEqual(output.read_text(), 'existing output')
            with self.assertRaises(SystemExit):
                export_main(base + ['--chromosome', '22,MISSING'])
            self.assertEqual(output.read_text(), 'existing output')
            export_main(base + ['--region-set', 'CHR', '--start', '2', '--end', '15'])
            self.assertIn('##sequence-region 22 2 15', output.read_text())

    def test_artifacts_vega_and_parent_duplication_are_independent(self):
        self.model(1, gene_type='artifact', tx_type='artifact')
        self.model(2, tx_type='artifact')
        self.model(3, gene_type='artifact', tx_type='artifact')
        self.model(4)
        self.model(5)
        self.insert('gene_attrib', gene_id=3, attrib_type_id=903, value='0')
        self.insert('transcript_attrib', transcript_id=4, attrib_type_id=902, value='not for VEGA')
        self.insert('gene_attrib', gene_id=5, attrib_type_id=902, value='not for VEGA')
        for artifacts, vega, expected in [(False, False, {'T3','T4','T5'}),
                                          (False, True, {'T3'}),
                                          (True, False, {'T1','T2','T3','T4','T5'}),
                                          (True, True, {'T1','T2','T3'})]:
            with self.subTest(artifacts=artifacts, vega=vega):
                records, _ = self.features(allow_artifacts=artifacts, vega_only=vega)
                self.assertEqual({r.transcript_stable_id for r in records}, expected)

    def test_gene_only_records_and_traversal_without_placeholder_transcripts(self):
        self.model(1, tx_type='artifact')
        self.model(2)
        self.model(3, with_tx=False)
        self.model(4)
        self.insert('transcript_attrib', transcript_id=4, attrib_type_id=902, value='not for VEGA')
        _, features = self.features(vega_only=True)
        self.assertEqual([(f.feature_type, f.attrs['gene_id']) for f in features],
                         [('gene','G1.1'), ('gene','G2.1'), ('transcript','G2.1'),
                          ('exon','G2.1'), ('gene','G3.1'), ('gene','G4.1')])
        rows = fetch_annotation_rows({}, transcript_ids={'T1'})
        self.assertEqual({r['gene_id'] for r in rows}, {'G1'})
        self.assertIsNone(rows[0]['transcript_db_id'])

    def test_remark_sources_allowlist_and_nonzero_flags(self):
        self.model(1)
        for attr, value in [(902,'CAGE supported TSS'), (904,'RNA-Seq supported only'),
                            (905,"alternative 5' UTR"), (905,'arbitrary private remark'),
                            (904,'No organism supported'), (907,'0'), (908,'1')]:
            self.insert('transcript_attrib', transcript_id=1, attrib_type_id=attr, value=value)
        records, _ = self.features()
        self.assertEqual(set(records[0].tags), {'CAGE_supported_TSS', 'RNA_Seq_supported_only',
                                              'alternative_5_UTR', 'not_organism_supported', 'mRNA_end_NF'})
        self.assertFalse(records[0].cds_start_nf)
        self.assertTrue(records[0].mrna_end_nf)

    def test_names_use_display_or_external_name_not_name_attribute(self):
        self.model(1)
        self.insert('transcript_attrib', transcript_id=1, attrib_type_id=906, value='AA-C720C749L.1-012')
        records, _ = self.features()
        self.assertEqual(records[0].transcript_name, 'T1')
        self.insert('xref', xref_id=1, display_label='AA-C720C749L.1-012')
        self.connection.exec_driver_sql('UPDATE transcript SET display_xref_id=1')
        records, _ = self.features()
        self.assertEqual(records[0].transcript_name, 'AA-C720C749L.1-012')
        rows = fetch_annotation_rows({})
        rows[0]['transcript_display_name'] = None
        rows[0]['transcript_external_name'] = 'external'
        self.assertEqual(rows_to_records(rows)[0].transcript_name, 'external')

    def test_gene_name_attribute_propagates_to_all_features_in_both_formats(self):
        self.model(1, tx_type='protein_coding')
        self.insert('gene_attrib', gene_id=1, attrib_type_id=906, value='PRR5')
        self.insert('transcript_attrib', transcript_id=1, attrib_type_id=906, value='private-transcript-name')
        self.insert('transcript_attrib', transcript_id=1, attrib_type_id=905, value='Ensembl_canonical_extended')
        self.insert('translation', translation_id=1, transcript_id=1, stable_id='P1', version=1,
                    start_exon_id=1, end_exon_id=1, seq_start=3, seq_end=8)
        rows = fetch_annotation_rows({})
        ga, ta = fetch_attribute_rows({}, {1}, {1})
        records = rows_to_records(rows, gene_attrs=ga, transcript_attrs=ta)
        self.assertEqual(records[0].gene_name, 'PRR5')
        self.assertEqual(records[0].transcript_name, 'T1')
        for fmt in ('gff3', 'gtf'):
            for genes in (None, gene_features_from_rows(rows, ga, fmt)):
                features = records_to_features(records, fmt, genes=genes)
                self.assertTrue({'gene', 'transcript', 'exon', 'CDS'} <= {f.feature_type for f in features})
                self.assertEqual({f.attrs['gene_name'] for f in features}, {'PRR5'})
                self.assertTrue(all('Ensembl_canonical_extended' in f.attrs['tag'].split(',')
                                    for f in features if f.feature_type != 'gene'))

    def test_gene_name_precedence_and_gene_only_output(self):
        self.model(1, tx_type='artifact')
        self.insert('gene_attrib', gene_id=1, attrib_type_id=906, value='PRR5')
        rows = fetch_annotation_rows({})
        ga, _ = fetch_attribute_rows({}, {1}, set())
        for display, external, attrs, expected in [
            ('display', 'external', ga, 'display'),
            (None, 'external', ga, 'external'),
            (None, None, ga, 'PRR5'),
            (None, None, {}, 'G1'),
        ]:
            with self.subTest(expected=expected):
                rows[0]['gene_display_name'] = display
                rows[0]['gene_external_name'] = external
                genes = gene_features_from_rows(rows, attrs, 'gff3')
                features = records_to_features(rows_to_records(rows), 'gff3', genes=genes)
                self.assertEqual(len(features), 1)
                self.assertEqual(features[0].attrs['gene_name'], expected)

    def test_canonical_extended_tag_is_retained_from_each_remark_source(self):
        self.model(1)
        for attr_id in (902, 904, 905):
            with self.subTest(attr_id=attr_id):
                self.connection.exec_driver_sql('DELETE FROM transcript_attrib')
                self.insert('transcript_attrib', transcript_id=1, attrib_type_id=attr_id,
                            value='Ensembl_canonical_extended')
                records, features = self.features()
                self.assertIn('Ensembl_canonical_extended', records[0].tags)
                self.assertTrue(all(f.attrs['tag'] == 'Ensembl_canonical_extended'
                                    for f in features if f.feature_type != 'gene'))

    def test_utr_biotype_rule_does_not_change_cds(self):
        self.model(1)
        self.insert('translation', translation_id=1, transcript_id=1, stable_id='P1', version=1,
                    start_exon_id=1, end_exon_id=1, seq_start=3, seq_end=8)
        rows = fetch_annotation_rows({})
        for biotype, expect_utrs in [('ig_gene',False), ('IG_V_gene',True), ('TR_V_gene',True),
                                     ('protein_coding',True), ('lncRNA',False)]:
            rows[0]['transcript_biotype'] = biotype
            features = records_to_features(rows_to_records(rows), 'gff3')
            self.assertEqual(any('UTR' in f.feature_type for f in features), expect_utrs)
            self.assertEqual([(f.start,f.end) for f in features if f.feature_type == 'CDS'], [(5,10)])

    def test_hidden_remarks_live_export_opt_in_and_ownership(self):
        self.model(1)
        self.model(2, with_tx=False)
        self.model(3)
        paragraphs = ['z last', '  First, paragraph; = & % "quoted" \\ path\tline\r\nnext \u03b1  ',
                      'z last', '', None, '   ']
        with tempfile.TemporaryDirectory() as tmp, patch('gxf_dump.exporter.load_secrets', return_value={}):
            output = Path(tmp) / 'result.gff3'
            base = ['--db-secrets', 'unused', '--output', str(output)]
            defaults = {}
            for fmt in ['gff3', 'gtf']:
                export_main(base + ['--format', fmt])
                defaults[fmt] = output.read_text()
            for value in paragraphs:
                # Synthetic attribute ID 904 deliberately differs from the usual 123.
                self.insert('gene_attrib', gene_id=1, attrib_type_id=904, value=value)
                self.insert('transcript_attrib', transcript_id=1, attrib_type_id=904, value=value)
            self.insert('gene_attrib', gene_id=2, attrib_type_id=904, value='gene without transcripts')
            self.insert('transcript_attrib', transcript_id=3, attrib_type_id=904, value=None)
            self.insert('transcript_attrib', transcript_id=3, attrib_type_id=904, value='')
            for fmt in ['gff3', 'gtf']:
                export_main(base + ['--format', fmt])
                self.assertEqual(output.read_text(), defaults[fmt])
                export_main(base + ['--format', fmt, '--include-hidden-remarks'])
                content = output.read_text()
                lines = [line.split('\t') for line in content.splitlines() if not line.startswith('#')]
                self.assertTrue(all(len(fields) == 9 for fields in lines))
                self.assertEqual(len(lines), len([line for line in defaults[fmt].splitlines() if not line.startswith('#')]))
                owners = [fields for fields in lines if 'hidden_remark' in fields[8]]
                self.assertEqual([fields[2] for fields in owners], ['gene', 'transcript', 'gene'])
                self.assertIn('G1.1', owners[0][8])
                self.assertIn('T1.1', owners[1][8])
                self.assertIn('G2.1', owners[2][8])
                if fmt == 'gff3':
                    for fields in owners[:2]:
                        encoded = next(part.split('=', 1)[1] for part in fields[8].split(';')
                                       if part.startswith('hidden_remark='))
                        self.assertEqual([unquote(value) for value in encoded.split(',')],
                                         sorted(value for value in paragraphs if value not in (None, '')))
                else:
                    self.assertEqual(owners[0][8].count('hidden_remark "'), 4)
                    self.assertEqual(owners[1][8].count('hidden_remark "'), 4)

    def test_hidden_remarks_do_not_change_child_features_tags_or_codons(self):
        self.model(1, tx_type='protein_coding')
        self.insert('translation', translation_id=1, transcript_id=1, stable_id='P1', version=1,
                    start_exon_id=1, end_exon_id=1, seq_start=3, seq_end=8)
        self.insert('gene_attrib', gene_id=1, attrib_type_id=904, value='gene note')
        self.insert('transcript_attrib', transcript_id=1, attrib_type_id=904, value='CAGE supported TSS')
        rows = fetch_annotation_rows({})
        ga, ta = fetch_attribute_rows({}, {1}, {1})
        records = rows_to_records(rows, gene_attrs=ga, transcript_attrs=ta, exon_sequences={'1': 'AAATGTAACC'})
        for fmt in ['gff3', 'gtf']:
            baseline = records_to_features(records, fmt, genes=gene_features_from_rows(rows, ga, fmt))
            enabled = records_to_features(
                records, fmt, genes=gene_features_from_rows(rows, ga, fmt, include_hidden_remarks=True),
                include_hidden_remarks=True)
            self.assertEqual(len(baseline), len(enabled))
            expected_types = {'exon', 'CDS', 'start_codon', 'stop_codon'}
            expected_types.update({'UTR'} if fmt == 'gtf' else {'five_prime_UTR', 'three_prime_UTR'})
            self.assertTrue(expected_types <= {feature.feature_type for feature in enabled})
            for original, feature in zip(baseline, enabled):
                if feature.feature_type in {'gene', 'transcript'}:
                    self.assertIn('hidden_remark', feature.attrs)
                    feature.attrs.pop('hidden_remark')
                else:
                    self.assertNotIn('hidden_remark', feature.attrs)
                self.assertEqual(feature, original)

    def test_hidden_remarks_preserved_on_gene_with_filtered_transcripts(self):
        self.model(1, tx_type='artifact')
        self.insert('gene_attrib', gene_id=1, attrib_type_id=904, value='kept gene note')
        self.insert('transcript_attrib', transcript_id=1, attrib_type_id=904, value='excluded transcript note')
        with tempfile.TemporaryDirectory() as tmp, patch('gxf_dump.exporter.load_secrets', return_value={}):
            output = Path(tmp) / 'result.gff3'
            export_main(['--db-secrets', 'unused', '--include-hidden-remarks', '--output', str(output)])
            self.assertIn('hidden_remark=kept gene note', output.read_text())
            self.assertNotIn('excluded', output.read_text())
            self.assertNotIn('\ttranscript\t', output.read_text())

    def test_live_export_headers_bounds_and_gene_only_gtf(self):
        self.model(1, tx_type='artifact')
        with tempfile.TemporaryDirectory() as tmp, patch('gxf_dump.exporter.load_secrets', return_value={}):
            output = Path(tmp) / 'result.gff3'
            args = ['--db-secrets','unused','--chromosome','22','--output',str(output), '--add-chr']
            export_main(args)
            content = output.read_text()
            self.assertTrue(content.startswith('##gff-version 3\n##sequence-region chr22 1 100\n'))
            self.assertEqual(content.count('\tgene\t'), 1)
            self.assertNotIn('\ttranscript\t', content)
            export_main(args + ['--start','2','--end','15'])
            self.assertIn('##sequence-region chr22 2 15', output.read_text())
            export_main(args + ['--format','gtf'])
            self.assertNotIn('##', output.read_text())
            self.assertIn('gene_id "G1.1";', output.read_text())

    def test_empty_named_region_keeps_header_and_invalid_bounds_preserve_output(self):
        with tempfile.TemporaryDirectory() as tmp, patch('gxf_dump.exporter.load_secrets', return_value={}):
            output = Path(tmp) / 'empty.gff3'
            args = ['--db-secrets','unused','--chromosome','22','--output',str(output)]
            export_main(args)
            expected = '##gff-version 3\n##sequence-region 22 1 100\n'
            self.assertEqual(output.read_text(), expected)
            with self.assertRaises(SystemExit):
                export_main(args + ['--start','0'])
            self.assertEqual(output.read_text(), expected)


    def test_leap_mane_remarks_and_selection_tags(self):
        cases = [('MANE_select', 'MANE_Select'), ('MANE_plus_clinical', 'MANE_Plus_Clinical'),
                 ('upstream atg', 'upstream_ATG'), ('cage supported tss', 'CAGE_supported_TSS'),
                 ('rp supported tis', 'RP_supported_TIS'), ('ensembl_canonical_extended', 'Ensembl_canonical_extended')]
        for n, (remark, expected) in enumerate(cases, 1):
            self.model(n)
            self.insert('transcript_attrib', transcript_id=n, attrib_type_id=902, value=remark)
        records, features = self.features()
        for n, (remark, expected) in enumerate(cases, 1):
            record = next(r for r in records if r.transcript_stable_id == f'T{n}')
            self.assertIn(expected, record.tags)
            children = [f for f in features if f.attrs.get('transcript_id') == f'T{n}.1']
            self.assertTrue(children)
            self.assertTrue(all(expected in f.attrs.get('tag','').split(',') for f in children))
        self.model(100)
        self.insert('transcript_attrib', transcript_id=100, attrib_type_id=904, value='MANE_select')
        records, _ = self.features()
        self.assertNotIn('MANE_Select', next(r for r in records if r.transcript_stable_id == 'T100').tags)
        self.insert('attrib_type', attrib_type_id=999, code='MANE_Select')
        self.insert('transcript_attrib', transcript_id=100, attrib_type_id=999, value='1')
        records, _ = self.features()
        self.assertIn('MANE_Select', next(r for r in records if r.transcript_stable_id == 'T100').tags)

    def test_leap_exclusion_queries_same_connection_default_assembly_and_remark_code(self):
        import get_current_annotation as prepare
        self.connection.exec_driver_sql('ALTER TABLE transcript ADD COLUMN is_current INTEGER DEFAULT 1')
        self.connection.exec_driver_sql('ALTER TABLE transcript ADD COLUMN seq_region_id INTEGER DEFAULT 10')
        self.model(1, tx_type='protein_coding')
        self.model(2)
        self.model(3, tx_type='protein_coding')
        self.model(4, region=20, tx_type='protein_coding')
        self.model(5, tx_type='protein_coding')
        self.connection.exec_driver_sql('UPDATE transcript SET seq_region_id=20 WHERE transcript_id=4')
        self.connection.exec_driver_sql('UPDATE transcript SET is_current=0 WHERE transcript_id=5')
        for n in (1,2,4,5):
            self.insert('transcript_attrib', transcript_id=n, attrib_type_id=902, value='CAGE supported TSS')
        self.insert('transcript_attrib', transcript_id=3, attrib_type_id=902, value='readthrough')
        readthrough, blacklist = prepare.exclusion_ids({}, '22')
        self.assertEqual(readthrough, {'G3'})
        self.assertEqual(blacklist, {'T1'})

if __name__ == '__main__':
    unittest.main()
