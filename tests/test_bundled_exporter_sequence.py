import sqlite3
import sys
import unittest
from contextlib import contextmanager
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'bin'))

from gxf_dump.sequence import SQLSequenceStore, SequenceError, SequenceResolver
from gxf_dump.exporter import (
    ExonRecord, attach_exon_sequences, fetch_exon_sequences,
    sequence_validation_available,
)


class SQLiteStore(SQLSequenceStore):
    def rows(self, sql, **params):
        return [dict(row) for row in self.connection.execute(sql, params)]


class SequenceTests(unittest.TestCase):
    def setUp(self):
        self.db = sqlite3.connect(':memory:')
        self.addCleanup(self.db.close)
        self.db.row_factory = sqlite3.Row
        self.db.executescript('''
            CREATE TABLE coord_system (coord_system_id INTEGER, species_id INTEGER,
                name TEXT, version TEXT, rank INTEGER, attrib TEXT);
            INSERT INTO coord_system VALUES
                (1, 1, 'chromosome', 'GRCh38', 1, 'default_version'),
                (2, 1, 'scaffold', 'GRCh38', 2, 'default_version'),
                (3, 1, 'contig', NULL, 3, 'default_version,sequence_level'),
                (4, 1, 'chromosome', 'GRCh37', 4, ''),
                (5, 1, 'chromosome', 'NCBI35', 5, ''),
                (6, 1, 'clone', NULL, 6, 'default_version');
            CREATE TABLE meta (meta_id INTEGER PRIMARY KEY, species_id INTEGER, meta_key TEXT, meta_value TEXT);
            INSERT INTO meta (species_id, meta_key, meta_value) VALUES
                (1, 'assembly.mapping', 'chromosome:GRCh38#contig'),
                (1, 'assembly.mapping', 'scaffold:GRCh38#contig'),
                (1, 'assembly.mapping', 'chromosome:GRCh37#contig'),
                (1, 'assembly.mapping', 'chromosome:NCBI35#contig');
            CREATE TABLE seq_region (seq_region_id INTEGER, name TEXT, length INTEGER, coord_system_id INTEGER);
            CREATE TABLE dna (seq_region_id INTEGER, sequence TEXT);
            CREATE TABLE assembly (asm_seq_region_id INTEGER, cmp_seq_region_id INTEGER,
                asm_start INTEGER, asm_end INTEGER, cmp_start INTEGER, cmp_end INTEGER, ori INTEGER);
            CREATE TABLE assembly_exception (seq_region_id INTEGER, seq_region_start INTEGER,
                seq_region_end INTEGER, exc_type TEXT, exc_seq_region_id INTEGER,
                exc_seq_region_start INTEGER, exc_seq_region_end INTEGER, ori INTEGER);
        ''')

    def region(self, id, length, dna=None, cs=1, name=None):
        self.db.execute('INSERT INTO seq_region VALUES (?, ?, ?, ?)', (id, name or str(id), length, cs))
        if dna is not None:
            self.db.execute('INSERT INTO dna VALUES (?, ?)', (id, dna))

    def assembly(self, source, target, start, end, target_start=1, ori=1):
        self.db.execute('INSERT INTO assembly VALUES (?, ?, ?, ?, ?, ?, ?)',
                        (source, target, start, end, target_start, target_start+end-start, ori))

    def exception(self, source, start, end, target, target_start, target_end, kind='PATCH_FIX', ori=1):
        self.db.execute('INSERT INTO assembly_exception VALUES (?, ?, ?, ?, ?, ?, ?, ?)',
                        (source, start, end, kind, target, target_start, target_end, ori))

    def resolver(self):
        return SequenceResolver(SQLiteStore(self.db))

    def mappings(self, *paths):
        self.db.execute('DELETE FROM meta')
        self.db.executemany("INSERT INTO meta (species_id, meta_key, meta_value) VALUES (1, 'assembly.mapping', ?)",
                            [(path,) for path in paths])

    def test_direct_and_truncated_dna(self):
        self.region(1, 8, 'aaccggtt', cs=3)
        self.assertEqual(self.resolver().fetch(1, 2, 6), 'ACCGG')
        self.assertEqual(self.resolver().fetch(1, -1, 2), 'NNAA')
        self.db.execute("UPDATE dna SET sequence='AAC'")
        with self.assertRaisesRegex(SequenceError, 'Truncated'):
            self.resolver().fetch(1, 1, 8)

    def test_recursive_reverse_orientation_and_gaps(self):
        self.mappings('chromosome:GRCh38|scaffold:GRCh38#contig')
        self.region(1, 12)
        self.region(2, 8, cs=2)
        self.region(3, 8, 'AACCGTTA', cs=3)
        self.assembly(1, 2, 3, 10, ori=-1)
        self.assembly(2, 3, 1, 8)
        self.assertEqual(self.resolver().fetch(1, 1, 12), 'NNTAACGGTTNN')

    def test_missing_component_dna_is_not_a_gap(self):
        self.region(1, 8)
        self.region(2, 8, cs=3)
        self.assembly(1, 2, 1, 8)
        with self.assertRaisesRegex(SequenceError, 'No DNA available'):
            self.resolver().fetch(1, 1, 8)

    def test_alternative_routes_are_not_concatenated(self):
        self.region(1, 8)
        self.region(2, 8, cs=2)
        self.region(3, 8, 'AACCGGTT', cs=3)
        self.assembly(1, 2, 1, 8)
        self.assembly(1, 3, 1, 8)
        self.assembly(2, 3, 1, 8)
        self.assertEqual(self.resolver().fetch(1, 1, 8), 'AACCGGTT')
        self.db.execute("INSERT INTO dna VALUES (2, 'TTTTTTTT')")
        self.assertEqual(self.resolver().fetch(1, 1, 8), 'AACCGGTT')

    def test_declared_recursive_route_ignores_other_direct_route(self):
        self.mappings('chromosome:GRCh38|scaffold:GRCh38#contig')
        self.region(1, 8)
        self.region(2, 8, cs=2)
        self.region(3, 8, cs=3)
        self.region(4, 8, 'AACCGGTT', cs=3)
        self.assembly(1, 3, 1, 8)
        self.assembly(1, 2, 1, 8)
        self.assembly(2, 4, 1, 8)
        self.assertEqual(self.resolver().fetch(1, 1, 8), 'AACCGGTT')

    def test_patch_flanks_and_boundary_crossing(self):
        self.region(1, 10)
        self.region(2, 12, 'AAAACCCCGGGG', cs=3)
        self.region(3, 2, 'TT', cs=3)
        self.exception(1, 5, 6, 2, 5, 8)
        self.assembly(1, 3, 5, 6)
        resolver = self.resolver()
        self.assertEqual(resolver.fetch(1, 1, 4), 'AAAA')
        self.assertEqual(resolver.fetch(1, 7, 10), 'GGGG')
        self.assertEqual(resolver.fetch(1, 3, 8), 'AATTGG')

    def test_multiple_patch_exceptions(self):
        self.region(1, 10)
        self.region(2, 12, 'AACCGGTTAACC', cs=3)
        self.region(3, 10, 'TTTTTTTTTT', cs=3)
        self.exception(1, 3, 4, 2, 3, 6)
        self.exception(1, 7, 8, 2, 9, 10)
        self.assembly(1, 3, 1, 10)
        self.assertEqual(self.resolver().fetch(1, 1, 10), 'AATTTTTTCC')

    def test_par_projects_interval_not_flanks(self):
        self.region(1, 8, 'AAAAAAAA', cs=3)
        self.region(2, 8, 'CCGTTGCC', cs=3)
        self.exception(1, 3, 6, 2, 3, 6, kind='PAR', ori=-1)
        self.assertEqual(self.resolver().fetch(1, 1, 8), 'AACAACAA')

    def test_cycle_overlap_and_invalid_flank_fail(self):
        self.mappings('chromosome:GRCh38|scaffold:GRCh38|chromosome:GRCh38|contig')
        self.region(1, 8)
        self.region(2, 8)
        self.assembly(1, 2, 1, 8)
        self.assembly(2, 1, 1, 8)
        with self.assertRaisesRegex(SequenceError, 'Cyclic'):
            self.resolver().fetch(1, 1, 8)
        self.assembly(1, 2, 2, 4)
        with self.assertRaisesRegex(SequenceError, 'Cyclic'):
            self.resolver().fetch(1, 1, 8)
        self.exception(1, 4, 5, 2, 3, 4)
        with self.assertRaisesRegex(SequenceError, 'inconsistent lengths'):
            self.resolver().fetch(1, 1, 2)

    def test_real_patch_coordinates_use_reference_and_169_offset(self):
        self.region(2006292501, 133275140, name='CHR_HG2063_PATCH')
        self.region(12, 133275309, cs=3)
        self.exception(2006292501, 86768210, 87083819, 12, 86768210, 87083988)
        store = SQLiteStore(self.db)
        requests = []

        def dna(region, start, end):
            requests.append((region, start, end))
            return 'A' * (end-start+1)

        store.dna = dna
        resolver = SequenceResolver(store)
        for start, end in [(86049674, 86049723), (85989400, 85989552),
                           (85983523, 85983670), (85955666, 85980430)]:
            self.assertEqual(len(resolver.fetch(2006292501, start, end)), end-start+1)
            self.assertEqual(requests[-1], (12, start, end))
        resolver.fetch(2006292501, 87083820, 87083829)
        self.assertEqual(requests[-1], (12, 87083989, 87083998))

    def test_exporter_orients_once_and_checks_full_exon_length(self):
        self.region(1, 12)
        self.region(2, 12, 'TTATTTCATCCC', cs=3)
        self.exception(1, 10, 12, 2, 10, 12)
        exon = ExonRecord('10', 'EX1.1', 1, 9, 1, seq_region_id='1')
        record = SimpleNamespace(transcript_id='TX1.1', transcript_strand='-',
            seq_start=1, seq_end=9, translation_start_exon_id='10',
            translation_end_exon_id='10', exons=[exon])

        @contextmanager
        def connection(secrets):
            yield self.db

        with patch('gxf_dump.exporter.sql_connection', connection), patch('gxf_dump.sequence.SQLSequenceStore', SQLiteStore):
            sequences = fetch_exon_sequences({}, [record])
            self.assertEqual(sequences['10'], 'TTATTTCAT')
            attach_exon_sequences([record], sequences)
            self.assertEqual(exon.sequence, 'ATGAAATAA')
            self.assertTrue(sequence_validation_available(record))
            exon.sequence = exon.sequence[:-1]
            self.assertFalse(sequence_validation_available(record))
            self.db.execute('DELETE FROM dna')
            with self.assertRaisesRegex(ValueError, r'TX1.1, exon EX1.1, region 1:1-9.*No DNA'):
                fetch_exon_sequences({}, [record])

    def test_terminal_patch_flank_padding(self):
        self.region(1, 12)
        self.region(2, 10, 'AACCGGTTAA', cs=3)
        self.exception(1, 3, 4, 2, 3, 4)
        self.assertEqual(self.resolver().fetch(1, 9, 12), 'AANN')

    def test_overlapping_components_rejected(self):
        self.mappings('chromosome:GRCh38|contig')
        self.region(1, 8)
        self.region(2, 8, 'AAAAAAAA', cs=3)
        self.assembly(1, 2, 1, 6)
        self.assembly(1, 2, 5, 8)
        with self.assertRaisesRegex(SequenceError, 'Overlapping components'):
            self.resolver().fetch(1, 1, 8)

    def test_invalid_component_bounds_are_not_padded(self):
        self.region(1, 8)
        self.region(2, 4, 'AAAA', cs=3)
        self.assembly(1, 2, 1, 8)
        with self.assertRaisesRegex(SequenceError, 'exceeds target'):
            self.resolver().fetch(1, 1, 8)

    def test_non_overlapping_route_does_not_conflict(self):
        self.region(1, 100)
        self.region(2, 8, 'AACCGGTT', cs=3)
        self.region(3, 8, cs=2)
        self.assembly(1, 2, 1, 8)
        self.assembly(1, 3, 50, 57)
        self.assertEqual(self.resolver().fetch(1, 1, 8), 'AACCGGTT')

    def test_gap_only_interval_on_declared_route_is_padded(self):
        self.region(1, 100)
        self.region(2, 8, 'AACCGGTT', cs=3)
        self.assembly(1, 2, 50, 57)
        self.assertEqual(self.resolver().fetch(1, 1, 8), 'N' * 8)

    def test_partial_alternative_compares_only_covered_bases(self):
        self.region(1, 8)
        self.region(2, 4, 'CCGG', cs=2)
        self.region(3, 8, 'AACCGGTT', cs=3)
        self.assembly(1, 2, 3, 6)
        self.assembly(1, 3, 1, 8)
        self.assertEqual(self.resolver().fetch(1, 1, 8), 'AACCGGTT')

    def test_declared_recursive_path_precedes_undeclared_direct_route(self):
        self.mappings('chromosome:GRCh38|scaffold:GRCh38#contig')
        self.region(1, 8)
        self.region(2, 8, cs=2)
        self.region(3, 4, 'CCGG', cs=3)
        self.region(4, 8, 'AACCGGTT', cs=3)
        self.assembly(1, 3, 3, 6)
        self.assembly(1, 2, 1, 8)
        self.assembly(2, 4, 1, 8)
        self.assertEqual(self.resolver().fetch(1, 1, 8), 'AACCGGTT')

    def test_recursive_padding_is_not_mapped_coverage(self):
        self.mappings('chromosome:GRCh38#scaffold:GRCh38#contig')
        self.region(1, 8)
        self.region(2, 8, cs=2)
        self.region(3, 8, cs=2)
        self.region(4, 4, 'CCGG', cs=3)
        self.region(5, 8, 'AACCGGTT', cs=3)
        self.assembly(1, 2, 1, 8)
        self.assembly(1, 3, 1, 8)
        self.assembly(2, 4, 3, 6)
        self.assembly(3, 5, 1, 8)
        self.assertEqual(self.resolver().fetch(1, 1, 8), 'AACCGGTT')

    def test_recursive_coverage_is_reversed_with_sequence(self):
        self.mappings('chromosome:GRCh38#scaffold:GRCh38#contig')
        self.region(1, 8)
        self.region(2, 8, cs=2)
        self.region(3, 8, cs=2)
        self.region(4, 3, 'AAC', cs=3)
        self.region(5, 8, 'AACCGGTT', cs=3)
        self.assembly(1, 2, 1, 8, ori=-1)
        self.assembly(1, 3, 1, 8)
        self.assembly(2, 4, 1, 3)
        self.assembly(3, 5, 1, 8)
        self.assertEqual(self.resolver().fetch(1, 1, 8), 'AACCGGTT')

    def test_dna_backed_n_is_not_treated_as_padding(self):
        self.mappings('chromosome:GRCh38#contig')
        self.region(1, 8)
        self.region(2, 4, 'NCGG', cs=3)
        self.region(3, 8, 'AACCGGTT', cs=3)
        self.assembly(1, 2, 3, 6)
        self.assembly(1, 3, 1, 8)
        with self.assertRaisesRegex(SequenceError, r'first conflict at seq_region_id=1:3 \(C vs N\)'):
            self.resolver().fetch(1, 1, 8)

    def test_conflict_reports_coordinate_systems_components_and_intervals(self):
        self.region(1, 8)
        self.region(2, 4, 'CTGG', cs=3)
        self.region(3, 8, 'AACCGGTT', cs=3)
        self.assembly(1, 2, 3, 6)
        self.assembly(1, 3, 1, 8)
        with self.assertRaises(SequenceError) as raised:
            self.resolver().fetch(1, 1, 8)
        message = str(raised.exception)
        self.assertIn('first conflict at seq_region_id=1:4 (C vs T)', message)
        self.assertIn('chromosome:GRCh38 (coord_system_id=1) # contig (coord_system_id=3)', message)
        self.assertIn('3-6 -> component 2:1-4 (ori=1)', message)
        self.assertIn('1-8 -> component 3:1-8 (ori=1)', message)

    def test_partial_routes_are_not_spliced_together(self):
        self.region(1, 8)
        self.region(2, 4, 'TTTT', cs=2)
        self.region(3, 4, 'AAAA', cs=3)
        self.assembly(1, 2, 5, 8)
        self.assembly(1, 3, 1, 4)
        self.assertEqual(self.resolver().fetch(1, 1, 8), 'AAAANNNN')

    def test_missing_dna_in_partial_route_is_not_padded(self):
        self.region(1, 8)
        self.region(2, 4, cs=3)
        self.assembly(1, 2, 3, 6)
        with self.assertRaisesRegex(SequenceError, 'No DNA available'):
            self.resolver().fetch(1, 1, 8)

    def test_nested_conflict_is_not_hidden_by_another_valid_route(self):
        self.mappings('chromosome:GRCh38|scaffold:GRCh38#contig')
        self.region(1, 8)
        self.region(2, 8, cs=2)
        self.region(3, 8, 'AACCGGTT', cs=3)
        self.region(4, 8, 'TTTTTTTT', cs=3)
        self.assembly(1, 2, 1, 8)
        self.assembly(1, 3, 1, 8)
        self.assembly(2, 3, 1, 8)
        self.assembly(2, 4, 1, 8)
        with self.assertRaisesRegex(SequenceError, 'first conflict at seq_region_id=2:1'):
            self.resolver().fetch(1, 1, 8)

    def test_reported_patch_interval_with_non_overlapping_reference_route(self):
        # Synthetic DNA and exception placement, using the reported source/target
        # exon coordinates. This tests the failure shape, not live DB parity.
        self.region(2006446536, 130000000, name='CHR_HG2047_PATCH')
        self.region(131546, 130000005, name='12')
        self.region(2, 145, 'A' * 145, cs=3)
        self.region(3, 10, cs=2)
        self.exception(2006446536, 100, 109, 131546, 100, 114)
        self.assembly(131546, 2, 121669339, 121669483)
        self.assembly(131546, 3, 1, 10)
        self.assertEqual(self.resolver().fetch(2006446536, 121669334, 121669478), 'A' * 145)

    def test_server_assembly_versions_do_not_enter_sequence_route(self):
        self.db.execute('DELETE FROM coord_system')
        self.db.executemany('INSERT INTO coord_system VALUES (?, 1, ?, ?, ?, ?)', [
            (4, 'chromosome', 'GRCh38', 1, 'default_version'),
            (12, 'chromosome', 'GRCh37', 5, ''),
            (1001, 'chromosome', 'NCBI35', 7, ''),
            (1, 'contig', None, 4, 'default_version,sequence_level'),
        ])
        self.mappings('chromosome:GRCh38#contig', 'chromosome:GRCh37#contig',
                      'chromosome:NCBI35#contig', 'chromosome:GRCh38#chromosome:GRCh37',
                      'chromosome:GRCh37#chromosome:NCBI35')
        self.region(131550, 248956422, cs=4, name='1')
        self.region(292305, 249250621, cs=12, name='1')
        self.region(1002619423, 245522847, cs=1001, name='1')
        self.region(7023, 112159, 'A' * 112159, cs=1, name='AL162741.44')
        self.region(8000, 117, 'C' * 117, cs=1)
        self.assembly(131550, 8000, 1196600, 1196716)
        self.assembly(131550, 292305, 1196600, 1196716, target_start=1131980)
        self.assembly(292305, 7023, 1131980, 1132096, target_start=6712)
        self.assembly(292305, 1002619423, 1131980, 1132096, target_start=1171903)
        # A broken historical route must not affect a modern sequence request.
        store = SQLiteStore(self.db)
        with patch.object(store, 'region', wraps=store.region) as regions, patch.object(store, 'dna', wraps=store.dna) as dna:
            self.assertEqual(SequenceResolver(store).fetch(131550, 1196600, 1196716), 'C' * 117)
            self.assertEqual({call.args[0] for call in regions.call_args_list}, {131550, 8000})
            dna.assert_called_once_with(8000, 1, 117)
        # Explicit GRCh37 requests use GRCh37's own contigs, not NCBI35.
        self.assertEqual(self.resolver().fetch(292305, 1131980, 1132096), 'A' * 117)

    def test_reverse_assembly_traversal_and_orientation(self):
        self.mappings('clone|chromosome:GRCh38', 'clone|contig')
        self.region(1, 4)
        self.region(2, 8, cs=6)
        self.region(3, 8, 'AACGTTAA', cs=3)
        self.assembly(2, 1, 3, 6, ori=-1)
        self.assembly(2, 3, 1, 8)
        self.assertEqual(self.resolver().fetch(1, 1, 4), 'AACG')

    def test_explicit_multistep_path_is_not_reselected_at_intermediate_region(self):
        self.mappings('chromosome:GRCh38|scaffold:GRCh38|clone|contig',
                      'scaffold:GRCh38|contig')
        self.region(1, 4)
        self.region(2, 4, cs=2)
        self.region(3, 4, cs=6)
        self.region(4, 4, 'AAAA', cs=3)
        self.region(5, 4, 'TTTT', cs=3)
        self.assembly(1, 2, 1, 4)
        self.assembly(2, 3, 1, 4)
        self.assembly(3, 4, 1, 4)
        self.assembly(2, 5, 1, 4)
        self.assertEqual(self.resolver().fetch(1, 1, 4), 'AAAA')
        self.assertEqual(self.resolver().fetch(2, 1, 4), 'TTTT')

    def test_patch_projection_selects_destination_coordinate_system_path(self):
        self.region(1, 8)
        self.region(2, 8, cs=4)
        self.region(3, 8, 'AACCGGTT', cs=3)
        self.exception(1, 7, 8, 2, 7, 8)
        self.assembly(2, 3, 1, 8)
        self.assertEqual(self.resolver().fetch(1, 1, 6), 'AACCGG')

    def test_multiple_parts_of_same_contig_are_preserved(self):
        self.region(1, 10)
        self.region(2, 8, 'AACCGGTT', cs=3)
        self.assembly(1, 2, 1, 4)
        self.assembly(1, 2, 7, 10, target_start=5, ori=-1)
        self.assertEqual(self.resolver().fetch(1, 1, 10), 'AACCNNAACC')

    def test_consistent_overlapping_components_do_not_duplicate_bases(self):
        self.region(1, 8)
        self.region(2, 8, 'AACCGGTT', cs=3)
        self.region(3, 4, 'CCGG', cs=3)
        self.assembly(1, 2, 1, 8)
        self.assembly(1, 3, 3, 6)
        self.assertEqual(self.resolver().fetch(1, 1, 8), 'AACCGGTT')

    def test_missing_selected_route_does_not_fall_back(self):
        self.region(1, 8)
        self.region(2, 8, cs=4)
        self.region(3, 8, 'AACCGGTT', cs=3)
        self.mappings('chromosome:GRCh38#contig', 'chromosome:GRCh38#chromosome:GRCh37',
                      'chromosome:GRCh37#contig')
        self.assembly(1, 2, 1, 8)
        self.assembly(2, 3, 1, 8)
        with self.assertRaisesRegex(SequenceError, 'Selected path chromosome:GRCh38.*No assembly components available'):
            self.resolver().fetch(1, 1, 8)

    def test_metadata_is_scoped_to_species(self):
        self.db.execute("INSERT INTO coord_system VALUES (50, 2, 'contig', NULL, 1, 'sequence_level,default_version')")
        self.db.execute("INSERT INTO meta (species_id, meta_key, meta_value) VALUES (2, 'assembly.mapping', 'invalid|path')")
        self.region(1, 4)
        self.region(2, 4, 'ACGT', cs=3)
        self.assembly(1, 2, 1, 4)
        self.assertEqual(self.resolver().fetch(1, 1, 4), 'ACGT')


if __name__ == '__main__':
    unittest.main()
