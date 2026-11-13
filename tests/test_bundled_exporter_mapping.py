import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'bin'))
from gxf_dump.mapping import MappingError, SequenceMappingPaths


class MappingTests(unittest.TestCase):
    def setUp(self):
        # Coordinate-system IDs and versions from the supplied database queries.
        self.systems = [
            dict(coord_system_id=4, species_id=1, name='chromosome', version='GRCh38', attrib='default_version'),
            dict(coord_system_id=3, species_id=1, name='scaffold', version='GRCh38', attrib='default_version'),
            dict(coord_system_id=2, species_id=1, name='clone', version=None, attrib='default_version'),
            dict(coord_system_id=1, species_id=1, name='contig', version=None, attrib='default_version,sequence_level'),
            dict(coord_system_id=12, species_id=1, name='chromosome', version='GRCh37', attrib=''),
            dict(coord_system_id=117, species_id=1, name='chromosome', version='NCBI36', attrib=''),
            dict(coord_system_id=1001, species_id=1, name='chromosome', version='NCBI35', attrib=''),
            dict(coord_system_id=10001, species_id=1, name='chromosome', version='NCBI34', attrib=''),
            dict(coord_system_id=10007, species_id=1, name='lrg', version=None, attrib='default_version'),
            dict(coord_system_id=10010, species_id=1, name='chromosome', version='OtterArchive', attrib=None),
        ]

    def paths(self, *declarations):
        return SequenceMappingPaths(self.systems, declarations)

    def test_reported_database_declarations(self):
        declarations = [
            'chromosome:GRCh37#chromosome:NCBI34', 'chromosome:GRCh37#chromosome:NCBI35',
            'chromosome:GRCh37#chromosome:NCBI36', 'chromosome:GRCh37#contig',
            'chromosome:GRCh38#chromosome:GRCh37', 'chromosome:GRCh38#chromosome:NCBI34',
            'chromosome:GRCh38#chromosome:NCBI35', 'chromosome:GRCh38#chromosome:NCBI36',
            'chromosome:GRCh38#contig', 'chromosome:GRCh38#contig|clone',
            'chromosome:GRCh38#lrg', 'chromosome:GRCh38|scaffold:GRCh38',
            'chromosome:NCBI34#contig', 'chromosome:NCBI35#contig',
            'chromosome:NCBI36#chromosome:NCBI35', 'chromosome:NCBI36#contig',
            'chromosome:OtterArchive#contig', 'clone|contig', 'lrg#contig',
            'scaffold:GRCh38#contig', 'scaffold:GRCh38#contig|clone',
        ]
        paths = self.paths(*declarations)
        for source in (4, 12, 1001, 10001, 117, 10010, 3, 10007):
            route = paths.to_sequence(source)
            self.assertEqual([(step.source, step.target, step.multiple) for step in route], [(source, 1, True)])
        self.assertEqual(paths.to_sequence(1), ())

    def test_default_and_explicit_versions(self):
        paths = self.paths('chromosome#contig', 'chromosome:GRCh37#contig')
        self.assertEqual(paths.to_sequence(4)[0].source, 4)
        self.assertEqual(paths.to_sequence(12)[0].source, 12)

    def test_reverse_declared_path(self):
        paths = self.paths('contig|clone|chromosome:GRCh38')
        self.assertEqual([(s.source, s.target) for s in paths.to_sequence(4)], [(4, 2), (2, 1)])

    def test_unambiguous_two_step_inference(self):
        paths = self.paths('chromosome:GRCh38|scaffold:GRCh38', 'scaffold:GRCh38#contig')
        self.assertEqual([(s.source, s.target) for s in paths.to_sequence(4)], [(4, 3), (3, 1)])

    def test_ambiguous_two_step_inference_fails(self):
        paths = self.paths('chromosome:GRCh38|scaffold:GRCh38', 'scaffold:GRCh38#contig',
                           'chromosome:GRCh38|clone', 'clone|contig')
        with self.assertRaisesRegex(MappingError, 'Ambiguous.*chromosome:GRCh38'):
            paths.to_sequence(4)

    def test_no_unrestricted_graph_search(self):
        paths = self.paths('chromosome:GRCh38|scaffold:GRCh38', 'scaffold:GRCh38|clone', 'clone|contig')
        with self.assertRaisesRegex(MappingError, 'No declared sequence path'):
            paths.to_sequence(4)

    def test_explicit_path_wins_over_inferred_historical_paths(self):
        paths = self.paths('chromosome:GRCh38#contig', 'chromosome:GRCh38#chromosome:GRCh37',
                           'chromosome:GRCh37#contig')
        self.assertEqual(len(paths.to_sequence(4)), 1)

    def test_shortest_explicit_path_and_duplicates(self):
        paths = self.paths('chromosome:GRCh38|contig', 'chromosome:GRCh38|contig',
                           'chromosome:GRCh38|scaffold:GRCh38|contig')
        self.assertEqual(len(paths.to_sequence(4)), 1)

    def test_equal_length_explicit_paths_are_ambiguous(self):
        paths = self.paths('chromosome:GRCh38|scaffold:GRCh38|contig', 'chromosome:GRCh38|clone|contig')
        with self.assertRaisesRegex(MappingError, 'Ambiguous'):
            paths.to_sequence(4)

    def test_missing_and_unknown_paths(self):
        with self.assertRaisesRegex(MappingError, 'No declared'):
            self.paths().to_sequence(4)
        with self.assertRaisesRegex(MappingError, 'Unknown or ambiguous'):
            self.paths('chromosome:NotAnAssembly#contig')

    def test_unknown_default_fails_when_versions_are_ambiguous(self):
        self.systems[0]['attrib'] = ''
        with self.assertRaisesRegex(MappingError, 'ambiguous coordinate system'):
            self.paths('chromosome#contig')

    def test_zero_or_multiple_sequence_levels_fail(self):
        self.systems[0]['attrib'] = 'sequence_level'
        with self.assertRaisesRegex(MappingError, 'Expected one sequence_level'):
            self.paths()
        self.systems[0]['attrib'] = ''
        self.systems[3]['attrib'] = ''
        with self.assertRaisesRegex(MappingError, 'Expected one sequence_level'):
            self.paths()

    def test_malformed_and_cyclic_declarations(self):
        for declaration in ('contig', 'chromosome||contig', 'chromosome|clone|chromosome|contig'):
            with self.subTest(declaration=declaration), self.assertRaises(MappingError):
                self.paths(declaration)


if __name__ == '__main__':
    unittest.main()
