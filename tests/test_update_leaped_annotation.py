"""Exercise the Perl uploader without loading Otter or opening a database."""

import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / 'assets/update_LEAPed_annotation.pl'


@unittest.skipUnless(shutil.which('perl'), 'requires Perl')
class CanonicalExtendedUploadTests(unittest.TestCase):
    def test_planning_and_actual_editing_with_mock_database_objects(self):
        result = subprocess.run(
            ['perl', str(ROOT / 'tests/perl/test_canonical_extended.t'), str(SCRIPT)],
            capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertNotIn('not ok', result.stdout)

    def test_complete_script_compiles_with_otter_import_stubs(self):
        # Only dependency loading is stubbed for this compile-only check; no
        # runtime upload code executes under perl -c.
        with tempfile.TemporaryDirectory() as directory:
            modules = (
                'Bio::Otter::Lace::Defaults', 'Bio::Otter::Server::Config',
                'Bio::Vega::Exon', 'Bio::Otter::ServerAction::Script::Region',
                'Bio::Otter::Server::Support::Local',
            )
            for module in modules:
                path = Path(directory).joinpath(*module.split('::')).with_suffix('.pm')
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text('package ' + module + '; sub import {} 1;\n')
            result = subprocess.run(['perl', '-I', directory, '-c', str(SCRIPT)],
                                    capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn('syntax OK', result.stderr)


if __name__ == '__main__':
    unittest.main()
