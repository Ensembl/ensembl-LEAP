"""Regression tests for removing OTT models while retaining OTT exon IDs."""
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / 'bin/filter_out_ott_ids.py'


def feature(kind, gene, transcript=None, parent=None, exon=None):
    attrs = ['gene_id=' + gene]
    if transcript:
        attrs.append('transcript_id=' + transcript)
    if parent:
        attrs.append('Parent=' + parent)
    if exon:
        attrs.append('exon_id=' + exon)
    return '11\thavana\t{}\t100\t200\t.\t-\t.\t{}\n'.format(kind, ';'.join(attrs))


class OttFilterTests(unittest.TestCase):
    def test_removes_ott_models_and_preserves_ensembl_children(self):
        retained = [feature('gene', 'ENSG1'),
                    feature('transcript', 'ENSG1', 'ENST1', 'ENSG1'),
                    feature('exon', 'ENSG1', 'ENST1', 'ENST1', 'OTTHUME1'),
                    feature('CDS', 'ENSG1', 'ENST1', 'ENST1', 'OTTHUME1')]
        excluded = [feature('transcript', 'ENSG1', 'OTTHUMT00000505311.1', 'ENSG1'),
                    feature('exon', 'ENSG1', 'OTTHUMT00000505311.1', 'OTTHUMT00000505311.1', 'ENSE1'),
                    feature('CDS', 'ENSG1', 'OTTHUMT00000505311.1', 'OTTHUMT00000505311.1'),
                    feature('gene', 'OTTHUMG1'),
                    feature('transcript', 'OTTHUMG1', 'ENST2', 'OTTHUMG1'),
                    feature('exon', 'OTTHUMG1', 'ENST2', 'ENST2', 'ENSE2')]
        # Children preceding their transcript must still be removed.
        for records in (retained + excluded, list(reversed(retained + excluded))):
            with self.subTest(reversed=records[0] != retained[0]):
                with tempfile.TemporaryDirectory() as directory:
                    source, output = Path(directory) / 'input.gff3', Path(directory) / 'out.gff3'
                    source.write_text(''.join(records))
                    result = subprocess.run([sys.executable, str(SCRIPT), str(source), str(output)],
                                            capture_output=True, text=True)
                    self.assertEqual(result.returncode, 0, result.stderr)
                    self.assertEqual(output.read_text(), ''.join(record for record in records if record in retained))


if __name__ == '__main__':
    unittest.main()
