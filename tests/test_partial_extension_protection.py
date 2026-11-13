"""Partials may restrict original candidates, never create or lengthen them."""
import csv
import gzip
import importlib.util
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from test_partial_transcript_processing import model, parse

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / 'bin/protect_partial_extensions.py'


class PartialProtectionTests(unittest.TestCase):
    def inputs(self, folder, partial='', full=None, strand='+', natural=''):
        annotation = folder / 'full.gff'
        annotation.write_text(model('source.1', strand=strand) if full is None else full)
        reference = folder / 'partial.gff'
        reference.write_text(partial)
        pairs = folder / 'natural.txt'
        pairs.write_text(natural)
        fe, le = folder / 'original_FE.bed', folder / 'original_LE.bed'
        if strand == '+':
            fe.write_text('1\t150\t250\t185/200\t1\t+\t.\t.\t3\tA\n')
            le.write_text('1\t400\t500\t470/450\t1\t+\t.\t.\t4\tA\n')
        else:
            fe.write_text('1\t400\t500\t470/450\t1\t-\t.\t.\t4\tA\n')
            le.write_text('1\t150\t250\t185/200\t1\t-\t.\t.\t3\tA\n')
        return [fe, le, annotation, reference, pairs]

    def protect(self, folder, inputs, ok=True):
        prefix = folder / 'protected'
        result = subprocess.run([sys.executable, str(SCRIPT), *map(str, inputs), str(prefix)],
                                capture_output=True, text=True)
        if ok:
            self.assertEqual(result.returncode, 0, result.stderr)
        paths = [folder / 'protected_FE_protected.bed', folder / 'protected_LE_protected.bed']
        report = folder / 'protected_partial_protection.tsv'
        rows = []
        if report.exists():
            with report.open() as handle:
                rows = list(csv.DictReader(handle, delimiter='\t'))
        return result, paths, rows

    def test_empty_partial_file_is_exact_noop(self):
        for strand in ['+', '-']:
            with self.subTest(strand=strand), tempfile.TemporaryDirectory() as tmp:
                folder = Path(tmp)
                inputs = self.inputs(folder, strand=strand)
                # Preserve stale longer alternatives and provenance when unchanged.
                inputs[0].write_bytes(inputs[0].read_bytes().replace(b'\n', b'\r\n'))
                _, outputs, audit = self.protect(folder, inputs)
                self.assertEqual(audit, [])
                for source, output in zip(inputs[:2], outputs):
                    self.assertEqual(source.read_bytes(), output.read_bytes())

    def test_cds_protection_selects_shorter_original_read_on_both_ends_and_strands(self):
        for strand in ['+', '-']:
            for kind in ['FE', 'LE']:
                with self.subTest(strand=strand, kind=kind), tempfile.TemporaryDirectory() as tmp:
                    folder = Path(tmp)
                    left = (kind == 'FE') == (strand == '+')
                    exons = ((100, 170), (300, 350)) if left else ((480, 520), (600, 650))
                    cds = ((160, 170),) if left else ((485, 495),)
                    inputs = self.inputs(folder, model('partial.1', 'A', strand,
                        exons=exons, cds=cds, five=False, three=False), strand=strand)
                    _, outputs, audit = self.protect(folder, inputs)
                    index = 0 if kind == 'FE' else 1
                    result = outputs[index].read_text().strip().split('\t')
                    self.assertEqual(result[1 if left else 2], '185' if left else '470')
                    self.assertEqual(result[8], inputs[index].read_text().strip().split('\t')[8])
                    self.assertEqual(audit[0]['end'], kind)
                    self.assertIn('partial_cds:partial.1', audit[0]['blocked_by'])
                    self.assertEqual(outputs[1-index].read_bytes(), inputs[1-index].read_bytes())

    def test_annotation_only_extension_is_unchanged_even_with_partial_contacts(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            partial = model('partial.1', 'B', exons=((150, 190),), cds=((160, 180),), five=False, three=False)
            inputs = self.inputs(folder, partial)
            inputs[0].write_text('1\t150\t250\t.\t0\t+\t.\t.\t3\tA\n')
            _, outputs, audit = self.protect(folder, inputs)
            self.assertEqual(outputs[0].read_bytes(), inputs[0].read_bytes())
            self.assertEqual(audit, [])

    def test_annotation_baseline_is_kept_when_all_read_candidates_fail(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            partial = model('partial.1', exons=((150, 195), (300, 350)), cds=((160, 195),), five=False, three=False)
            inputs = self.inputs(folder, partial)
            _, outputs, audit = self.protect(folder, inputs)
            fields = outputs[0].read_text().strip().split('\t')
            self.assertEqual(fields[1:5], ['200', '250', '.', '0'])
            self.assertEqual(audit[0]['outcome'], 'annotation_fallback')

    def test_longer_stale_options_never_revive(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            partial = model('partial.1', exons=((150, 190),), cds=((180, 190),), five=False, three=False)
            inputs = self.inputs(folder, partial)
            # 100 survived in metadata after an earlier terminal fix selected 185.
            inputs[0].write_text('1\t185\t250\t100/195/200\t2\t+\t.\t.\t3\tA\n')
            _, outputs, _ = self.protect(folder, inputs)
            fields = outputs[0].read_text().strip().split('\t')
            self.assertEqual(fields[1], '195')
            self.assertEqual(fields[3], '200')

    def test_shared_splice_partial_is_not_an_endpoint_donor(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            inputs = self.inputs(folder, model('partial.1', exons=((100, 250), (400, 450)),
                                                cds=((220, 240),), five=False, three=False))
            inputs[0].write_text('1\t200\t250\t.\t0\t+\t.\t.\t0\tA\n')
            _, outputs, audit = self.protect(folder, inputs)
            self.assertEqual(outputs[0].read_bytes(), inputs[0].read_bytes())
            self.assertEqual(audit, [])

    def test_opposite_strand_or_other_chromosome_cds_does_not_block(self):
        for extra in [dict(strand='-'), dict(chromosome='2')]:
            with self.subTest(extra=extra), tempfile.TemporaryDirectory() as tmp:
                folder = Path(tmp)
                partial = model('partial.1', exons=((100, 190),), cds=((160, 190),), five=False, three=False, **extra)
                inputs = self.inputs(folder, partial)
                _, outputs, audit = self.protect(folder, inputs)
                self.assertEqual(outputs[0].read_bytes(), inputs[0].read_bytes())
                self.assertEqual(audit, [])

    def test_unknown_partial_terminal_contact_blocks_but_original_natural_pair_allows(self):
        for strand in ['+', '-']:
            for natural in ['', 'A\tB\n', 'B\tA\n']:
                with self.subTest(strand=strand, natural=natural), tempfile.TemporaryDirectory() as tmp:
                    folder = Path(tmp)
                    # Left extension is FE on + and LE on -; the ordered pair reverses.
                    partial = model('partial.1', 'B', strand, exons=((160, 180),), cds=(), five=False, three=False)
                    inputs = self.inputs(folder, partial, strand=strand, natural=natural)
                    _, outputs, audit = self.protect(folder, inputs)
                    index = 0 if strand == '+' else 1
                    allowed = natural == ('A\tB\n' if strand == '+' else 'B\tA\n')
                    self.assertEqual(outputs[index].read_text().split('\t')[1], '150' if allowed else '185')
                    self.assertEqual(bool(audit), not allowed)

    def test_multi_exon_terminal_contact_only_uses_opposite_end(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            partial = model('partial.1', 'B', exons=((50, 90), (160, 180)), cds=(), five=False, three=False)
            inputs = self.inputs(folder, partial)
            _, outputs, audit = self.protect(folder, inputs)
            self.assertEqual(outputs[0].read_text().split('\t')[1], '185')
            self.assertIn('partial_terminal_LE', audit[0]['blocked_by'])

    def test_raw_terminal_contact_outside_added_read_region_does_not_block(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            partial = model('partial.1', 'B', exons=((225, 235),), cds=(), five=False, three=False)
            inputs = self.inputs(folder, partial)
            _, outputs, audit = self.protect(folder, inputs)
            self.assertEqual(outputs[0].read_bytes(), inputs[0].read_bytes())
            self.assertEqual(audit, [])

    def test_same_kind_terminal_contact_is_not_a_new_blanket_exon_rule(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            partial = model('partial.1', 'B', exons=((160, 180), (300, 350)), cds=(), five=False, three=False)
            inputs = self.inputs(folder, partial)
            _, outputs, audit = self.protect(folder, inputs)
            self.assertEqual(outputs[0].read_bytes(), inputs[0].read_bytes())
            self.assertEqual(audit, [])

    def test_gene_groups_retain_original_ordered_zip_and_shared_gene_allowances(self):
        for genes, pairs, allowed in [('X,A', 'A\tB\n', False),
                                      ('X,A', 'X\tB\n', True),
                                      ('B,A', '', True)]:
            with self.subTest(genes=genes, pairs=pairs), tempfile.TemporaryDirectory() as tmp:
                folder = Path(tmp)
                partial = model('partial.1', 'B', exons=((160, 180),), cds=(), five=False, three=False)
                inputs = self.inputs(folder, partial, natural=pairs)
                inputs[0].write_text('1\t150\t250\t185/200\t1\t+\t.\t.\t3\t'+genes+'\n')
                _, outputs, audit = self.protect(folder, inputs)
                self.assertEqual(outputs[0].read_text().split('\t')[1], '150' if allowed else '185')
                self.assertEqual(bool(audit), not allowed)

    def test_single_cds_classification_is_frozen_and_unique_partial_flanks_are_protected(self):
        spec = importlib.util.spec_from_file_location('protection', SCRIPT)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            full = model('source.1', cds=((410, 420),))
            partial = model('partial.1', exons=((100, 170), (180, 190), (400, 450)),
                            cds=((410, 415),), five=False, three=False)
            inputs = self.inputs(folder, partial, full=full, natural='A\tB\n')
            protection = module.Protection(inputs[2], inputs[3], {('A', 'B')})
            self.assertEqual(protection.special, {('1', '+', 'A')})
            _, outputs, audit = self.protect(folder, inputs)
            self.assertEqual(outputs[0].read_text().split('\t')[1], '200')
            self.assertIn('partial_single_cds_flank', audit[0]['blocked_by'])

    def test_cds_boundary_touch_matches_existing_guard_semantics(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            partial = model('partial.1', exons=((200, 210),), cds=((200, 210),), five=False, three=False)
            inputs = self.inputs(folder, partial)
            _, outputs, audit = self.protect(folder, inputs)
            self.assertEqual(outputs[0].read_text().split('\t')[1], '200')
            self.assertTrue(audit)

    def test_terminal_boundary_touch_is_not_positive_length_contact(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            partial = model('partial.1', 'B', exons=((200, 210),), cds=(), five=False, three=False)
            inputs = self.inputs(folder, partial)
            _, outputs, audit = self.protect(folder, inputs)
            self.assertEqual(outputs[0].read_bytes(), inputs[0].read_bytes())
            self.assertEqual(audit, [])

    def test_gzip_and_input_errors(self):
        for invalid in [False, True]:
            with self.subTest(invalid=invalid), tempfile.TemporaryDirectory() as tmp:
                folder = Path(tmp)
                inputs = self.inputs(folder)
                zipped = folder / 'partial.gff.gz'
                with gzip.open(zipped, 'wt') as handle:
                    handle.write('invalid\n' if invalid else '')
                inputs[3] = zipped
                result, outputs, audit = self.protect(folder, inputs, ok=not invalid)
                self.assertEqual(result.returncode, 2 if invalid else 0)
                if invalid:
                    self.assertIn('expected 9 GFF columns', result.stderr)
                    self.assertFalse(outputs[0].exists())

    def test_builder_changes_only_the_protected_end_and_its_support(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            partial = model('partial.1', exons=((150, 195),), cds=((160, 195),), five=False, three=False)
            inputs = self.inputs(folder, partial)
            _, outputs, _ = self.protect(folder, inputs)
            blacklist = folder / 'blacklist'; blacklist.touch()
            output, stats = folder / 'new.gff', folder / 'new.tsv'
            # The final writer derives strand from the FE filename.
            fe = folder / '1_positive_FE.bed'; fe.write_bytes(outputs[0].read_bytes())
            result = subprocess.run([sys.executable, str(ROOT/'bin/build_chromosome_gff_and_stats.py'),
                str(fe), str(outputs[1]), str(inputs[2]), str(blacklist), str(output), str(stats),
                '--partial-annotation', str(inputs[3])], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            with stats.open() as handle:
                row = next(csv.DictReader(handle, delimiter='\t'))
            self.assertEqual((row['five_prime_extension'], row['five_prime_support'], row['five_prime_read_count']), ('0', 'none', '0'))
            self.assertEqual((row['three_prime_extension'], row['three_prime_support'], row['three_prime_read_count']), ('50', 'reads', '4'))
            self.assertEqual({a['transcript_id'] for f,a in parse(output) if f[2]!='gene'}, {'source.1'})


if __name__ == '__main__':
    unittest.main()
