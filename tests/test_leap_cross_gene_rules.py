"""LEAP clipping and terminal-safety predicates, independently of old exon audits."""

from collections import Counter
import csv
import importlib.util
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from leap_test_context import write_context, check

ROOT = Path(__file__).resolve().parents[1]


def model(
    tid, gene, exons, cds=(), utr5=(), utr3=(), strand='+', biotype='protein_coding'
):
    features = [('transcript', (exons[0][0], exons[-1][1]))]
    features += [('exon', x) for x in exons] + [('CDS', x) for x in cds]
    features += [('five_prime_UTR', x) for x in utr5] + [
        ('three_prime_UTR', x) for x in utr3
    ]
    return ''.join(
        '1\ttest\t{}\t{}\t{}\t.\t{}\t.\ttranscript_id={};gene_id={};transcript_type={};{}\n'.format(
            kind,
            *interval,
            strand,
            tid,
            gene,
            biotype,
            (
                'exon_number={};'.format(
                    exons.index(interval) + 1
                    if strand == '+'
                    else len(exons) - exons.index(interval)
                )
                if kind == 'exon'
                else ''
            ),
        )
        for kind, interval in features
    )


def source(start=200, strand='+'):
    return model(
        'source',
        'A',
        [(start, 250), (400, 450)],
        [(221, 250), (400, 429)],
        [(start, 220)],
        [(430, 450)],
        strand,
    )


def terminal_target(end=180, strand='+'):
    return model(
        'target',
        'B',
        [(80, 100), (140, end)],
        [(90, 100), (140, 149)],
        [(80, 89)],
        [(150, end)] if end >= 150 else [(end, end)],
        strand,
    )


class LeapRulesTests(unittest.TestCase):
    def evaluate(self, raw, extended=None, natural='', edit=None, stats=None):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            r, e = directory / 'raw.gff3', directory / 'extended.gff3'
            r.write_text(raw)
            e.write_text(extended or source(150))
            context = write_context(directory, r, e)
            for p in (context / 'natural_overlaps').glob('*'):
                p.write_text(natural)
            if edit:
                edit(context)
            original, new = check.read_annotation(r), check.read_annotation(e)
            run = check.LeapContext(original, new, context)
            counts = Counter()
            rows = list(
                check.compare(
                    original, new, counts, extension_stats=stats, leap_context=run
                )
            )
            return rows, counts

    def test_clipping_eligibility_and_independent_exon_audit(self):
        full = model(
            'target',
            'B',
            [(140, 180), (300, 350)],
            [(160, 170), (300, 330)],
            [(140, 159)],
            [(331, 350)],
        )
        rows, counts = self.evaluate(source() + full)
        failure = next(r for r in rows if r['failing_intervals'])
        self.assertEqual(
            (
                failure['check_rule'],
                failure['target_feature'],
                failure['failing_intervals'],
            ),
            ('full_length_cds', 'CDS', '160-170'),
        )
        self.assertEqual(failure['gene_relationship'], 'cross_gene')
        self.assertEqual(counts['new_cross_gene_cds_overlap'], 1)
        self.assertTrue(any(r['status'] == 'allowed_leap_overlap' for r in rows))
        partial = full.replace('three_prime_UTR', 'five_prime_UTR')
        rows, _ = self.evaluate(source() + partial)
        self.assertEqual({r['check_rule'] for r in rows}, {'ineligible_partial_model'})
        self.assertEqual(rows[0]['target_is_partial'], 'true')
        self.assertEqual(rows[0]['target_missing_utr'], 'three_prime_UTR')
        npc = model('target', 'B', [(160, 180)], biotype='lncRNA')
        rows, _ = self.evaluate(source() + npc)
        self.assertEqual(rows[0]['check_rule'], 'noncoding_exon')
        self.assertEqual(rows[0]['failing_intervals'], '160-180')
        se = model(
            'target', 'B', [(140, 190)], [(160, 170)], [(140, 159)], [(171, 190)]
        )
        rows, _ = self.evaluate(source() + se)
        self.assertEqual(
            {r['check_rule'] for r in rows if r['failing_intervals']},
            {'single_exon_cds', 'terminal_single_exon:unknown_pair'},
        )

    def test_terminal_pairs_order_and_shared_groups_on_both_strands(self):
        for strand in ('+', '-'):
            raw = source(strand=strand) + terminal_target(strand=strand)
            ext = source(150, strand) + terminal_target(strand=strand)
            pair = 'A\tB\n' if strand == '+' else 'B\tA\n'
            for natural, expected in [
                ('', False),
                (pair, True),
                ('B\tA\n' if strand == '+' else 'A\tB\n', False),
            ]:
                with self.subTest(strand=strand, natural=natural):
                    rows, counts = self.evaluate(raw, ext, natural)
                    self.assertEqual(counts['unsafe_terminal_overlap'] == 0, expected)
                    if expected:
                        self.assertTrue(
                            any(
                                'recorded_natural_pair' in r['check_rule'] for r in rows
                            )
                        )

            def groups(context):
                for path in (context / 'fixed_extensions').glob('*'):
                    path.write_text(
                        path.read_text()
                        .replace('\tA\n', '\tA,C\n')
                        .replace('\tB\n', '\tB,C\n')
                    )

            rows, counts = self.evaluate(raw, ext, edit=groups)
            self.assertEqual(counts['unsafe_terminal_overlap'], 0)
            self.assertTrue(any('shared_gene_group' in r['check_rule'] for r in rows))

        def zipped(context):
            for path in (context / 'fixed_extensions').glob('*'):
                path.write_text(
                    path.read_text()
                    .replace('\tA\n', '\tA,C\n')
                    .replace('\tB\n', '\tD,B\n')
                )

        _, counts = self.evaluate(
            source() + terminal_target(),
            source(150) + terminal_target(),
            'A\tB\n',
            zipped,
        )
        self.assertEqual(
            counts['unsafe_terminal_overlap'], 1
        )  # zip, not Cartesian product.

    def test_annotation_bypasses_clipping_but_not_terminal_safety(self):
        def annotation(context):
            for path in (context / 'fixed_extensions').glob('*FE*'):
                fields = [line.split('\t') for line in path.read_text().splitlines()]
                for f in fields:
                    if f[-1] == 'A':
                        f[3] = '.'
                path.write_text(''.join('\t'.join(f) + '\n' for f in fields))

        full = model(
            'target',
            'B',
            [(140, 180), (300, 350)],
            [(160, 170), (300, 330)],
            [(140, 159)],
            [(331, 350)],
        )
        rows, counts = self.evaluate(source() + full, edit=annotation)
        self.assertEqual(counts['new_cross_gene_cds_overlap'], 0)
        self.assertTrue(any(r['status'] == 'allowed_current_annotation' for r in rows))
        rows, counts = self.evaluate(
            source() + terminal_target(),
            source(150) + terminal_target(),
            edit=annotation,
        )
        self.assertEqual(counts['unsafe_terminal_overlap'], 1)
        self.assertTrue(
            any(r['check_rule'] == 'terminal_fe_le:unknown_pair' for r in rows)
        )

    def test_bed_touch_and_inclusive_cds_guard_and_grouped_original_boundary(self):
        _, counts = self.evaluate(
            source() + terminal_target(150), source(150) + terminal_target(150)
        )
        self.assertEqual(counts['unsafe_terminal_overlap'], 0)
        npc = model('target', 'B', [(200, 200)], biotype='lncRNA')
        rows, _ = self.evaluate(source() + npc)
        self.assertEqual(
            rows[0]['failing_intervals'], '200-200'
        )  # clipping includes the old boundary.

        def baseline(context):
            path = next((context / 'fixed_extensions').glob('*FE*'))
            path.write_text(path.read_text().replace('\t200\t0\t', '\t170/190\t0\t'))

        rows, counts = self.evaluate(source() + npc, edit=baseline)
        self.assertEqual(counts['new_cross_gene_exon_overlap'], 0)
        self.assertEqual(
            rows, []
        )  # 200 was already annotation, beyond the grouped candidate guard.

    def test_terminal_safety_checks_whole_candidate_exon(self):
        # Even a contact already inside this isoform's exon must pass the separate terminal predicate.
        raw_source = model(
            'source',
            'A',
            [(160, 250), (400, 450)],
            [(221, 250), (400, 429)],
            [(160, 220)],
            [(430, 450)],
        )
        rows, counts = self.evaluate(
            raw_source + terminal_target(), source(150) + terminal_target()
        )
        self.assertEqual(counts['unsafe_terminal_overlap'], 1)
        self.assertEqual(
            next(r for r in rows if r['failing_intervals'])['failing_intervals'],
            '150-180',
        )

    def test_same_terminal_kind_and_internal_exons_are_not_terminal_tests(self):
        target = model(
            'target',
            'B',
            [(140, 180), (300, 350), (500, 550)],
            [(310, 330)],
            [(140, 180)],
            [(520, 550)],
        )
        _, counts = self.evaluate(source() + target, source(150) + target)
        self.assertEqual(counts['unsafe_terminal_overlap'], 0)
        target = model(
            'target',
            'B',
            [(50, 80), (160, 180), (500, 550)],
            [(510, 520)],
            [(50, 80)],
            [(521, 550)],
        )
        rows, counts = self.evaluate(source() + target, source(150) + target)
        self.assertEqual(counts['unsafe_terminal_overlap'], 0)
        self.assertEqual(rows[0]['check_rule'], 'outside_leap_protected_regions')

    def test_special_single_cds_gene_flanks_and_exact_group_trigger(self):
        # Source gene A itself must be a single-CDS gene appearing in a natural pair.
        a = model(
            'source',
            'A',
            [(200, 250), (400, 450)],
            [(410, 420)],
            [(200, 250)],
            [(421, 450)],
        )
        ae = model(
            'source',
            'A',
            [(150, 250), (400, 450)],
            [(410, 420)],
            [(150, 250)],
            [(421, 450)],
        )
        b = model(
            'target',
            'B',
            [(160, 180), (300, 350), (500, 550)],
            [(310, 330)],
            [(160, 180)],
            [(500, 550)],
        )
        rows, counts = self.evaluate(a + b, ae, 'A\tB\n')
        self.assertEqual(
            next(r for r in rows if r['failing_intervals'])['check_rule'],
            'single_cds_gene_flank',
        )
        for natural in ('', 'B\tB\n'):
            _, counts = self.evaluate(a + b, ae, natural)
            self.assertEqual(counts['new_cross_gene_exon_overlap'], 0)

        def grouped(context):
            for path in (context / 'fixed_extensions').glob('*'):
                path.write_text(path.read_text().replace('\tA\n', '\tA,C\n'))

        _, counts = self.evaluate(a + b, ae, 'A\tB\n', grouped)
        self.assertEqual(counts['new_cross_gene_exon_overlap'], 0)
        # Directly inspect both strands' FE/LE flank choices against transcript order.
        for strand in ('+', '-'):
            with tempfile.TemporaryDirectory() as directory:
                r, e = Path(directory) / 'r', Path(directory) / 'e'

                def numbered(text, count):
                    text = text.replace('\t+\t', f'\t{strand}\t')
                    if strand == '-':
                        import re

                        text = re.sub(
                            r'exon_number=(\d+);',
                            lambda m: 'exon_number={};'.format(count + 1 - int(m[1])),
                            text,
                        )
                    return text

                r.write_text(numbered(a, 2) + numbered(b, 3))
                e.write_text(numbered(ae, 2))
                context = write_context(directory, r, e)
                next((context / 'natural_overlaps').glob('*')).write_text('A\tB\n')
                ctx = check.LeapContext(
                    check.read_annotation(r), check.read_annotation(e), context
                )
                self.assertIn(
                    (
                        *((160, 180) if strand == '+' else (500, 550)),
                        'target',
                        'exon',
                        'single_cds_gene_flank',
                    ),
                    list(ctx.flanks[('1', strand, 'FE')].query(1, 1000)),
                )
                self.assertIn(
                    (
                        *((500, 550) if strand == '+' else (160, 180)),
                        'target',
                        'exon',
                        'single_cds_gene_flank',
                    ),
                    list(ctx.flanks[('1', strand, 'LE')].query(1, 1000)),
                )

    def test_only_partial_cds_has_the_extra_within_gene_rule(self):
        target = model(
            'target',
            'A',
            [(140, 180), (300, 350)],
            [(160, 170), (300, 330)],
            [(140, 159)],
            [(331, 350)],
        )
        # A read extension into eligible same-gene CDS is tested by LEAP clipping.
        rows, counts = self.evaluate(source() + target)
        failing = [r for r in rows if r['failing_intervals']]
        self.assertEqual(len(failing), 1)
        self.assertEqual(failing[0]['check_rule'], 'full_length_cds')
        self.assertEqual(failing[0]['intentional_partial_cds_check'], 'false')
        self.assertEqual(counts['new_cds_overlap'], 1)

        # LEAP clips only beyond the grouped annotation boundary, not the shorter isoform boundary.
        def grouped_boundary(context):
            path = next((context / 'fixed_extensions').glob('*FE*'))
            path.write_text(path.read_text().replace('\t200\t0\t', '\t140\t0\t'))

        rows, counts = self.evaluate(source() + target, edit=grouped_boundary)
        self.assertEqual(counts['new_cds_overlap'], 0)
        self.assertTrue(
            any(r['check_rule'] == 'outside_leap_extension_guard' for r in rows)
        )
        # A recorded boundary alone cannot excuse a partial-CDS hit: no eligible
        # raw terminal donor shares this source splice site at the claimed bound.
        partial = target.replace('three_prime_UTR', 'five_prime_UTR')
        rows, counts = self.evaluate(source() + partial, edit=grouped_boundary)
        self.assertEqual(counts['new_cds_overlap'], 1)
        self.assertEqual(rows[0]['check_rule'], 'within_gene_partial_cds')
        self.assertEqual(rows[0]['intentional_partial_cds_check'], 'true')

    def grouped_partial_case(self, side, strand, cds, donor_kind='eligible', selected=None, recorded=None):
        def coding(tid, start, end, gene='A', biotype='protein_coding'):
            low, high = [(start, 219)], [(431, end)]
            return model(tid, gene, [(start, 250), (400, end)], [(220, 240), (410, 430)],
                         low if strand == '+' else high, high if strand == '+' else low,
                         strand, biotype)
        raw_source = coding('source', 200, 450)
        boundary = 170 if side == 'left' else 480
        endpoint = selected if selected is not None else (150 if side == 'left' else 500)
        extended = coding('source', endpoint if side == 'left' else 200,
                          endpoint if side == 'right' else 450)
        donor = coding('donor', 170 if side == 'left' else 200, 480 if side == 'right' else 450)
        if donor_kind == 'partial':
            donor = donor.replace('three_prime_UTR', 'five_prime_UTR')
        elif donor_kind == 'noncoding':
            donor = donor.replace('protein_coding', 'nonsense_mediated_decay')
        elif donor_kind == 'other_gene':
            donor = donor.replace('gene_id=A;', 'gene_id=B;')
        elif donor_kind == 'other_strand':
            donor = donor.replace('\t' + strand + '\t', '\t' + ('-' if strand == '+' else '+') + '\t')
        elif donor_kind == 'other_chromosome':
            donor = donor.replace('1\ttest', '2\ttest')
        elif donor_kind == 'other_site':
            donor = donor.replace('\t250\t', '\t251\t') if side == 'left' else donor.replace('\t400\t', '\t399\t')
        elif donor_kind == 'source_splice_change':
            extended = extended.replace('\t250\t', '\t251\t')
            donor = donor.replace('\t250\t', '\t251\t')
        elif donor_kind == 'internal_exon':
            donor = model('donor', 'A', [(100, 120), (170, 250), (600, 650)],
                          [(180, 200)], [(100, 120)], [(630, 650)], strand)
        elif donor_kind == 'single_exon':
            donor = model('donor', 'A', [(170, 250)], [(200, 220)], [(170, 199)], [(221, 250)], strand)
        target = model('target', 'A', [(cds[0]-1, cds[1]+1)], [cds], strand=strand)
        kind = check.terminal_kind(strand, side)
        def edit(context):
            path = next((context/'fixed_extensions').glob('*_'+kind+'_*'))
            f = path.read_text().rstrip().split('\t')
            f[3] = str(boundary if recorded is None else recorded)
            path.write_text('\t'.join(f)+'\n')
        start, end = (endpoint,450) if side == 'left' else (200,endpoint)
        five,three = (200-start,end-450) if strand == '+' else (end-450,200-start)
        stats = {('1','source'): dict(gene_id='A',transcript_location='1:{}-{}'.format(start,end),
                 five_prime_extension=five,three_prime_extension=three,
                 five_prime_support='reads,current_annotation' if five else 'none',
                 three_prime_support='reads,current_annotation' if three else 'none')}
        rows, counts = self.evaluate(raw_source+donor+target, extended, edit=edit, stats=stats)
        row = next(r for r in rows if r['target_transcript']=='target' and r['check_rule']=='within_gene_partial_cds')
        return row, counts

    def test_grouped_annotation_partial_cds_allowance_both_strands_and_ends(self):
        for strand in ('+','-'):
            for side, cds, expected in [('left',(170,185),'170-185'), ('right',(465,480),'465-480')]:
                with self.subTest(strand=strand, side=side):
                    row,counts = self.grouped_partial_case(side,strand,cds)
                    self.assertEqual(row['status'],'allowed_grouped_annotation')
                    self.assertEqual(row['grouped_annotation_allowed_intervals'],expected)
                    self.assertEqual(row['allowed_intervals'],expected)
                    self.assertEqual(row['failing_intervals'],'')
                    self.assertEqual(row['annotation_allowed_intervals'],'')
                    self.assertEqual(row['grouped_annotation_donors'],'donor')
                    prime = 'five_prime' if check.terminal_kind(strand,side)=='FE' else 'three_prime'
                    self.assertEqual(row['grouped_annotation_boundaries'],prime+(':'+'170' if side=='left' else ':480'))
                    self.assertEqual(counts['partial_cds_failures'],0)

    def test_partial_cds_crossing_grouped_boundary_keeps_read_portion_failing(self):
        for strand in ('+','-'):
            for side, cds, allowed, failing in [('left',(165,185),'170-185','165-169'),
                                               ('right',(465,485),'465-480','481-485')]:
                row,counts = self.grouped_partial_case(side,strand,cds)
                self.assertEqual(row['status'],'new_cds_overlap')
                self.assertEqual(row['grouped_annotation_allowed_intervals'],allowed)
                self.assertEqual(row['failing_intervals'],failing)
                self.assertEqual(row['failing_bases'],5)
                self.assertEqual(row['added_overlap_bases'],21)
                self.assertEqual(counts['partial_cds_failures'],1)
            for side, cds in [('left',(169,169)),('right',(481,481))]:
                row,counts = self.grouped_partial_case(side,strand,cds)
                self.assertEqual(row['failing_bases'],1)
                self.assertEqual(row['grouped_annotation_allowed_intervals'],'')

    def test_only_matching_eligible_raw_terminal_donors_can_allow_partial_contacts(self):
        for donor in ('partial','noncoding','single_exon','internal_exon','other_gene',
                      'other_strand','other_chromosome','other_site','source_splice_change'):
            with self.subTest(donor=donor):
                row,counts = self.grouped_partial_case('left','+',(170,185),donor_kind=donor)
                self.assertEqual(row['status'],'new_cds_overlap')
                self.assertEqual(row['grouped_annotation_allowed_intervals'],'')
                self.assertEqual(counts['partial_cds_failures'],1)

    def test_group_allowance_is_capped_by_selected_output_and_recorded_baseline(self):
        for strand in ('+','-'):
            row,_ = self.grouped_partial_case('left',strand,(180,185),selected=180)
            self.assertEqual(row['grouped_annotation_allowed_intervals'],'180-185')
            row,_ = self.grouped_partial_case('right',strand,(465,470),selected=470)
            self.assertEqual(row['grouped_annotation_allowed_intervals'],'465-470')
        # Annotation fallback removes the stored alternatives; raw donor evidence
        # still identifies the grouped boundary without exempting read sequence.
        row,_ = self.grouped_partial_case('left','+',(170,185),recorded='.')
        self.assertEqual(row['status'],'allowed_grouped_annotation')
        row,counts = self.grouped_partial_case('left','+',(165,169),recorded='.')
        self.assertEqual(row['failing_intervals'],'165-169')
        row,counts = self.grouped_partial_case('left','+',(170,185),recorded=190)
        self.assertEqual(row['failing_intervals'],'170-185')
        self.assertEqual(row['grouped_annotation_allowed_intervals'],'')
        self.assertEqual(counts['partial_cds_failures'],1)

    def test_run_candidate_controls_clipping_not_stats_support(self):
        target = model(
            'target',
            'B',
            [(140, 180), (300, 350)],
            [(160, 170), (300, 330)],
            [(140, 159)],
            [(331, 350)],
        )
        stats = {
            ('1', 'source'): dict(
                gene_id='A',
                transcript_location='1:150-450',
                five_prime_extension=50,
                three_prime_extension=0,
                five_prime_support='current_annotation',
                three_prime_support='none',
            )
        }
        # Stats cannot turn the original run's read-derived candidate into an annotation candidate.
        rows, counts = self.evaluate(source() + target, stats=stats)
        self.assertEqual(counts['new_cross_gene_cds_overlap'], 1)
        self.assertEqual(
            next(r for r in rows if r['failing_intervals'])['check_rule'],
            'full_length_cds',
        )

    def test_leap_noncoding_guard_applies_inside_the_source_gene_too(self):
        target = model('target', 'A', [(160, 180)], biotype='processed_transcript')
        rows, counts = self.evaluate(source() + target)
        self.assertEqual(counts['new_within_gene_exon_overlap'], 1)
        self.assertEqual(rows[0]['check_rule'], 'noncoding_exon')
        self.assertEqual(rows[0]['gene_relationship'], 'within_gene')
        self.assertEqual(rows[0]['intentional_partial_cds_check'], 'false')

    def test_annotation_candidate_bypasses_eligible_same_gene_cds(self):
        target = model(
            'target',
            'A',
            [(140, 180), (300, 350)],
            [(160, 170), (300, 330)],
            [(140, 159)],
            [(331, 350)],
        )

        def annotation(context):
            path = next((context / 'fixed_extensions').glob('*FE*'))
            path.write_text(path.read_text().replace('\t200\t0\t', '\t.\t0\t'))

        rows, counts = self.evaluate(source() + target, edit=annotation)
        self.assertEqual(counts['new_cds_overlap'], 0)
        self.assertTrue(
            any(
                r['status'] == 'allowed_current_annotation'
                and r['check_rule'] == 'full_length_cds'
                for r in rows
            )
        )

    def test_saved_scnna_grouped_annotation_cli_and_missing_raw_donor(self):
        fixture = ROOT/'tests/data/leap_cross_gene/grouped_annotation_partial_cds'
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            raw = root/'raw.gff3';raw.write_text((fixture/'raw.gff3').read_text())
            output = root/'report.csv'
            command = [sys.executable,str(ROOT/'assets/check_cds_overlaps.py'),str(raw),
                       str(fixture/'extended.gff3'),str(output),'--leap-output-dir',str(fixture/'leap_run'),
                       '--extension-stats',str(fixture/'stats.tsv')]
            result = subprocess.run(command,capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stderr)
            with output.open() as handle:
                row = next(r for r in csv.DictReader(handle) if r['target_transcript']=='ENST00000536788.1')
            self.assertEqual(row['status'],'allowed_grouped_annotation')
            self.assertEqual(row['grouped_annotation_allowed_intervals'],'6374804-6374837')
            self.assertEqual(row['grouped_annotation_boundaries'],'five_prime:6375224')
            self.assertEqual(row['grouped_annotation_donors'],'ENST00000360168.8')
            self.assertEqual(row['failing_bases'],'0')
            self.assertIn('allowed_grouped_annotation: 1',result.stdout)
            self.assertIn('partial_cds_failures: 0',result.stdout)
            # The saved BED baseline alone cannot create the new exemption.
            raw.write_text(''.join(line for line in raw.read_text().splitlines(True)
                                   if 'transcript_id=ENST00000360168.8;' not in line))
            result = subprocess.run(command,capture_output=True,text=True)
            self.assertEqual(result.returncode,1,result.stderr)
            self.assertIn('partial_cds_failures: 1',result.stdout)

    def test_saved_cli_examples_and_seven_partial_cds_regression(self):
        data = ROOT / 'tests/data/leap_cross_gene'
        expected = {
            'full_length_cds': 1,
            'noncoding_exon': 1,
            'single_exon': 1,
            'unknown_terminal_pair': 1,
            'natural_terminal_pair': 0,
            'annotation_unsafe_terminal': 1,
            'bed_touch': 0,
            'seven_partial_cds_failures': 1,
            'shared_first_exon_different_genes': 0,
        }
        for name, code in expected.items():
            with self.subTest(name=name), tempfile.TemporaryDirectory() as directory:
                fixture = data / name
                output = Path(directory) / 'report.csv'
                command = [
                    sys.executable,
                    str(ROOT / 'assets/check_cds_overlaps.py'),
                    str(fixture / 'raw.gff3'),
                    str(fixture / 'extended.gff3'),
                    str(output),
                    '--leap-output-dir',
                    str(fixture / 'leap_run'),
                ]
                if (fixture / 'stats.tsv').exists():
                    command += ['--extension-stats', str(fixture / 'stats.tsv')]
                result = subprocess.run(command, capture_output=True, text=True)
                self.assertEqual(result.returncode, code, result.stderr)
                with output.open() as handle:
                    rows = list(csv.DictReader(handle))
                if name == 'seven_partial_cds_failures':
                    failures = [r for r in rows if r['failing_intervals']]
                    self.assertEqual(len(failures), 7)
                    self.assertEqual(
                        {r['target_transcript'] for r in failures},
                        {'ENST00000534890.1', 'ENST00000585355.5'},
                    )
                    self.assertEqual(
                        {r['target_missing_utr'] for r in failures}, {'three_prime_UTR'}
                    )
                    self.assertEqual(
                        {r['gene_relationship'] for r in failures}, {'within_gene'}
                    )
                    self.assertIn('new_cds_overlap: 7', result.stdout)
                    self.assertIn('leap_rule_failures: 0', result.stdout)
                    self.assertIn('partial_cds_failures: 7', result.stdout)
                    self.assertEqual(
                        {r['check_rule'] for r in failures}, {'within_gene_partial_cds'}
                    )
                    self.assertEqual(
                        {r['intentional_partial_cds_check'] for r in failures}, {'true'}
                    )
                if name == 'shared_first_exon_different_genes':
                    self.assertTrue(rows)
                    self.assertEqual(
                        {r['gene_relationship'] for r in rows}, {'cross_gene'}
                    )
                    self.assertTrue(all(not r['failing_intervals'] for r in rows))
                    self.assertTrue(
                        any(
                            r['source_transcript'] == 'ENST00000241356.5'
                            and r['target_transcript'] == 'ENST00001035321.1'
                            for r in rows
                        )
                    )
                if name == 'full_length_cds':
                    self.assertIn('new_cross_gene_cds_overlap: 1', result.stdout)
                if name == 'annotation_unsafe_terminal':
                    self.assertEqual(
                        rows[0]['check_rule'], 'terminal_fe_le:unknown_pair'
                    )

    def test_decisions_match_production_leap_predicates(self):
        import ast

        # Execute the real CDS clipping predicate with a tiny frame adapter; no pandas dependency.
        tree = ast.parse((ROOT / 'bin/stop_cds_overlap.py').read_text())
        node = next(
            n
            for n in tree.body
            if isinstance(n, ast.FunctionDef)
            and n.name == 'overlap_region_with_feature'
        )
        namespace = {}
        exec(
            compile(
                ast.Module(body=[node], type_ignores=[]),
                'LEAP clipping predicate',
                'exec',
            ),
            namespace,
        )

        class Vector(list):
            def __lt__(self, n):
                return Vector(x < n for x in self)

            def __gt__(self, n):
                return Vector(x > n for x in self)

            def __or__(self, other):
                return Vector(a or b for a, b in zip(self, other))

        class Frame:
            def __init__(self, intervals):
                self.intervals = intervals

            def __len__(self):
                return len(self.intervals)

            def __getitem__(self, key):
                if isinstance(key, str):
                    return Vector(x[0 if key == 'start' else 1] for x in self.intervals)
                return Frame([x for x, keep in zip(self.intervals, key) if keep])

        for strand in ('positive', 'negative'):
            for kind in ('FE', 'LE'):
                left = (strand == 'positive') == (kind == 'FE')
                fields = [
                    '1',
                    '150' if left else '400',
                    '250' if left else '500',
                    '200' if left else '450',
                    '0',
                    '+' if strand == 'positive' else '-',
                    '.',
                    '.',
                    '2',
                    'A',
                ]
                guard = (150, 200) if left else (450, 500)
                for target in [
                    (guard[0] - 1, guard[0] - 1),
                    (guard[0], guard[0]),
                    (guard[1], guard[1]),
                    (guard[1] + 1, guard[1] + 1),
                ]:
                    for annotation in (False, True):
                        candidate = fields.copy()
                        if annotation:
                            candidate[3] = '.'
                        result = namespace['overlap_region_with_feature'](
                            Frame([target]), candidate, kind, strand
                        )
                        overlaps = target[0] <= guard[1] and target[1] >= guard[0]
                        clipped = (int(result[1]), int(result[2])) != (
                            int(fields[1]),
                            int(fields[2]),
                        )
                        self.assertEqual(clipped, overlaps and not annotation)
        # Execute LEAP's actual terminal checker, including its ordered groups and SE field positions.
        cases = [
            (('A',), ('B',), set()),
            (('A',), ('B',), {('A', 'B')}),
            (('A',), ('B',), {('B', 'A')}),
            (('A', 'C'), ('B', 'C'), set()),
            (('A', 'C'), ('D', 'B'), {('A', 'B')}),
            (('A', 'C'), ('D', 'B'), {('C', 'B')}),
        ]
        for fe, le, pairs in cases:
            for single in (False, True):
                with tempfile.TemporaryDirectory() as directory:
                    d = Path(directory)
                    natural = d / 'natural'
                    overlap = d / 'overlap'
                    unsafe = d / 'unsafe'
                    natural.write_text(''.join(a + '\t' + b + '\n' for a, b in pairs))
                    fields = [
                        '1',
                        '100',
                        '200',
                        'single_exon' if single else '.',
                        '0',
                        '+',
                        ','.join(fe) if single else '.',
                        '.',
                        '0',
                        ','.join(fe),
                    ]
                    fields += [
                        '1',
                        '150',
                        '250',
                        '.',
                        '0',
                        '+',
                        '.',
                        '.',
                        '0',
                        ','.join(le),
                    ]
                    overlap.write_text('\t'.join(fields) + '\n')
                    for name in ('fe', 'le'):
                        (d / (name + '.bed')).write_text('')
                    result = subprocess.run(
                        [
                            sys.executable,
                            str(ROOT / 'bin/overlap_safety_checker.py'),
                            str(natural),
                            str(overlap),
                            str(d / 'fe'),
                            str(d / 'le'),
                            str(unsafe),
                        ],
                        capture_output=True,
                        text=True,
                    )
                    self.assertEqual(result.returncode, 0, result.stderr)
                    self.assertEqual(
                        unsafe.stat().st_size == 0,
                        bool(check.terminal_safe(fe, le, pairs)),
                    )

    def test_checker_clipping_matches_leap_for_both_gene_relationships(self):
        import ast

        tree = ast.parse((ROOT / 'bin/stop_cds_overlap.py').read_text())
        node = next(
            n
            for n in tree.body
            if isinstance(n, ast.FunctionDef)
            and n.name == 'overlap_region_with_feature'
        )
        ns = {}
        exec(
            compile(
                ast.Module(body=[node], type_ignores=[]),
                'production LEAP clipping',
                'exec',
            ),
            ns,
        )

        class Vector(list):
            def __lt__(self, n):
                return Vector(x < n for x in self)

            def __gt__(self, n):
                return Vector(x > n for x in self)

            def __or__(self, other):
                return Vector(a or b for a, b in zip(self, other))

        class Frame:
            def __init__(self, intervals):
                self.intervals = intervals

            def __len__(self):
                return len(self.intervals)

            def __getitem__(self, key):
                if isinstance(key, str):
                    return Vector(x[0 if key == 'start' else 1] for x in self.intervals)
                return Frame([x for x, keep in zip(self.intervals, key) if keep])

        for strand in ('+', '-'):
            for side in ('left', 'right'):
                kind = check.terminal_kind(strand, side)
                raw_source = source(strand=strand)
                new_source = (
                    source(150, strand)
                    if side == 'left'
                    else model(
                        'source',
                        'A',
                        [(200, 250), (400, 500)],
                        [(221, 250), (400, 429)],
                        [(200, 220)],
                        [(430, 500)],
                        strand,
                    )
                )
                new_start, new_end = (150, 250) if side == 'left' else (400, 500)
                baseline = 200 if side == 'left' else 450
                for same_gene in (False, True):
                    for annotation in (False, True):
                        for feature in ('CDS', 'noncoding_exon', 'single_exon_cds'):
                            for position in (
                                [150, 200, 201] if side == 'left' else [449, 450, 500]
                            ):
                                with self.subTest(
                                    strand=strand,
                                    side=side,
                                    same_gene=same_gene,
                                    annotation=annotation,
                                    feature=feature,
                                    position=position,
                                ):
                                    gene = 'A' if same_gene else 'B'
                                    if feature == 'noncoding_exon':
                                        target = model(
                                            'target',
                                            gene,
                                            [(position, position)],
                                            strand=strand,
                                            biotype='lncRNA',
                                        )
                                    elif feature == 'single_exon_cds':
                                        target = model(
                                            'target',
                                            gene,
                                            [(position - 10, position + 10)],
                                            [(position, position)],
                                            [(position - 10, position - 1)],
                                            [(position + 1, position + 10)],
                                            strand,
                                        )
                                    else:
                                        target = model(
                                            'target',
                                            gene,
                                            [
                                                (position - 10, position + 10),
                                                (700, 750),
                                            ],
                                            [(position, position)],
                                            [(position - 10, position - 1)],
                                            [(740, 750)],
                                            strand,
                                        )

                                    def metadata(context):
                                        path = next(
                                            (context / 'fixed_extensions').glob(
                                                '*_' + kind + '_*'
                                            )
                                        )
                                        fields = path.read_text().strip().split('\t')
                                        fields[3] = '.' if annotation else str(baseline)
                                        path.write_text('\t'.join(fields) + '\n')

                                    rows, _ = self.evaluate(
                                        raw_source + target, new_source, edit=metadata
                                    )
                                    rule = {
                                        'CDS': 'full_length_cds',
                                        'noncoding_exon': 'noncoding_exon',
                                        'single_exon_cds': 'single_exon_cds',
                                    }[feature]
                                    actual = any(
                                        r['failing_intervals']
                                        for r in rows
                                        if r['check_rule'] == rule
                                    )
                                    fields = [
                                        '1',
                                        str(new_start),
                                        str(new_end),
                                        '.' if annotation else str(baseline),
                                        '0',
                                        strand,
                                        '.',
                                        '.',
                                        '2',
                                        'A',
                                    ]
                                    result = ns['overlap_region_with_feature'](
                                        Frame([(position, position)]),
                                        fields.copy(),
                                        kind,
                                        'positive' if strand == '+' else 'negative',
                                    )
                                    clipped = (int(result[1]), int(result[2])) != (
                                        new_start,
                                        new_end,
                                    )
                                    self.assertEqual(actual, clipped)

    def test_report_cannot_overwrite_annotation_stats_or_run_context(self):
        with tempfile.TemporaryDirectory() as directory:
            d = Path(directory)
            r, e = d / 'raw', d / 'extended'
            r.write_text(source())
            e.write_text(source(150))
            context = write_context(d, r, e)
            stats = d / 'stats.tsv'
            stats.write_text(
                'gene_id\ttranscript_id\ttranscript_location\tfive_prime_extension\tthree_prime_extension\tfive_prime_support\tthree_prime_support\nA\tsource\t1:150-450\t50\t0\treads\tnone\n'
            )
            for path in [r, e, stats, *context.rglob('*.bed'), *context.rglob('*.txt')]:
                content = path.read_bytes()
                result = subprocess.run(
                    [
                        sys.executable,
                        str(ROOT / 'assets/check_cds_overlaps.py'),
                        str(r),
                        str(e),
                        str(path),
                        '--leap-output-dir',
                        str(context),
                        '--extension-stats',
                        str(stats),
                    ],
                    capture_output=True,
                    text=True,
                )
                self.assertEqual(result.returncode, 2, result.stderr)
                self.assertIn('must not overwrite', result.stderr)
                self.assertEqual(path.read_bytes(), content)

    def test_duplicate_candidates_use_last_support_on_both_strands_and_ends(self):
        for strand in ('+','-'):
            for side in ('left','right'):
                for last_annotation in (False,True):
                    with self.subTest(strand=strand,side=side,last_annotation=last_annotation):
                        raw = model('source','A',[(200,250),(400,450)],[(221,250),(400,429)],
                                    [(200,220)],[(430,450)],strand)
                        ext = model('source','A',[(150 if side=='left' else 200,250),
                                                (400,500 if side=='right' else 450)],
                                    [(221,250),(400,429)],[(150 if side=='left' else 200,220)],
                                    [(430,500 if side=='right' else 450)],strand)
                        a,b = (160,170) if side=='left' else (480,490)
                        target = model('target','B',[(a-10,b+10),(600,650)],[(a,b)],
                                       [(a-10,a-1)],[(631,650)],strand)
                        kind = check.terminal_kind(strand,side)
                        def duplicate(context):
                            path = next((context/'fixed_extensions').glob('*_'+kind+'_*'))
                            read = path.read_text().rstrip().split('\t')
                            annotation = read.copy();annotation[3]='.'
                            rows = [read,annotation] if last_annotation else [annotation,read]
                            path.write_text(''.join('\t'.join(row)+'\n' for row in rows))
                        rows,counts = self.evaluate(raw+target,ext,edit=duplicate)
                        self.assertEqual(counts['new_cross_gene_cds_overlap'],0 if last_annotation else 1)
                        hit = next(r for r in rows if r['check_rule']=='full_length_cds')
                        self.assertEqual(hit['status'],'allowed_current_annotation' if last_annotation else 'new_cross_gene_cds_overlap')

    def test_duplicate_candidates_keep_last_ordered_gene_group(self):
        for final_groups,safe in [('A,Z',True),('Z,A',False)]:
            def duplicate(context):
                path = next((context/'fixed_extensions').glob('*FE*'))
                rows = []
                for line in path.read_text().splitlines():
                    row = line.split('\t')
                    if row[9]=='A':
                        first = row.copy();first[9]='A,Z'
                        last = row.copy();last[9]=final_groups
                        rows.extend([first,last])
                    else:
                        rows.append(row)
                path.write_text(''.join('\t'.join(row)+'\n' for row in rows))
            rows,counts = self.evaluate(source()+terminal_target(),source(150)+terminal_target(),
                                        natural='A\tB\n',edit=duplicate)
            self.assertEqual(counts['unsafe_terminal_overlap'],0 if safe else 1)

    def test_last_candidate_is_validated_without_falling_back(self):
        for problem in ('gene','geometry'):
            def duplicate(context):
                path = next((context/'fixed_extensions').glob('*FE*'))
                original = path.read_text()
                last = original.rstrip().split('\t')
                if problem=='gene':last[9]='B'
                else:last[1]='170'  # Cannot contain the actual output starting at 150.
                path.write_text(original+'\t'.join(last)+'\n')
            with self.assertRaisesRegex(ValueError, 'splice_site=250.*'+('omits source gene' if problem=='gene' else 'does not contain output exon')):
                self.evaluate(source(),source(150),edit=duplicate)

    def test_processed_end_can_be_contained_and_unchanged_end_can_exceed_candidate(self):
        partial = model('target','A',[(180,195)],[(185,190)])
        def selected(context):
            path = next((context/'fixed_extensions').glob('*FE*'))
            row = path.read_text().rstrip().split('\t');row[1]='150'
            path.write_text('\t'.join(row)+'\n')
        rows,counts = self.evaluate(source()+partial,source(180),edit=selected)
        self.assertEqual(counts['partial_cds_failures'],1)
        self.assertEqual(next(r for r in rows if r['check_rule']=='within_gene_partial_cds')['failing_intervals'],'185-190')
        rows,counts = self.evaluate(source(100),source(100),edit=selected)
        self.assertEqual(rows,[])
        self.assertEqual(counts['leap_rule_failures'],0)

    def test_saved_ddx41_duplicate_annotation_read_candidates(self):
        import shutil
        fixture = ROOT/'tests/data/leap_cross_gene/ddx41_duplicate_candidates'
        with tempfile.TemporaryDirectory() as directory:
            context = Path(directory)/'leap_run'
            shutil.copytree(fixture/'leap_run',context)
            output = Path(directory)/'report.csv'
            command = [sys.executable,str(ROOT/'assets/check_cds_overlaps.py'),str(fixture/'raw.gff3'),
                       str(fixture/'extended.gff3'),str(output),'--extension-stats',str(fixture/'stats.tsv'),
                       '--leap-output-dir',str(context)]
            result = subprocess.run(command,capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stderr)
            self.assertIn('partial_cds_failures: 0',result.stdout)
            raw,new = check.read_annotation(fixture/'raw.gff3'),check.read_annotation(fixture/'extended.gff3')
            source_model = new['ENST00000330503.12']
            ctx = check.LeapContext(raw,new,context)
            self.assertEqual(ctx.candidate(source_model,source_model['exon'][0],'LE'),
                             (('ENSG00000183258.18',),False,177511569))
            bed = context/'fixed_extensions/5_negative_LE_protected.bed'
            bed.write_text(''.join(reversed(bed.read_text().splitlines(True))))
            ctx = check.LeapContext(raw,new,context)
            self.assertEqual(ctx.candidate(source_model,source_model['exon'][0],'LE'),
                             (('ENSG00000183258.18',),True,None))

    def test_required_context_errors(self):
        for error in (
            'natural',
            'fixed_missing',
            'fixed_duplicate',
            'candidate_missing',
            'candidate_wrong_gene',
            'type_missing',
        ):
            with self.subTest(error=error), tempfile.TemporaryDirectory() as directory:
                r, e = Path(directory) / 'r', Path(directory) / 'e'
                r.write_text(source())
                e.write_text(source(150))
                context = write_context(directory, r, e)
                fe = next((context / 'fixed_extensions').glob('*FE*'))
                if error == 'natural':
                    next((context / 'natural_overlaps').glob('*')).unlink()
                if error == 'fixed_missing':
                    fe.unlink()
                if error == 'fixed_duplicate':
                    fe.with_name('duplicate_' + fe.name).write_text(fe.read_text())
                if error == 'candidate_missing':
                    fe.write_text('')
                if error == 'candidate_wrong_gene':
                    fe.write_text(fe.read_text() + fe.read_text().replace('\tA\n', '\tB\n'))
                if error == 'type_missing':
                    r.write_text(
                        r.read_text().replace('transcript_type=protein_coding;', '')
                    )
                output = Path(directory) / 'report.csv'
                result = subprocess.run(
                    [
                        sys.executable,
                        str(ROOT / 'assets/check_cds_overlaps.py'),
                        str(r),
                        str(e),
                        str(output),
                        '--leap-output-dir',
                        str(context),
                    ],
                    capture_output=True,
                    text=True,
                )
                self.assertEqual(result.returncode, 2, result.stderr)
                self.assertIn('ERROR:', result.stderr)
                self.assertFalse(output.exists())


if __name__ == '__main__':
    unittest.main()
