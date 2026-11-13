"""Small original-run contexts for the standalone datacheck CLI tests."""

from collections import defaultdict
import csv
import importlib.util
from pathlib import Path

spec = importlib.util.spec_from_file_location(
    'context_check',
    Path(__file__).resolve().parents[1] / 'assets/check_cds_overlaps.py',
)
check = importlib.util.module_from_spec(spec)
spec.loader.exec_module(check)


def write_context(root, raw_path, extended_path, stats=None):
    root = Path(root) / 'leap_run'
    (root / 'natural_overlaps').mkdir(parents=True)
    (root / 'fixed_extensions').mkdir()
    try:
        raw = check.read_annotation(raw_path)
        extended = check.read_annotation(extended_path)
    except ValueError:
        return root  # Invalid annotation tests must fail in the annotation parser.
    entries = {
        row['transcript_id']: row
        for row in csv.DictReader((stats or '').splitlines(), delimiter='\t')
    }
    rows = defaultdict(list)
    for tid, model in extended.items():
        gene, chromosome, strand = model['location']
        if not model['exon']:
            continue
        for side, exon in [('left', model['exon'][0]), ('right', model['exon'][-1])]:
            kind = check.terminal_kind(strand, side)
            support_end = 'five_prime' if kind == 'FE' else 'three_prime'
            support = entries.get(tid, {}).get(support_end + '_support')
            old = raw.get(tid, model)['exon']
            original = old[0][0] if side == 'left' else old[-1][1]
            fields = [
                chromosome,
                str(exon[0]),
                str(exon[1]),
                '.' if support == 'current_annotation' else str(original),
                '0',
                strand,
                '.',
                '.',
                '0',
                gene,
            ]
            rows[(chromosome, strand, kind)].append('\t'.join(fields) + '\n')
    for chromosome, strand in {m['location'][1:] for m in extended.values()}:
        label = 'positive' if strand == '+' else 'negative'
        (root / 'natural_overlaps' / f'{chromosome}_{label}_NO.txt').write_text('')
        for kind in ('FE', 'LE'):
            (
                root
                / 'fixed_extensions'
                / f'test_{chromosome}_{label}_{kind}_fixed.bed'
            ).write_text(''.join(rows[(chromosome, strand, kind)]))
    return root
