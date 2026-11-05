#!/usr/bin/env python3
"""Use partial models to restrict already selected LEAP read extensions.

Partials never supply endpoints, gene groups, natural pairs or source models.
This runs after the original overlap fixer and cannot lengthen its selections.
"""
import argparse
from bisect import bisect_right
from collections import Counter, defaultdict
import csv
import gzip
from pathlib import Path
import sys


REPORT_FIELDS = [
    'chromosome', 'strand', 'end', 'splice_site', 'gene_ids',
    'original_selected_endpoint', 'protected_endpoint', 'annotation_endpoint',
    'outcome', 'blocked_by',
]


def open_text(path):
    return (gzip.open(path, 'rt', encoding='utf-8', newline='')
            if str(path).endswith('.gz') else open(path, encoding='utf-8', newline=''))


def models(path):
    """Preserve exon input order, as production terminal extraction does."""
    result = {}
    with open_text(path) as handle:
        for number, line in enumerate(handle, 1):
            if line.strip() == '##FASTA':
                break
            if not line.strip() or line.startswith('#'):
                continue
            context = f'{path}:{number}'
            fields = line.rstrip('\r\n').split('\t')
            if len(fields) != 9:
                raise ValueError(context + ': expected 9 GFF columns')
            kind = fields[2]
            if kind not in ('transcript', 'exon', 'CDS'):
                continue
            attrs = {}
            for item in fields[8].split(';'):
                if not item:
                    continue
                name, separator, value = item.partition('=')
                if not separator or not name or name in attrs:
                    raise ValueError(context + ': malformed or duplicate attribute')
                attrs[name] = value
            tid, gene = attrs.get('transcript_id'), attrs.get('gene_id')
            if not tid or tid == '.' or not gene or fields[6] not in ('+', '-'):
                raise ValueError(context + ': transcript_id, gene_id and explicit strand required')
            try:
                start, end = int(fields[3]), int(fields[4])
                if not fields[0] or fields[0] == '.' or not 1 <= start <= end:
                    raise ValueError()
            except ValueError:
                raise ValueError(context + ': invalid coordinates') from None
            key = (fields[0], fields[6], tid)
            model = result.setdefault(key, dict(gene=gene, exons=[], cds=[], transcript=False))
            if model['gene'] != gene:
                raise ValueError(context + ': inconsistent gene for ' + tid)
            if kind == 'transcript':
                if model['transcript']:
                    raise ValueError(context + ': duplicate transcript ' + tid)
                model['transcript'] = True
            elif kind == 'exon':
                model['exons'].append((start, end, attrs.get('exon_number')))
            else:
                model['cds'].append((start, end, attrs.get('exon_number')))
    for key, model in result.items():
        if not model['transcript']:
            raise ValueError(f'{path}: orphan features for {key[2]}')
        ordered = sorted(model['exons'], key=lambda exon: exon[:2])
        if any(a[1] >= b[0] for a, b in zip(ordered, ordered[1:])):
            raise ValueError(f'{path}: overlapping exons for {key[2]}')
        if any(not any(x <= a <= b <= y for x, y, _ in ordered)
               for a, b, _ in model['cds']):
            raise ValueError(f'{path}: CDS without containing exon for {key[2]}')
    return result


class Intervals:
    def __init__(self, entries):
        self.entries = sorted(entries)
        self.starts = [entry[0] for entry in self.entries]
        self.maximum = []
        maximum = 0
        for entry in self.entries:
            maximum = max(maximum, entry[1])
            self.maximum.append(maximum)

    def query(self, start, end):
        index = bisect_right(self.starts, end) - 1
        while index >= 0 and self.maximum[index] >= start:
            entry = self.entries[index]
            if entry[1] >= start:
                yield entry
            index -= 1


def natural_pairs(path):
    pairs = set()
    with open_text(path) as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            fields = line.rstrip('\r\n').split('\t')
            if len(fields) != 2 or not all(fields):
                raise ValueError(f'{path}:{number}: expected two ordered gene IDs')
            pairs.add(tuple(fields))
    return pairs


class Protection:
    def __init__(self, annotation, partial_annotation, pairs):
        original, partials = models(annotation), models(partial_annotation)
        self.pairs = pairs
        natural_genes = {gene for a, b in pairs if a != b for gene in (a, b)}
        cds_by_gene = defaultdict(set)
        for (chrom, strand, _), model in original.items():
            if len(model['exons']) > 1:
                cds_by_gene[(chrom, strand, model['gene'])].update(
                    (a, b) for a, b, _ in model['cds'])
        # Freeze eligibility to the original full-length reference. Additional
        # partial CDS must not remove an existing single-CDS flank restriction.
        self.special = {key for key, cds in cds_by_gene.items()
                        if len(cds) == 1 and key[2] in natural_genes}
        coding, terminals, flanks = defaultdict(list), defaultdict(list), defaultdict(list)
        for (chrom, strand, tid), model in partials.items():
            gene = model['gene']
            coding[(chrom, strand)].extend((a, b, tid, gene) for a, b, _ in model['cds'])
            exons = model['exons']
            if len(exons) == 1:
                terminals[(chrom, strand, 'SE')].append((*exons[0][:2], tid, gene))
            elif exons:
                terminals[(chrom, strand, 'FE')].append((*exons[0][:2], tid, gene))
                terminals[(chrom, strand, 'LE')].append((*exons[-1][:2], tid, gene))
            if (chrom, strand, gene) not in self.special or not model['cds']:
                continue
            if any(n is None for _, _, n in exons + model['cds']):
                raise ValueError(f'{partial_annotation}: exon_number required for special-gene flanks: {tid}')
            try:
                # Match the existing pandas grouping's string order and its
                # first coding-exon choice, including multiple-CDS partials.
                coding_numbers = sorted({n for _, _, n in model['cds']})
                first_coding = float(coding_numbers[0])
                coding_exons = [exon for exon in exons if exon[2] in coding_numbers]
                for kind in ('FE', 'LE'):
                    selected = [exon for exon in exons
                                if (float(exon[2]) < first_coding if kind == 'FE'
                                    else float(exon[2]) > first_coding)] or coding_exons
                    flanks[(chrom, strand, kind)].extend((a, b, tid, gene) for a, b, _ in selected)
            except ValueError:
                raise ValueError(f'{partial_annotation}: invalid exon_number for {tid}') from None
        self.coding = {key: Intervals(value) for key, value in coding.items()}
        self.terminals = {key: Intervals(value) for key, value in terminals.items()}
        self.flanks = {key: Intervals(value) for key, value in flanks.items()}

    def blocked(self, row, kind, endpoint, baseline):
        chrom, strand = row[0], row[5]
        genes = tuple(row[9].split(','))
        start, end = sorted((endpoint, baseline))
        reasons = set()
        checks = [('partial_cds', self.coding.get((chrom, strand)))]
        if len(genes) == 1 and (chrom, strand, genes[0]) in self.special:
            checks.append(('partial_single_cds_flank', self.flanks.get((chrom, strand, kind))))
        for rule, index in checks:
            if index:
                for a, b, tid, _ in index.query(start, end):
                    # As in stop_cds_overlap.py, the guard includes the
                    # original grouped annotation boundary.
                    reasons.add(f'{rule}:{tid}:{a}-{b}')
        opposite = 'LE' if kind == 'FE' else 'FE'
        for target_kind in (opposite, 'SE'):
            index = self.terminals.get((chrom, strand, target_kind))
            if not index:
                continue
            for a, b, tid, gene in index.query(start, end):
                if max(a, start) >= min(b, end):
                    continue  # Existing terminal checks use positive-length BED contact.
                first, last = (genes, (gene,)) if kind == 'FE' else ((gene,), genes)
                safe = bool(set(first) & set(last)) or any(pair in self.pairs for pair in zip(first, last))
                if not safe:
                    reasons.add(f'partial_terminal_{target_kind}:{tid}:{a}-{b}')
        return reasons

    def protect(self, row, kind):
        if row[3] in ('.', 'single_exon'):
            return row, None
        left = (kind == 'FE') == (row[5] == '+')
        outer, inner = (1, 2) if left else (2, 1)
        selected = int(row[outer])
        alternatives = [int(value) for value in row[3].split('/')]
        baseline = alternatives[-1]
        outward = (lambda value: value < baseline) if left else (lambda value: value > baseline)
        if not outward(selected):
            return row, None
        # The old terminal fixer can leave longer options in field 4. Never
        # revive them, or let another end grow when this end is restricted.
        inward = (lambda value: selected < value < baseline) if left else (lambda value: baseline < value < selected)
        options = list(dict.fromkeys([selected] + [p for p in alternatives[:-1] if inward(p)]))
        reasons = set()
        chosen = baseline
        for endpoint in options:
            blocked = self.blocked(row, kind, endpoint, baseline)
            if not blocked:
                chosen = endpoint
                break
            reasons.update(blocked)
        if chosen == selected:
            return row, None  # Preserve every byte of unchanged BED records.
        result = row.copy()
        result[outer] = str(chosen)
        if chosen == baseline:
            result[3], result[4] = '.', '0'
        else:
            shorter = [p for p in options if (p > chosen if left else p < chosen)] + [baseline]
            result[3] = '/'.join(map(str, shorter))
            result[4] = str(len(shorter) - 1)
        audit = dict(chromosome=row[0], strand=row[5], end=kind,
                     splice_site=row[inner], gene_ids=row[9],
                     original_selected_endpoint=selected, protected_endpoint=chosen,
                     annotation_endpoint=baseline,
                     outcome='annotation_fallback' if chosen == baseline else 'shorter_read_candidate',
                     blocked_by=';'.join(sorted(reasons)))
        return result, audit


def bed_records(path):
    rows = []
    with open_text(path) as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            fields = line.rstrip('\r\n').split('\t')
            try:
                if len(fields) != 10 or fields[5] not in ('+', '-'):
                    raise ValueError()
                if (fields[0] in ('', '.') or not 0 <= int(fields[1]) <= int(fields[2])
                        or any(gene in ('', '.') for gene in fields[9].split(','))):
                    raise ValueError()
                if fields[3] not in ('.', 'single_exon'):
                    if any(int(value) < 1 for value in fields[3].split('/')):
                        raise ValueError()
            except ValueError:
                raise ValueError(f'{path}:{number}: invalid fixed-extension BED record') from None
            rows.append((fields, line))
    return rows


def run(args):
    inputs = [args.first_exons, args.last_exons, args.annotation,
              args.partial_annotation, args.natural_overlaps]
    destinations = [Path(str(args.output_prefix) + suffix)
                    for suffix in ('_FE_protected.bed', '_LE_protected.bed', '_partial_protection.tsv')]
    for path in destinations:
        if path.exists() or any(path.resolve() == source.resolve() for source in inputs):
            raise ValueError('Output already exists or would overwrite an input: ' + str(path))
    protection = Protection(args.annotation, args.partial_annotation, natural_pairs(args.natural_overlaps))
    outputs, audits = [], []
    for kind, path in [('FE', args.first_exons), ('LE', args.last_exons)]:
        lines = []
        for row, original in bed_records(path):
            updated, audit = protection.protect(row, kind)
            lines.append('\t'.join(updated) + '\n' if audit else original)
            if audit:
                audits.append(audit)
        outputs.append(lines)
    # Validate both inputs before opening any output.
    for path, lines in zip(destinations[:2], outputs):
        path.write_text(''.join(lines))
    with destinations[2].open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=REPORT_FIELDS, delimiter='\t')
        writer.writeheader()
        writer.writerows(audits)
    counts = Counter(row['outcome'] for row in audits)
    print(f'Partial protection: {len(audits)} changed candidate records; {dict(counts)}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('first_exons', 'last_exons', 'annotation', 'partial_annotation', 'natural_overlaps', 'output_prefix'):
        parser.add_argument(name, type=Path)
    args = parser.parse_args()
    try:
        run(args)
        return 0
    except (OSError, EOFError, UnicodeError, ValueError, csv.Error) as error:
        print('ERROR: ' + str(error), file=sys.stderr)
        return 2


if __name__ == '__main__':
    sys.exit(main())
