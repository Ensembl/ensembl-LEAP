#!/usr/bin/env python3
"""Independently audit LEAP's partial-protection candidate transition and evidence."""

import argparse
import csv
import gzip
import json
import math
import re
import shutil
import sys
import tempfile
import zlib
from bisect import bisect_right
from collections import Counter, defaultdict
from itertools import zip_longest
from pathlib import Path

AUDIT_FIELDS = ['chromosome', 'strand', 'end', 'splice_site', 'gene_ids',
                'original_selected_endpoint', 'protected_endpoint', 'annotation_endpoint',
                'outcome', 'blocked_by']
FLAGS = ['missing_candidate', 'unexpected_candidate', 'endpoint_mismatch', 'lengthened_endpoint',
         'invented_endpoint', 'protected_fields_changed', 'alternatives_mismatch',
         'unchanged_record_changed', 'candidate_order_mismatch', 'audit_mismatch']
FIELDS = ['row_type', 'chromosome', 'strand', 'end', 'candidate_position', 'splice_site',
          'gene_ids', 'status', 'failed', 'original_endpoint', 'expected_endpoint', 'actual_endpoint',
          'annotation_endpoint', 'expected_outcome', 'blocker_evidence'] + FLAGS + [
          'original_candidate', 'expected_candidate', 'actual_candidate', 'raw_path',
          'before_path', 'after_path', 'audit_path', 'natural_overlaps_path', 'details']
STRANDS = {'positive': '+', 'negative': '-'}


def open_text(path):
    return (gzip.open if str(path).endswith('.gz') else open)(path, 'rt', encoding='utf-8', newline='')


def compact(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'))


def attributes(value, context):
    result = {}
    for item in value.split(';'):
        if not item:
            continue
        key, sep, text = item.partition('=')
        if not sep or not key or key in result:
            raise ValueError(context + ': malformed or duplicate attribute')
        result[key] = text
    return result


def raw_context(path, chromosome, strands):
    """Classify independently; retain exon encounter order for LEAP's FE/LE identity."""
    models, saw_transcript = {}, False
    with open_text(path) as handle:
        for number, line in enumerate(handle, 1):
            if line.strip() == '##FASTA':
                break
            if not line.strip() or line.startswith('#'):
                continue
            f = line.rstrip('\r\n').split('\t')
            context = '{}:{}'.format(path, number)
            if len(f) >= 7:
                if f[0] == chromosome and f[2] == 'transcript':
                    saw_transcript = True
                if f[0] != chromosome or (f[6] in ('+', '-') and f[6] not in strands):
                    continue
            if len(f) != 9:
                raise ValueError(context + ': expected 9 GFF3 columns')
            try:
                start, end = int(f[3]), int(f[4])
                if not 1 <= start <= end or f[6] not in strands:
                    raise ValueError()
            except ValueError:
                raise ValueError(context + ': invalid interval or strand') from None
            if f[2] in ('', '.') or not f[1] or f[7] not in ('.', '0', '1', '2'):
                raise ValueError(context + ': invalid feature type, source or phase')
            try:
                if f[5] != '.' and not math.isfinite(float(f[5])):
                    raise ValueError()
            except ValueError:
                raise ValueError(context + ': invalid score') from None
            a = attributes(f[8], context)
            if f[2] == 'gene':
                continue
            if any(not a.get(k, '').strip() or a[k] == '.' for k in ('transcript_id', 'gene_id')):
                raise ValueError(context + ': missing transcript_id or gene_id')
            key = (f[0], f[6], a['transcript_id'])
            if key not in models:
                models[key] = dict(exons=[], cds=[], utr=set(), utr_intervals=[], genes=set(), types=set(),
                                   transcript=None, bounds=[start, end], context=context)
            m = models[key]
            m['bounds'][0] = min(m['bounds'][0], start)
            m['bounds'][1] = max(m['bounds'][1], end)
            m['genes'].add(a['gene_id'])
            if 'transcript_type' in a:
                if a['transcript_type'] in ('', '.'):
                    raise ValueError(context + ': missing transcript_type')
                m['types'].add(a['transcript_type'])
            if f[2] == 'transcript':
                if m['transcript'] is not None:
                    raise ValueError(context + ': duplicate transcript ' + key[2])
                if not a.get('transcript_type'):
                    raise ValueError(context + ': missing transcript_type')
                m['transcript'] = (start, end)
            elif f[2] == 'exon':
                m['exons'].append((start, end, a.get('exon_number')))
            elif f[2] == 'CDS':
                m['cds'].append((start, end, a.get('exon_number')))
            elif f[2] in ('five_prime_UTR', 'three_prime_UTR'):
                m['utr'].add(f[2])
                m['utr_intervals'].append((start, end))
    if not saw_transcript:
        raise ValueError(str(path) + ': no raw transcripts on selected chromosome')
    full, partial = {}, {}
    for key, m in models.items():
        context = m['context'] + ': ' + key[2]
        if m['transcript'] is None:
            raise ValueError(context + ': orphan child records')
        if not m['exons']:
            raise ValueError(context + ': transcript has no exons')
        ordered = sorted(m['exons'], key=lambda e: e[:2])
        if any(a[1] >= b[0] for a, b in zip(ordered, ordered[1:])):
            raise ValueError(context + ': duplicate or overlapping exons')
        if len(m['types']) != 1 or len(m['genes']) != 1:
            raise ValueError(context + ': inconsistent transcript_type or gene_id')
        if m['bounds'] != list(m['transcript']):
            raise ValueError(context + ': child outside transcript bounds')
        if any(not any(x <= a <= b <= y for x, y, _ in ordered) for a, b in [(a, b) for a, b, _ in m['cds']] + m['utr_intervals']):
            raise ValueError(context + ': CDS/UTR outside exons')
        m['gene'] = next(iter(m['genes']))
        if next(iter(m['types'])) != 'protein_coding':
            continue
        if m['utr'] == {'five_prime_UTR', 'three_prime_UTR'}:
            if len(m['exons']) > 1:
                full[key] = m
        else:
            partial[key] = m
    return full, partial


class IntervalLookup:
    """Independent sweep index; queries use inclusive genomic endpoints."""

    def __init__(self, records):
        self.records = sorted(records)
        self.starts = [r[0] for r in self.records]
        self.high_water = []
        for record in self.records:
            self.high_water.append(max(record[1], self.high_water[-1] if self.high_water else 0))

    def contacts(self, left, right):
        i = bisect_right(self.starts, right)
        while i:
            i -= 1
            if self.high_water[i] < left:
                break
            if self.records[i][1] >= left:
                yield self.records[i]


class Rules:
    def __init__(self, full, partial, pairs):
        self.pairs = pairs
        natural_genes = {g for a, b in pairs if a != b for g in (a, b)}
        gene_cds = defaultdict(set)
        for (chrom, strand, _), m in full.items():
            gene_cds[(chrom, strand, m['gene'])].update((a, b) for a, b, _ in m['cds'])
        self.single_cds_genes = {key for key, values in gene_cds.items()
                                 if len(values) == 1 and key[2] in natural_genes}
        cds, terminal, flank = defaultdict(list), defaultdict(list), defaultdict(list)
        for (chrom, strand, tid), m in partial.items():
            gene = m['gene']
            cds[(chrom, strand)].extend((a, b, tid, gene) for a, b, _ in m['cds'])
            exons = m['exons']
            if len(exons) == 1:
                terminal[(chrom, strand, 'SE')].append((*exons[0][:2], tid, gene))
            else:
                terminal[(chrom, strand, 'FE')].append((*exons[0][:2], tid, gene))
                terminal[(chrom, strand, 'LE')].append((*exons[-1][:2], tid, gene))
            if (chrom, strand, gene) not in self.single_cds_genes or not m['cds']:
                continue
            numbered = exons + m['cds']
            try:
                if any(n is None or not math.isfinite(float(n)) for _, _, n in numbered):
                    raise ValueError()
                # LEAP's special-gene predicate uses string-sorted CDS exon numbers.
                first_coding = float(min(n for _, _, n in m['cds']))
                coding_numbers = {n for _, _, n in m['cds']}
                coding_exons = [e for e in exons if e[2] in coding_numbers]
                for kind in ('FE', 'LE'):
                    selected = [e for e in exons if (float(e[2]) < first_coding if kind == 'FE'
                                                     else float(e[2]) > first_coding)] or coding_exons
                    flank[(chrom, strand, kind)].extend((a, b, tid, gene) for a, b, _ in selected)
            except (TypeError, ValueError):
                raise ValueError(m['context'] + ': exon_number required for special-gene flanks: ' + tid) from None
        self.cds = {key: IntervalLookup(rows) for key, rows in cds.items()}
        self.terminal = {key: IntervalLookup(rows) for key, rows in terminal.items()}
        self.flank = {key: IntervalLookup(rows) for key, rows in flank.items()}

    def blockers(self, row, kind, endpoint, baseline):
        location = (row[0], row[5])
        genes = tuple(row[9].split(','))
        left, right = sorted((endpoint, baseline))
        found = set()
        checks = [('partial_cds', self.cds.get(location))]
        if len(genes) == 1 and location + genes in self.single_cds_genes:
            checks.append(('partial_single_cds_flank', self.flank.get(location + (kind,))))
        for rule, lookup in checks:
            if lookup:
                for a, b, tid, _ in lookup.contacts(left, right):
                    found.add('{}:{}:{}-{}'.format(rule, tid, a, b))
        for target in (('LE' if kind == 'FE' else 'FE'), 'SE'):
            lookup = self.terminal.get(location + (target,))
            if not lookup:
                continue
            for a, b, tid, gene in lookup.contacts(left, right):
                if min(right, b) <= max(left, a):
                    continue
                first, last = (genes, (gene,)) if kind == 'FE' else ((gene,), genes)
                if set(first).intersection(last) or any(pair in self.pairs for pair in zip(first, last)):
                    continue
                found.add('partial_terminal_{}:{}:{}-{}'.format(target, tid, a, b))
        return found

    def expected(self, row, kind):
        """Derive a transition from source alternatives; never trust recorded audit fields."""
        left = (kind == 'FE') == (row[5] == '+')
        outer, inner = (1, 2) if left else (2, 1)
        selected = int(row[outer])
        if row[3] in ('.', 'single_exon'):
            return row[:], None
        alternatives = list(map(int, row[3].split('/')))
        baseline = alternatives[-1]
        if (left and selected >= baseline) or (not left and selected <= baseline):
            return row[:], None
        options = [selected]
        for value in alternatives[:-1]:
            if (selected < value < baseline if left else baseline < value < selected) and value not in options:
                options.append(value)
        evidence = set()
        chosen = baseline
        for endpoint in options:
            blocked = self.blockers(row, kind, endpoint, baseline)
            if not blocked:
                chosen = endpoint
                break
            evidence.update(blocked)
        if chosen == selected:
            return row[:], None
        result = row[:]
        result[outer] = str(chosen)
        if chosen == baseline:
            result[3:5] = ['.', '0']
        else:
            remaining = [p for p in options if (p > chosen if left else p < chosen)] + [baseline]
            result[3:5] = ['/'.join(map(str, remaining)), str(len(remaining) - 1)]
        audit = dict(chromosome=row[0], strand=row[5], end=kind, splice_site=row[inner], gene_ids=row[9],
                     original_selected_endpoint=str(selected), protected_endpoint=str(chosen),
                     annotation_endpoint=str(baseline),
                     outcome='annotation_fallback' if chosen == baseline else 'shorter_read_candidate',
                     blocked_by=';'.join(sorted(evidence)))
        return result, audit


def read_bed(path, chromosome, strand):
    records = []
    with open_text(path) as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            f = line.rstrip('\r\n').split('\t')
            context = '{}:{}'.format(path, number)
            try:
                if len(f) != 10 or f[0] != chromosome or f[5] != strand or not 0 <= int(f[1]) <= int(f[2]):
                    raise ValueError()
                if not re.fullmatch(r'[0-9]+', f[4]):
                    raise ValueError()
                if f[3] not in ('.', 'single_exon') and any(int(v) < 1 for v in f[3].split('/')):
                    raise ValueError()
                genes = f[6] if f[3] == 'single_exon' else f[9]
                if any(g in ('', '.') for g in genes.split(',')):
                    raise ValueError()
                if f[8] != '.' or f[3] != 'single_exon':
                    if not re.fullmatch(r'[0-9]+', f[8]):
                        raise ValueError()
            except ValueError:
                raise ValueError(context + ': invalid fixed-extension BED record or file scope') from None
            records.append((f, line))
    return records


def read_pairs(path):
    pairs = set()
    with open_text(path) as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            f = line.rstrip('\r\n').split('\t')
            if len(f) != 2 or any(v in ('', '.') for v in f):
                raise ValueError('{}:{}: expected two ordered gene IDs'.format(path, number))
            pairs.add(tuple(f))
    return pairs


def audit_tuple(row):
    return tuple(row[name] for name in AUDIT_FIELDS)


def read_audit(path):
    result = Counter()
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        header = reader.fieldnames or []
        if len(header) != len(set(header)) or not set(AUDIT_FIELDS) <= set(header):
            raise ValueError(str(path) + ': missing or duplicate audit columns')
        for row in reader:
            context = '{}:{}'.format(path, reader.line_num)
            if None in row or any(row.get(name) in (None, '') for name in header):
                raise ValueError(context + ': incomplete or extra audit fields')
            if row['strand'] not in ('+', '-') or row['end'] not in ('FE', 'LE') or any(
                    not re.fullmatch(r'[0-9]+', row[name]) for name in
                    ('splice_site', 'original_selected_endpoint', 'protected_endpoint', 'annotation_endpoint')):
                raise ValueError(context + ': invalid audit location or integer')
            # Reasons are a set; require no duplicate claims but accept their ordering.
            reasons = row['blocked_by'].split(';')
            if len(set(reasons)) != len(reasons) or any(not r for r in reasons):
                raise ValueError(context + ': duplicate or empty blocker reason')
            row['blocked_by'] = ';'.join(sorted(reasons))
            result[audit_tuple(row)] += 1
    return result


def one_file(root, predicate, label):
    if not root.is_dir():
        raise ValueError(str(root) + ': missing required context directory')
    paths = sorted(p for p in root.iterdir() if p.is_file() and predicate(p.name))
    if len(paths) != 1:
        raise ValueError('{}: expected exactly one {}; found {}'.format(root, label, len(paths)))
    return paths[0]


def discover(root, chromosomes, labels):
    result = {}
    for chrom in chromosomes:
        raw = one_file(root / 'raw_annotation', lambda n: n in [
            'current_annotation_' + chrom + suffix for suffix in ('.gff', '.gff3', '.gff.gz', '.gff3.gz')],
            'raw annotation for ' + chrom)
        for label in labels:
            slot = dict(raw=raw)
            stem = chrom + '_' + label
            for kind in ('FE', 'LE'):
                pattern = re.compile(r'(?:^|_)' + re.escape(stem + '_' + kind) + '_')
                for phase, folder in (('before', 'fixed_extensions_before_partial_protection'), ('after', 'fixed_extensions')):
                    slot[phase + '_' + kind] = one_file(root / folder,
                        lambda n, p=pattern: n.endswith(('.bed', '.bed.gz')) and bool(p.search(n)), phase + ' ' + stem + '_' + kind)
            slot['audit'] = one_file(root / 'partial_protection', lambda n: n in (
                stem + '_partial_protection.tsv', stem + '_partial_protection.tsv.gz'), stem + ' audit')
            slot['natural'] = one_file(root / 'natural_overlaps', lambda n: n in (
                stem + '_NO.txt', stem + '_NO.txt.gz'), stem + ' natural pairs')
            result[(chrom, STRANDS[label])] = slot
    return result


def source_fields(slot, kind):
    return dict(raw_path=str(slot['raw']), before_path=str(slot['before_' + kind]),
                after_path=str(slot['after_' + kind]), audit_path=str(slot['audit']),
                natural_overlaps_path=str(slot['natural']))


def inspect_partition(chrom, strand, slot, full, partial):
    rules = Rules(full, partial, read_pairs(slot['natural']))
    expected_audit = Counter()
    for kind in ('FE', 'LE'):
        before = read_bed(slot['before_' + kind], chrom, strand)
        after = read_bed(slot['after_' + kind], chrom, strand)
        expected_rows = []
        for position, pair in enumerate(zip_longest(before, after), 1):
            original, actual = pair
            row = dict(row_type='candidate', chromosome=chrom, strand=strand, end=kind,
                       candidate_position=position, **source_fields(slot, kind))
            for flag in FLAGS:
                row[flag] = False
            if original is None:
                row.update(unexpected_candidate=True, status='unexpected_candidate',
                           actual_candidate='\t'.join(actual[0]), failed=True,
                           details='No original candidate at this record position')
                yield row
                continue
            f, original_line = original
            expected, audit = rules.expected(f, kind)
            expected_rows.append(tuple(expected))
            if audit:
                expected_audit[audit_tuple(audit)] += 1
            left = (kind == 'FE') == (strand == '+')
            outer, inner = (1, 2) if left else (2, 1)
            selected = int(f[outer])
            baseline = int(f[3].split('/')[-1]) if f[3] not in ('.', 'single_exon') else selected
            row.update(splice_site=f[inner], gene_ids=f[6] if f[3] == 'single_exon' else f[9],
                       original_endpoint=selected, expected_endpoint=int(expected[outer]),
                       annotation_endpoint=baseline, expected_outcome=audit['outcome'] if audit else 'unchanged',
                       blocker_evidence=audit['blocked_by'] if audit else '', original_candidate='\t'.join(f),
                       expected_candidate='\t'.join(expected))
            if actual is None:
                row.update(missing_candidate=True, failed=True, status='missing_candidate', details='Original candidate is missing')
                yield row
                continue
            a, actual_line = actual
            endpoint = int(a[outer])
            row.update(actual_endpoint=endpoint, actual_candidate='\t'.join(a))
            row['endpoint_mismatch'] = a[outer] != expected[outer]
            row['lengthened_endpoint'] = endpoint < selected if left else endpoint > selected
            legal = {selected}
            if f[3] not in ('.', 'single_exon'):
                legal.update(v for v in map(int, f[3].split('/')) if
                             (selected <= v <= baseline if left else baseline <= v <= selected))
            row['invented_endpoint'] = endpoint not in legal
            row['protected_fields_changed'] = any(a[i] != f[i] for i in range(10) if i not in (outer, 3, 4))
            row['alternatives_mismatch'] = a[3:5] != expected[3:5]
            row['unchanged_record_changed'] = not audit and actual_line != original_line
            row['failed'] = any(row[flag] for flag in FLAGS)
            row['status'] = 'invalid_transition' if row['failed'] else (audit['outcome'] if audit else 'unchanged')
            row['details'] = compact(dict(flags=[flag for flag in FLAGS if row[flag]])) if row['failed'] else ''
            yield row
        if (len(before) == len(after) and Counter(expected_rows) == Counter(tuple(r[0]) for r in after)
                and expected_rows != [tuple(r[0]) for r in after]):
            yield dict(row_type='file', chromosome=chrom, strand=strand, end=kind,
                       status='candidate_order_mismatch', failed=True, candidate_order_mismatch=True,
                       details='Candidate order changed; final writer lookup is order-sensitive', **source_fields(slot, kind))
    observed = read_audit(slot['audit'])
    for status, difference in (('missing_audit', expected_audit - observed), ('unexpected_audit', observed - expected_audit)):
        for key, count in sorted(difference.items()):
            audit = dict(zip(AUDIT_FIELDS, key))
            yield dict(row_type='audit', chromosome=chrom, strand=strand, end=audit['end'],
                       splice_site=audit['splice_site'], gene_ids=audit['gene_ids'], status=status,
                       failed=True, audit_mismatch=True, blocker_evidence=audit['blocked_by'],
                       details=compact(dict(occurrences=count, audit=audit)), **source_fields(slot, audit['end']))


def run(args):
    if len(set(args.chromosomes)) != len(args.chromosomes) or any(
            not re.fullmatch(r'[A-Za-z0-9_.-]+', c) or c in ('.', '..') for c in args.chromosomes):
        raise ValueError('invalid or duplicate chromosome selection')
    root, output = Path(args.leap_output_dir), Path(args.output)
    labels = list(STRANDS) if args.strand == 'both' else [args.strand]
    slots = discover(root, args.chromosomes, labels)
    inputs = {path for slot in slots.values() for path in slot.values()}
    for folder in ('raw_annotation', 'fixed_extensions_before_partial_protection', 'fixed_extensions', 'partial_protection', 'natural_overlaps'):
        inputs.update(p for p in (root / folder).iterdir() if p.is_file())
    for source in inputs:
        if output.resolve() == source.resolve() or (output.exists() and output.samefile(source)):
            raise ValueError('report would overwrite input: ' + str(source))
    counts, flags = Counter(), Counter()
    with tempfile.TemporaryFile(mode='w+', newline='', encoding='utf-8') as spool:
        writer = csv.DictWriter(spool, fieldnames=FIELDS)
        writer.writeheader()
        for chrom in args.chromosomes:
            full, partial = raw_context(slots[(chrom, STRANDS[labels[0]])]['raw'], chrom,
                                        {STRANDS[label] for label in labels})
            for label in labels:
                strand = STRANDS[label]
                for row in inspect_partition(chrom, strand, slots[(chrom, strand)], full, partial):
                    writer.writerow(row)
                    counts[row['status']] += 1
                    counts['failing_' + row['row_type']] += int(row['failed'])
                    counts['candidate_records'] += int(row['row_type'] == 'candidate')
                    for flag in FLAGS:
                        flags[flag] += int(bool(row.get(flag)))
        spool.seek(0)
        with output.open('w', newline='', encoding='utf-8') as destination:
            shutil.copyfileobj(spool, destination)
    failed = sum(counts['failing_' + kind] for kind in ('candidate', 'audit', 'file'))
    print(('FAIL' if failed else 'PASS') + ': {} failing candidate(s), {} audit finding(s), {} file finding(s)'.format(
        counts['failing_candidate'], counts['failing_audit'], counts['failing_file']))
    print('partitions: {}'.format(len(slots)))
    for name, count in sorted(counts.items()):
        print('{}: {}'.format(name, count))
    for name in FLAGS:
        print('flag_{}: {}'.format(name, flags[name]))
    return int(bool(failed))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('leap_output_dir')
    parser.add_argument('output')
    parser.add_argument('--chromosomes', nargs='+', required=True)
    parser.add_argument('--strand', choices=('positive', 'negative', 'both'), default='both')
    args = parser.parse_args()
    try:
        return run(args)
    except (OSError, ValueError, csv.Error, EOFError, zlib.error) as error:
        print('ERROR: ' + str(error), file=sys.stderr)
        return 2


if __name__ == '__main__':
    sys.exit(main())
