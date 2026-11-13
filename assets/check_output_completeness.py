#!/usr/bin/env python3
"""Independently check published-annotation eligibility against GFF/stats membership."""

import argparse
import csv
import gzip
import json
import math
import re
import sys
import zlib
from collections import Counter
from pathlib import Path

STRANDS = {'positive': '+', 'negative': '-'}
FIELDS = ['row_type', 'chromosome', 'strand', 'transcript_id', 'status', 'failed',
          'expected_retention', 'exclusion_reasons', 'gff_present', 'stats_present',
          'transcript_type', 'exon_count', 'cds_count', 'raw_gff_id', 'extended_gff_id',
          'raw_path', 'gff_path', 'stats_path', 'details']


def open_text(path):
    return (gzip.open if str(path).endswith('.gz') else open)(path, 'rt', encoding='utf-8', newline='')


def attributes(text, context):
    result = {}
    for part in text.split(';'):
        if not part:
            continue
        key, sep, value = part.partition('=')
        if not sep or not key or key in result:
            raise ValueError(context + ': malformed or duplicate attribute')
        result[key] = value
    return result


def read_gff(path, scope, raw=False):
    """Aggregate models without relying on feature order or GFF primary ID suffixes."""
    models, foreign = {}, Counter()
    raw_chromosome_transcript = False
    chromosomes = {c for c, _ in scope}
    with open_text(path) as handle:
        for number, line in enumerate(handle, 1):
            if line.strip() == '##FASTA':
                break
            if not line.strip() or line.startswith('#'):
                continue
            f = line.rstrip('\r\n').split('\t')
            context = '{}:{}'.format(path, number)
            # Raw records outside the requested location are deliberately ignored.
            if raw and len(f) >= 7:
                if f[0] in chromosomes and f[2] == 'transcript':
                    raw_chromosome_transcript = True
                if f[0] not in chromosomes or (f[6] in ('+', '-') and (f[0], f[6]) not in scope):
                    continue
            if len(f) != 9:
                raise ValueError(context + ': expected 9 GFF3 columns')
            try:
                start, end = int(f[3]), int(f[4])
            except ValueError:
                raise ValueError(context + ': coordinates must be integers') from None
            if f[0] in ('', '.') or f[6] not in ('+', '-') or not 1 <= start <= end:
                raise ValueError(context + ': invalid chromosome, strand or interval')
            if f[2] in ('', '.') or not f[1] or f[7] not in ('.', '0', '1', '2'):
                raise ValueError(context + ': invalid feature type, source or phase')
            try:
                if f[5] != '.' and not math.isfinite(float(f[5])):
                    raise ValueError()
            except ValueError:
                raise ValueError(context + ': invalid score') from None
            if (f[0], f[6]) not in scope:
                foreign[(f[0], f[6])] += 1
            a = attributes(f[8], context)
            if f[2] == 'gene':
                continue
            for name in ('transcript_id', 'gene_id'):
                if not a.get(name, '').strip() or a[name] == '.':
                    raise ValueError(context + ': missing ' + name)
            key = (f[0], f[6], a['transcript_id'])
            m = models.setdefault(key, dict(transcript=None, exons=[], features=[], seen=set(),
                                           types=set(), genes=set(), utr=set(), cds=0, path=str(path)))
            fingerprint = tuple(f[:8]) + (tuple(sorted(a.items())),)
            if fingerprint in m['seen']:
                raise ValueError(context + ': duplicate feature record for ' + key[2])
            m['seen'].add(fingerprint)
            m['genes'].add(a['gene_id'])
            if 'transcript_type' in a:
                if a['transcript_type'] in ('', '.'):
                    raise ValueError(context + ': empty transcript_type')
                m['types'].add(a['transcript_type'])
            if f[2] == 'transcript':
                if m['transcript'] is not None:
                    raise ValueError(context + ': duplicate transcript record for ' + key[2])
                m['transcript'] = (start, end, a, context)
            else:
                m['features'].append((f[2], start, end, context))
                if f[2] == 'exon':
                    m['exons'].append((start, end))
                if f[2] in ('five_prime_UTR', 'three_prime_UTR'):
                    m['utr'].add(f[2])
                if f[2] == 'CDS':
                    m['cds'] += 1
    for key, m in models.items():
        if m['transcript'] is None:
            raise ValueError(str(path) + ': orphan child records for ' + str(key))
        start, end, a, context = m['transcript']
        if not m['exons']:
            raise ValueError(context + ': transcript has no exons: ' + key[2])
        m['exons'].sort()
        if any(x[1] >= y[0] for x, y in zip(m['exons'], m['exons'][1:])):
            raise ValueError(context + ': duplicate or overlapping exons: ' + key[2])
        if len(m['genes']) != 1 or len(m['types']) > 1:
            raise ValueError(context + ': inconsistent gene_id or transcript_type: ' + key[2])
        if raw and (not a.get('transcript_type') or len(m['types']) != 1):
            raise ValueError(context + ': missing/unclassifiable transcript_type: ' + key[2])
        for kind, left, right, location in m['features']:
            if not start <= left <= right <= end:
                raise ValueError(location + ': feature outside transcript: ' + key[2])
            if kind in ('CDS', 'five_prime_UTR', 'three_prime_UTR') and not any(
                    x <= left <= right <= y for x, y in m['exons']):
                raise ValueError(location + ': ' + kind + ' outside exons: ' + key[2])
        m['type'] = a.get('transcript_type', '')
        m['id'] = a.get('ID', '')
        m['exon_count'] = len(m['exons'])
        # Release large parsing-only collections before reading the next file.
        for name in ('seen', 'features', 'genes', 'types', 'transcript', 'exons'):
            del m[name]
    if raw and not raw_chromosome_transcript:
        raise ValueError(str(path) + ': no transcript records on selected raw chromosome')
    return models, foreign


def read_stats(path, strand):
    models = {}
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        header = reader.fieldnames or []
        if len(header) != len(set(header)) or not {'gene_id', 'transcript_id', 'transcript_location'} <= set(header):
            raise ValueError(str(path) + ': missing or duplicate stats columns')
        for row in reader:
            context = '{}:{}'.format(path, reader.line_num)
            if None in row or any(row.get(k) in (None, '') for k in header):
                raise ValueError(context + ': incomplete or extra stats fields')
            if any(row[k] == '.' or not row[k].strip() for k in ('gene_id', 'transcript_id')):
                raise ValueError(context + ': missing gene_id or transcript_id')
            match = re.fullmatch(r'([^:]+):([0-9]+)-([0-9]+)', row['transcript_location'])
            if not match or match[1] == '.' or not 1 <= int(match[2]) <= int(match[3]):
                raise ValueError(context + ': invalid transcript_location')
            key = (match[1], strand, row['transcript_id'])
            if key in models:
                raise ValueError(context + ': duplicate stats transcript: ' + str(key))
            models[key] = dict(path=str(path))
    return models


def exclusion_reasons(model, stage):
    reasons = []
    if model['type'] != 'protein_coding':
        reasons.append('non_protein_coding')
    else:
        for kind in ('five_prime_UTR', 'three_prime_UTR'):
            if kind not in model['utr']:
                reasons.append('missing_' + kind)
    if model['exon_count'] == 1:
        reasons.append('single_exon')
    if stage == 'processed' and not reasons and not model['cds']:
        reasons.append('processed_no_CDS')
    return reasons


def discover(root, stems):
    matches = [root / name for name in stems if (root / name).exists()]
    if len(matches) > 1:
        raise ValueError('ambiguous input files: ' + ', '.join(map(str, matches)))
    if matches and not matches[0].is_file():
        raise ValueError(str(matches[0]) + ': expected a regular file')
    return matches[0] if matches else None


def gff_names(stem):
    return [stem + suffix for suffix in ('.gff3', '.gff', '.gff3.gz', '.gff.gz')]


def protect_report(output, paths):
    for path in paths:
        if output.resolve() == path.resolve() or (output.exists() and path.exists() and output.samefile(path)):
            raise ValueError('report would overwrite input: ' + str(path))


def compare(args):
    raw_root, gff_root, stats_root = map(Path, (args.current_annotation, args.extended_annotation, args.extension_stats))
    directory = raw_root.is_dir()
    if directory:
        if any(p.exists() and not p.is_dir() for p in (gff_root, stats_root)):
            raise ValueError('directory mode requires raw, GFF and stats directories')
    elif not raw_root.is_file() or any(p.exists() and not p.is_file() for p in (gff_root, stats_root)):
        raise ValueError('file mode requires raw/GFF/stats file paths')
    if not directory and (len(args.chromosomes) != 1 or args.strand == 'both'):
        raise ValueError('single-file mode requires one chromosome and --strand positive or negative')
    labels = list(STRANDS) if args.strand == 'both' else [args.strand]
    scope = {(chrom, STRANDS[label]) for chrom in args.chromosomes for label in labels}
    raw_slots, slots, paths = {}, {}, []
    for chrom in args.chromosomes:
        raw_slots[chrom] = discover(raw_root, gff_names('current_annotation_' + chrom)) if directory else raw_root
        for label in labels:
            stem = chrom + '_' + label
            gff = discover(gff_root, gff_names(stem + '_extended')) if directory else (gff_root if gff_root.exists() else None)
            stats = discover(stats_root, [stem + '_stats.tsv', stem + '_stats.tsv.gz']) if directory else (stats_root if stats_root.exists() else None)
            slots[(chrom, STRANDS[label])] = (gff, stats)
    # Protect all supplied inputs, including unselected files and missing slots.
    if directory:
        for root in (raw_root, gff_root, stats_root):
            if root.exists():
                paths.extend(p for p in root.iterdir() if p.is_file())
        for chrom in args.chromosomes:
            paths.extend(raw_root / name for name in gff_names('current_annotation_' + chrom))
            for label in labels:
                stem = chrom + '_' + label
                paths.extend(gff_root / name for name in gff_names(stem + '_extended'))
                paths.extend(stats_root / (stem + suffix) for suffix in ('_stats.tsv', '_stats.tsv.gz'))
    else:
        paths.extend((raw_root, gff_root, stats_root))
    paths.extend(p for p in raw_slots.values() if p)
    paths.extend(p for pair in slots.values() for p in pair if p)
    output = Path(args.output)
    protect_report(output, paths)
    rows, raw, extended, stats = [], {}, {}, {}

    def finding(chrom, strand, status, details, **sources):
        rows.append(dict(row_type='file', chromosome=chrom, strand=strand, status=status,
                         failed=True, details=details, **sources))

    for chrom, path in raw_slots.items():
        if path is None:
            finding(chrom, '', 'missing_raw_file', 'Cannot derive expected membership without raw annotation',
                    raw_path=str(raw_root / ('current_annotation_' + chrom + '.gff3')))
        else:
            records, _ = read_gff(path, {(chrom, STRANDS[x]) for x in labels}, raw=True)
            raw.update(records)
    for (chrom, strand), (gff, tsv) in slots.items():
        for kind, path in (('gff', gff), ('stats', tsv)):
            if path is None:
                label = next(k for k, v in STRANDS.items() if v == strand)
                name = chrom + '_' + label + ('_extended.gff' if kind == 'gff' else '_stats.tsv')
                root = gff_root if kind == 'gff' else stats_root
                finding(chrom, strand, 'missing_' + kind + '_file', 'Required output file is absent',
                        **{kind + '_path': str(root / name if directory else root)})
                continue
            if kind == 'gff':
                records, foreign = read_gff(path, {(chrom, strand)})
                for (c, s), count in sorted(foreign.items()):
                    finding(c, s, 'gff_outside_file_scope', '{} records outside {} {}'.format(count, chrom, strand),
                            gff_path=str(path))
                target = extended
            else:
                records = read_stats(path, strand)
                for c in sorted({k[0] for k in records if k[0] != chrom}):
                    finding(c, strand, 'stats_outside_file_scope', 'Stats records outside chromosome ' + chrom,
                            stats_path=str(path))
                target = stats
            overlap = target.keys() & records.keys()
            if overlap:
                raise ValueError(str(path) + ': duplicate transcript across output files: ' + str(sorted(overlap)[0]))
            target.update(records)
    for key in sorted(raw.keys() | extended.keys() | stats.keys()):
        chrom, strand, tid = key
        r, g, s = raw.get(key), extended.get(key), stats.get(key)
        reasons = exclusion_reasons(r, args.stage) if r else []
        known = key[:2] in scope and raw_slots.get(chrom) is not None
        expected = not reasons if r else False
        details = []
        if not known:
            status = 'outside_selected_scope' if key[:2] not in scope else 'raw_unavailable'
            failed = True
        elif expected:
            if not g:
                details.append('missing_from_gff')
            if not s:
                details.append('missing_from_stats')
            status = 'missing_expected' if details else 'retained'
            failed = bool(details)
        else:
            if g:
                details.append('unexpected_in_gff')
            if s:
                details.append('unexpected_in_stats')
            status = 'unexpected_retained' if details else 'expected_exclusion'
            failed = bool(details)
            if not r:
                reasons = ['absent_from_raw_exact_id']
        slot_gff, slot_stats = slots.get(key[:2], (None, None))
        rows.append(dict(row_type='transcript', chromosome=chrom, strand=strand, transcript_id=tid,
                         status=status, failed=failed, expected_retention=expected if known else '',
                         exclusion_reasons=';'.join(reasons), gff_present=bool(g), stats_present=bool(s),
                         transcript_type=r['type'] if r else '', exon_count=r['exon_count'] if r else '',
                         cds_count=r['cds'] if r else '', raw_gff_id=r['id'] if r else '',
                         extended_gff_id=g['id'] if g else '', raw_path=str(raw_slots.get(chrom) or ''),
                         gff_path=g['path'] if g else str(slot_gff or ''),
                         stats_path=s['path'] if s else str(slot_stats or ''), details=';'.join(details)))
    with output.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    counts = Counter(row['status'] for row in rows)
    failures = sum(row['failed'] for row in rows)
    print(('FAIL' if failures else 'PASS') + ': {} failing transcript(s), {} file finding(s)'.format(
        sum(row['failed'] for row in rows if row['row_type'] == 'transcript'),
        sum(row['row_type'] == 'file' for row in rows)))
    print('expected_retained: ' + str(sum(row.get('expected_retention') is True for row in rows)))
    print('gff_transcripts: {}\nstats_transcripts: {}'.format(len(extended), len(stats)))
    for status, count in sorted(counts.items()):
        print('{}: {}'.format(status, count))
    print('exclusion_reasons: ' + json.dumps(dict(sorted(Counter(
        reason for row in rows for reason in row.get('exclusion_reasons', '').split(';') if reason).items()))))
    return int(bool(failures))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('current_annotation')
    parser.add_argument('extended_annotation')
    parser.add_argument('output')
    parser.add_argument('--extension-stats', required=True)
    parser.add_argument('--chromosomes', nargs='+', required=True, help='Explicit chromosome names, space-separated')
    parser.add_argument('--strand', choices=('positive', 'negative', 'both'), default='both')
    parser.add_argument('--stage', choices=('leap', 'processed'), default='leap')
    args = parser.parse_args()
    try:
        if len(set(args.chromosomes)) != len(args.chromosomes) or any(
                not re.fullmatch(r'[A-Za-z0-9_.-]+', c) or c in ('.', '..') for c in args.chromosomes):
            raise ValueError('invalid or duplicate chromosome selection')
        return compare(args)
    except (OSError, ValueError, csv.Error, EOFError, zlib.error) as error:
        print('ERROR: ' + str(error), file=sys.stderr)
        return 2


if __name__ == '__main__':
    sys.exit(main())
