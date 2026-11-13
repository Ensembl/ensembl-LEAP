#!/usr/bin/env python3
"""Check final LEAP GFF, extension stats, tags and optional candidate provenance."""

import argparse
import csv
import gzip
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path


SIDES = ('five_prime', 'three_prime')
STATS_FIELDS = ['gene_id', 'transcript_id', 'transcript_location'] + [
    side + '_' + field for side in SIDES for field in ('extension', 'support', 'read_count')]
FLAGS = ('missing_from_raw', 'missing_stats', 'extra_stats', 'location_mismatch',
         'gene_mismatch', 'shortened', 'exon_span_mismatch', 'extension_length_mismatch',
         'support_mismatch', 'read_count_mismatch', 'leap_tag_mismatch', 'candidate_mismatch')
STATUSES = ('consistent', 'inconsistent', 'missing_from_extended')
FIELDS = ['transcript_id', 'chromosome', 'status', 'failed', 'strand', 'gene_id',
          'raw_start', 'raw_end', 'extended_start', 'extended_end', 'stats_location',
          'observed_leap_tags', 'expected_leap_tags'] + [
    side + '_' + field for side in SIDES
    for field in ('expected_extension', 'reported_extension', 'support', 'read_count')
] + list(FLAGS) + ['candidate_provenance_checked', 'details']


def open_text(path):
    return (gzip.open if str(path).endswith('.gz') else open)(path, 'rt', encoding='utf-8', newline='')


def compact(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'))


def nonnegative(value, context):
    if not re.fullmatch(r'[0-9]+', value):
        raise ValueError(context + ': expected a nonnegative integer')
    return int(value)


def read_annotation(path, include_exons=False):
    """Raw coordinates need transcript rows; final annotation also needs exons."""
    transcripts, exons = {}, defaultdict(list)
    with open_text(path) as handle:
        for number, line in enumerate(handle, 1):
            if line.strip() == '##FASTA':
                break
            if not line.strip() or line.startswith('#'):
                continue
            context = '{}:{}'.format(path, number)
            f = line.rstrip('\r\n').split('\t')
            if len(f) != 9:
                raise ValueError(context + ': expected 9 GFF3 columns')
            if f[2] != 'transcript' and not (include_exons and f[2] == 'exon'):
                continue
            attrs = {}
            for item in f[8].split(';'):
                if not item:
                    continue
                key, separator, value = item.partition('=')
                if not separator or not key or key in attrs:
                    raise ValueError(context + ': malformed or duplicate attribute')
                attrs[key] = value
            if any(not attrs.get(k, '').strip() or attrs[k] == '.' for k in ('transcript_id', 'gene_id')):
                raise ValueError(context + ': missing transcript_id or gene_id')
            try:
                start, end = int(f[3]), int(f[4])
            except ValueError:
                raise ValueError(context + ': coordinates must be integers') from None
            if not 1 <= start <= end or f[0] in ('', '.') or f[6] not in ('+', '-'):
                raise ValueError(context + ': invalid coordinates, chromosome or strand')
            key = (f[0], attrs['transcript_id'])
            record = dict(chromosome=f[0], strand=f[6], gene_id=attrs['gene_id'],
                          start=start, end=end, leap_tag=attrs.get('leap_tag'), context=context)
            if f[2] == 'transcript':
                if key in transcripts:
                    raise ValueError(context + ': duplicate transcript record')
                transcripts[key] = record
            else:
                exons[key].append(record)
    if not transcripts:
        raise ValueError(str(path) + ': no transcript records found')
    if include_exons:
        for key in exons.keys() - transcripts.keys():
            raise ValueError(exons[key][0]['context'] + ': orphan exon record')
        for key, transcript in transcripts.items():
            records = sorted(exons[key], key=lambda r: (r['start'], r['end']))
            if not records:
                raise ValueError(transcript['context'] + ': transcript has no exons: ' + key[1])
            for exon in records:
                if (exon['strand'], exon['gene_id']) != (transcript['strand'], transcript['gene_id']):
                    raise ValueError(exon['context'] + ': inconsistent exon location/gene')
            if any(a['end'] >= b['start'] for a, b in zip(records, records[1:])):
                raise ValueError(transcript['context'] + ': overlapping or duplicate exons')
            transcript['exons'] = records
    return transcripts


def stats_paths(path):
    path = Path(path)
    paths = sorted(p for p in path.iterdir() if p.is_file() and p.name.endswith(('.tsv', '.tsv.gz'))) if path.is_dir() else [path]
    if not paths:
        raise ValueError(str(path) + ': no stats TSV files found')
    return paths


def read_stats(path):
    entries = {}
    for source in stats_paths(path):
        with open_text(source) as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            header = reader.fieldnames or []
            if len(header) != len(set(header)) or not set(STATS_FIELDS) <= set(header):
                raise ValueError(str(source) + ': missing or duplicate stats columns')
            for row in reader:
                context = '{}:{}'.format(source, reader.line_num)
                if None in row or any(row.get(k) in (None, '') for k in header):
                    raise ValueError(context + ': incomplete or extra stats fields')
                for field in ('gene_id', 'transcript_id'):
                    if not row[field].strip() or row[field] == '.':
                        raise ValueError(context + ': missing ' + field)
                try:
                    chromosome, interval = row['transcript_location'].rsplit(':', 1)
                    start, end = map(int, interval.split('-'))
                    if chromosome in ('', '.') or not 1 <= start <= end:
                        raise ValueError()
                except ValueError:
                    raise ValueError(context + ': invalid transcript_location') from None
                key = (chromosome, row['transcript_id'])
                if key in entries:
                    raise ValueError(context + ': duplicate stats transcript key ' + str(key))
                for side in SIDES:
                    row[side + '_extension'] = nonnegative(row[side + '_extension'], context)
                    support = row[side + '_support'].split(',')
                    reads = row[side + '_read_count'].split(',')
                    if any(v not in ('none', 'reads', 'current_annotation', 'blacklisted') for v in support):
                        raise ValueError(context + ': unknown support category')
                    if len(support) != len(reads):
                        raise ValueError(context + ': support/read-count lists have different lengths')
                    for value in reads:
                        if value != 'blacklisted':
                            nonnegative(value, context)
                    row[side + '_read_count'] = ','.join(
                        value if value == 'blacklisted' else str(int(value)) for value in reads)
                entries[key] = row
    if not entries:
        raise ValueError(str(path) + ': no stats rows found')
    return entries


class CandidateContext:
    """Replay the final writer's ordered, last-row-per-splice-site lookup."""

    def __init__(self, root, extended):
        self.records, self.paths = {}, []
        for chromosome, strand in sorted({(r['chromosome'], r['strand']) for r in extended.values()}):
            label = 'positive' if strand == '+' else 'negative'
            for kind in ('FE', 'LE'):
                pattern = re.compile(r'(?:^|_)' + re.escape(chromosome + '_' + label + '_' + kind) + '_')
                files = [p for p in (Path(root) / 'fixed_extensions').iterdir()
                         if p.is_file() and p.name.endswith(('.bed', '.bed.gz')) and pattern.search(p.name)]
                if len(files) != 1:
                    raise ValueError('Expected exactly one fixed {} BED for {} {}'.format(kind, chromosome, strand))
                self.paths.extend(files)
                with open_text(files[0]) as handle:
                    for number, line in enumerate(handle, 1):
                        if not line.strip():
                            continue
                        f = line.rstrip('\r\n').split('\t')
                        context = '{}:{}'.format(files[0], number)
                        try:
                            if len(f) != 10 or f[0] != chromosome or f[5] != strand:
                                raise ValueError()
                            start, end = int(f[1]), int(f[2])
                            if not 0 <= start <= end:
                                raise ValueError()
                            # Single-exon context rows have no splice-site read count.
                            count = None if f[3] == 'single_exon' and f[8] == '.' else nonnegative(f[8], context)
                            if f[3] not in ('.', 'single_exon') and any(int(v) < 1 for v in f[3].split('/')):
                                raise ValueError()
                        except ValueError:
                            raise ValueError(context + ': invalid fixed-extension BED record') from None
                        left = (kind == 'FE') == (strand == '+')
                        inner, outer = (end, start) if left else (start, end)
                        genes = f[6] if f[3] == 'single_exon' else f[9]
                        if any(g in ('', '.') for g in genes.split(',')):
                            raise ValueError(context + ': missing candidate gene IDs')
                        self.records[(chromosome, strand, kind, inner)] = dict(
                            endpoint=outer, support='current_annotation' if f[3] == '.' else 'reads',
                            read_count=count, genes=genes.split(','), context=context)

    def candidate(self, transcript, side):
        kind = 'FE' if side == 'five_prime' else 'LE'
        left = (side == 'five_prime') == (transcript['strand'] == '+')
        exon = transcript['exons'][0 if left else -1]
        inner = exon['end'] if left else exon['start']
        return self.records.get((transcript['chromosome'], transcript['strand'], kind, inner))


def compare(raw, extended, stats, candidates=None):
    for key in sorted(raw.keys() | extended.keys() | stats.keys()):
        old, new, entry = raw.get(key), extended.get(key), stats.get(key)
        row = dict(transcript_id=key[1], chromosome=key[0], candidate_provenance_checked=False)
        details = []
        for flag in FLAGS:
            row[flag] = False

        def issue(flag, **detail):
            row[flag] = True
            details.append(dict(issue=flag, **detail))

        if old:
            row.update(raw_start=old['start'], raw_end=old['end'])
        if new:
            row.update(strand=new['strand'], gene_id=new['gene_id'], extended_start=new['start'], extended_end=new['end'])
        if entry:
            row['stats_location'] = entry['transcript_location']
            for side in SIDES:
                for field in ('extension', 'support', 'read_count'):
                    name = 'reported_extension' if field == 'extension' else field
                    row[side + '_' + name] = entry[side + '_' + field]
        if not new:
            if entry:
                issue('extra_stats', reason='stats row has no final transcript')
        else:
            if old is None:
                issue('missing_from_raw')
            elif (old['strand'], old['chromosome']) != (new['strand'], new['chromosome']):
                issue('location_mismatch', current_strand=old['strand'], extended_strand=new['strand'])
            if entry is None:
                issue('missing_stats')
            elif entry['gene_id'] != new['gene_id']:
                issue('gene_mismatch', reported=entry['gene_id'], expected=new['gene_id'])
            if old and old['gene_id'] != new['gene_id']:
                issue('gene_mismatch', current=old['gene_id'], extended=new['gene_id'])
            location = '{}:{}-{}'.format(new['chromosome'], new['start'], new['end'])
            if entry and entry['transcript_location'] != location:
                issue('location_mismatch', reported=entry['transcript_location'], expected=location)
            exon_span = (new['exons'][0]['start'], new['exons'][-1]['end'])
            if (new['start'], new['end']) != exon_span:
                issue('exon_span_mismatch', transcript=[new['start'], new['end']], exons=exon_span)
            expected_tags = set()
            if old and old['strand'] == new['strand']:
                left, right = old['start'] - new['start'], new['end'] - old['end']
                if min(left, right) < 0:
                    issue('shortened', left_delta=left, right_delta=right)
                lengths = (left, right) if new['strand'] == '+' else (right, left)
                for side, length, tag in zip(SIDES, lengths, ('LEAP_5', 'LEAP_3')):
                    row[side + '_expected_extension'] = length
                    if length > 0:
                        expected_tags.add(tag)
                    if not entry:
                        continue
                    if entry[side + '_extension'] != length:
                        issue('extension_length_mismatch', end=side, expected=length, reported=entry[side + '_extension'])
                    support = entry[side + '_support'].split(',')
                    reads = entry[side + '_read_count'].split(',')
                    blacklisted = support == ['blacklisted']
                    if blacklisted:
                        if side != 'five_prime' or length != 0:
                            issue('support_mismatch', end=side, support=support, expected_extension=length)
                        if reads != ['blacklisted']:
                            issue('read_count_mismatch', end=side, expected='blacklisted', reported=reads)
                        if side == 'five_prime':
                            expected_tags.add('LEAP_blacklisted')
                    elif length == 0:
                        if support != ['none']:
                            issue('support_mismatch', end=side, expected='none', reported=support)
                        if reads != ['0']:
                            issue('read_count_mismatch', end=side, expected='0', reported=reads)
                    elif length > 0:
                        if any(v not in ('reads', 'current_annotation') for v in support):
                            issue('support_mismatch', end=side, reported=support, expected='reads and/or current_annotation')
                        if any(v == 'blacklisted' or (s == 'reads' and int(v) == 0) for s, v in zip(support, reads)):
                            issue('read_count_mismatch', end=side, reported=reads, support=support)
                        if candidates:
                            row['candidate_provenance_checked'] = True
                            candidate = candidates.candidate(new, side)
                            endpoint = new['start'] if (side == 'five_prime') == (new['strand'] == '+') else new['end']
                            if (candidate is None or candidate['endpoint'] != endpoint
                                    or new['gene_id'] not in candidate['genes']
                                    or any(s != candidate['support'] for s in support)
                                    or any(v != str(candidate['read_count']) for v in reads)):
                                issue('candidate_mismatch', end=side, final_endpoint=endpoint,
                                      reported_support=support, reported_read_count=reads, candidate=candidate)
            # Tags also supply an independently checkable blacklist declaration.
            tags = new['leap_tag'].split(',') if new['leap_tag'] is not None else []
            row['observed_leap_tags'] = ','.join(tags)
            if old and old['strand'] == new['strand']:
                if not expected_tags:
                    expected_tags.add('LEAP_none')
                row['expected_leap_tags'] = ','.join(sorted(expected_tags))
                if len(tags) != len(set(tags)) or set(tags) != expected_tags:
                    issue('leap_tag_mismatch', observed=tags, expected=sorted(expected_tags))
        row['failed'] = any(row[flag] for flag in FLAGS)
        row['status'] = 'inconsistent' if row['failed'] else 'consistent' if new else 'missing_from_extended'
        row['details'] = compact(details)
        yield row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('current_annotation', help='Raw GFF3 (optionally .gz)')
    parser.add_argument('extended_annotation', help='Final GFF3 (optionally .gz)')
    parser.add_argument('output', help='CSV report filename')
    parser.add_argument('--extension-stats', required=True, help='Corresponding TSV file or directory, optionally gzipped')
    parser.add_argument('--leap-output-dir', help='Optional original run directory: verify provenance against fixed_extensions/')
    args = parser.parse_args()
    try:
        output = Path(args.output)
        sources = [Path(args.current_annotation), Path(args.extended_annotation)] + stats_paths(args.extension_stats)
        if args.leap_output_dir:
            sources.extend(p for p in (Path(args.leap_output_dir) / 'fixed_extensions').glob('*') if p.is_file())
        for source in sources:
            if output.resolve() == source.resolve() or (output.exists() and output.samefile(source)):
                raise ValueError('output must not overwrite an input annotation, stats or candidate file')
        raw = read_annotation(args.current_annotation)
        extended = read_annotation(args.extended_annotation, include_exons=True)
        stats = read_stats(args.extension_stats)
        candidates = CandidateContext(args.leap_output_dir, extended) if args.leap_output_dir else None
        counts, flags, failed = Counter(), Counter(), 0
        with output.open('w', newline='', encoding='utf-8') as handle:
            writer = csv.DictWriter(handle, fieldnames=FIELDS); writer.writeheader()
            for row in compare(raw, extended, stats, candidates):
                writer.writerow(row); counts[row['status']] += 1; failed += row['failed']
                for flag in FLAGS:
                    flags[flag] += row[flag]
        print('{}: {} failing transcript key(s)'.format('FAIL' if failed else 'PASS', failed))
        print('candidate_provenance: ' + ('enabled' if candidates else 'disabled; support/read counts checked for consistency only'))
        for status in STATUSES:
            print('{}: {}'.format(status, counts[status]))
        for flag in FLAGS:
            print('flag_{}: {}'.format(flag, flags[flag]))
        return int(bool(failed))
    except (OSError, EOFError, UnicodeError, ValueError, csv.Error) as error:
        print('ERROR: ' + str(error), file=sys.stderr)
        return 2


if __name__ == '__main__':
    sys.exit(main())
