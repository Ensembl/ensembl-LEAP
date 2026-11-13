#!/usr/bin/env python3
"""Check that extended GFF3 transcripts preserve exon counts and splice boundaries."""

import argparse
import csv
import gzip
import sys
from collections import Counter
from pathlib import Path

STATUSES = ('unchanged_internal_structure', 'changed_boundaries', 'exon_count_mismatch',
            'location_mismatch', 'missing_from_extended', 'missing_from_raw')
FLAGS = ('chromosome_mismatch', 'strand_mismatch', 'exon_count_mismatch', 'changed_boundaries')
FIELDS = ['transcript_id', 'status'] + [
    '{}_{}'.format(field, source) for source in ('current', 'extended')
    for field in ('chromosome', 'strand', 'exon_count', 'exon_coordinates')
] + list(FLAGS) + ['changed_boundary_details']


def read_annotation(path):
    """Load transcript and exon records, allowing arbitrary feature ordering."""
    transcripts, exons = {}, {}
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', encoding='utf-8') as handle:
        for number, line in enumerate(handle, 1):
            if line.strip() == '##FASTA':
                break
            if not line.strip() or line.startswith('#'):
                continue
            context = '{}:{}'.format(path, number)
            fields = line.rstrip('\r\n').split('\t')
            if len(fields) != 9:
                raise ValueError('{}: expected 9 GFF3 columns'.format(context))
            if fields[2] not in ('transcript', 'exon'):
                continue
            attributes = {}
            for item in fields[8].split(';'):
                if not item:
                    continue
                key, separator, value = item.partition('=')
                if not separator or not key or key in attributes:
                    raise ValueError('{}: malformed or duplicate attribute'.format(context))
                attributes[key] = value
            identifier = attributes.get('transcript_id', '')
            if not identifier.strip() or identifier == '.':
                raise ValueError('{}: missing transcript_id'.format(context))
            try:
                start, end = int(fields[3]), int(fields[4])
            except ValueError:
                raise ValueError('{}: coordinates must be integers'.format(context)) from None
            if start < 1 or end < start:
                raise ValueError('{}: invalid coordinates'.format(context))
            if not fields[0] or fields[0] == '.' or fields[6] not in ('+', '-', '.', '?'):
                raise ValueError('{}: invalid chromosome or strand'.format(context))
            record = dict(chromosome=fields[0], strand=fields[6], start=start, end=end,
                          context=context)
            if fields[2] == 'transcript':
                if identifier in transcripts:
                    raise ValueError('{}: duplicate transcript_id {}'.format(context, identifier))
                transcripts[identifier] = record
            else:
                exons.setdefault(identifier, []).append(record)
    for identifier, records in exons.items():
        if identifier not in transcripts:
            raise ValueError('{}: orphan exon for {}'.format(records[0]['context'], identifier))
    if not transcripts:
        raise ValueError('{}: no transcript records found'.format(path))
    for identifier, transcript in transcripts.items():
        records = exons.get(identifier, [])
        if not records:
            raise ValueError('{}: transcript {} has no exons'.format(transcript['context'], identifier))
        records.sort(key=lambda record: (record['start'], record['end']))
        previous_end = 0
        for record in records:
            if any(record[key] != transcript[key] for key in ('chromosome', 'strand')):
                raise ValueError('{}: inconsistent exon location for {}'.format(record['context'], identifier))
            if record['start'] < transcript['start'] or record['end'] > transcript['end']:
                raise ValueError('{}: exon outside transcript interval for {}'.format(record['context'], identifier))
            if record['start'] <= previous_end:
                raise ValueError('{}: overlapping or duplicate exons for {}'.format(record['context'], identifier))
            previous_end = record['end']
        transcript['exons'] = [(record['start'], record['end']) for record in records]
    return transcripts


def compare_transcripts(current, extended):
    for identifier in sorted(current.keys() | extended.keys()):
        raw, new = current.get(identifier), extended.get(identifier)
        row = dict(transcript_id=identifier)
        for source, record in (('current', raw), ('extended', new)):
            if record is not None:
                row.update({
                    'chromosome_' + source: record['chromosome'],
                    'strand_' + source: record['strand'],
                    'exon_count_' + source: len(record['exons']),
                    'exon_coordinates_' + source: ';'.join('{}-{}'.format(*exon) for exon in record['exons'])})
        if raw is None:
            row['status'] = 'missing_from_raw'
        elif new is None:
            row['status'] = 'missing_from_extended'
        else:
            row['chromosome_mismatch'] = raw['chromosome'] != new['chromosome']
            row['strand_mismatch'] = raw['strand'] != new['strand']
            row['exon_count_mismatch'] = len(raw['exons']) != len(new['exons'])
            details = []
            # Equal counts give an unambiguous genomic-position pairing, even when
            # another mismatch is also present. Never pair differing exon counts.
            if not row['exon_count_mismatch']:
                count = len(raw['exons'])
                for index, (old_exon, new_exon) in enumerate(zip(raw['exons'], new['exons'])):
                    for boundary, offset, protected in (('start', 0, index > 0),
                                                        ('end', 1, index < count - 1)):
                        if protected and old_exon[offset] != new_exon[offset]:
                            details.append('exon_{}_{}:{}->{}'.format(
                                index + 1, boundary, old_exon[offset], new_exon[offset]))
                row['changed_boundaries'] = bool(details)
                row['changed_boundary_details'] = ';'.join(details)
            if row['chromosome_mismatch'] or row['strand_mismatch']:
                row['status'] = 'location_mismatch'
            elif row['exon_count_mismatch']:
                row['status'] = 'exon_count_mismatch'
            elif details:
                row['status'] = 'changed_boundaries'
            else:
                row['status'] = 'unchanged_internal_structure'
        yield row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('current_annotation', help='Raw GFF3 annotation (optionally .gz)')
    parser.add_argument('extended_annotation', help='Extended GFF3 annotation (optionally .gz)')
    parser.add_argument('output', help='Comparison CSV output path')
    args = parser.parse_args()
    try:
        output = Path(args.output)
        for source in (args.current_annotation, args.extended_annotation):
            if output.resolve() == Path(source).resolve() or (output.exists() and output.samefile(source)):
                raise ValueError('output must not overwrite an input annotation')
        current = read_annotation(args.current_annotation)
        extended = read_annotation(args.extended_annotation)
        counts, flag_counts = Counter(), Counter()
        failed = 0
        with output.open('w', newline='', encoding='utf-8') as handle:
            writer = csv.DictWriter(handle, fieldnames=FIELDS)
            writer.writeheader()
            for row in compare_transcripts(current, extended):
                writer.writerow(row)
                counts[row['status']] += 1
                failed += any(row.get(flag) is True for flag in FLAGS)
                for flag in FLAGS:
                    flag_counts[flag] += row.get(flag) is True
        print('{}: {} failing transcript(s)'.format('FAIL' if failed else 'PASS', failed))
        for status in STATUSES:
            print('{}: {}'.format(status, counts[status]))
        for flag in FLAGS:
            print('flag_{}: {}'.format(flag, flag_counts[flag]))
        return 1 if failed else 0
    except (OSError, EOFError, UnicodeError, ValueError, csv.Error) as error:
        print('ERROR: {}'.format(error), file=sys.stderr)
        return 2


if __name__ == '__main__':
    sys.exit(main())
