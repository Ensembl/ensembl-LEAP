#!/usr/bin/env python3
"""Check that LEAP preserves every transcript feature except terminal extensions."""

import argparse
import csv
import gzip
import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import NamedTuple


IGNORED_ATTRIBUTES = {'leap_tag', 'transcript_class'}
EXON_FEATURES = {'CDS', 'start_codon', 'stop_codon', 'five_prime_UTR', 'three_prime_UTR'}
STATUSES = ('unchanged', 'extended', 'changed_features', 'missing_from_extended', 'missing_from_raw')
FLAGS = ('location_mismatch', 'shortened', 'feature_count_mismatch',
         'inconsistent_extension', 'feature_records_differ')
FIELDS = ['transcript_id', 'status', 'failed'] + [
    source + '_' + field for source in ('current', 'extended')
    for field in ('chromosome', 'strand', 'start', 'end', 'feature_counts')
] + list(FLAGS) + ['extension_details', 'feature_differences']


class Feature(NamedTuple):
    chromosome: str
    source: str
    kind: str
    start: int
    end: int
    score: str
    strand: str
    phase: str
    attributes: tuple


def compact(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'))


def read_annotation(path):
    """Read all transcript-associated features, without relying on their order."""
    models = {}
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
                raise ValueError(context + ': expected 9 GFF3 columns')
            if fields[2] == 'gene':
                continue
            attributes = {}
            for item in fields[8].split(';'):
                if not item:
                    continue
                key, separator, value = item.partition('=')
                if not separator or not key or key in attributes:
                    raise ValueError(context + ': malformed or duplicate attribute')
                attributes[sys.intern(key)] = sys.intern(value)
            for key in ('transcript_id', 'gene_id'):
                if not attributes.get(key, '').strip() or attributes[key] == '.':
                    raise ValueError(context + ': missing ' + key)
            try:
                start, end = int(fields[3]), int(fields[4])
            except ValueError:
                raise ValueError(context + ': coordinates must be integers') from None
            if not 1 <= start <= end:
                raise ValueError(context + ': invalid coordinates')
            if (fields[0] in ('', '.') or fields[6] not in ('+', '-')
                    or fields[1] == '' or fields[2] in ('', '.')):
                raise ValueError(context + ': invalid chromosome, source, feature type or strand')
            if fields[7] not in ('.', '0', '1', '2'):
                raise ValueError(context + ': invalid phase')
            if fields[5] != '.':
                try:
                    if not math.isfinite(float(fields[5])):
                        raise ValueError()
                except ValueError:
                    raise ValueError(context + ': invalid score') from None
            feature = Feature(sys.intern(fields[0]), sys.intern(fields[1]), sys.intern(fields[2]),
                              start, end, sys.intern(fields[5]), fields[6], fields[7],
                              tuple(sorted((k, v) for k, v in attributes.items()
                                           if k not in IGNORED_ATTRIBUTES)))
            identifier = attributes['transcript_id']
            model = models.setdefault(identifier, dict(transcript=None, features=set(),
                                                       context=context, gene=attributes['gene_id']))
            if feature.kind == 'transcript':
                if model['transcript'] is not None:
                    raise ValueError(context + ': duplicate transcript record for ' + identifier)
                model['transcript'] = feature
            if feature in model['features']:
                raise ValueError(context + ': duplicate feature record for ' + identifier)
            model['features'].add(feature)
            if model['gene'] != attributes['gene_id']:
                raise ValueError(context + ': inconsistent gene_id for ' + identifier)
    if not models:
        raise ValueError(str(path) + ': no transcript records found')
    for identifier, model in models.items():
        transcript = model['transcript']
        context = model['context']
        if transcript is None:
            raise ValueError(context + ': orphan child records for ' + identifier)
        exons = sorted((f for f in model['features'] if f.kind == 'exon'),
                       key=lambda f: (f.start, f.end))
        if not exons:
            raise ValueError(context + ': transcript ' + identifier + ' has no exons')
        parent_ids = {dict(f.attributes).get('ID') for f in model['features']} | {identifier}
        for previous, following in zip(exons, exons[1:]):
            if previous.end >= following.start:
                raise ValueError(context + ': overlapping or duplicate exon intervals for ' + identifier)
        for feature in model['features']:
            parent = dict(feature.attributes).get('Parent')
            if feature.kind != 'transcript' and parent and any(
                    value not in parent_ids for value in parent.split(',')):
                raise ValueError(context + ': orphan Parent reference for ' + identifier)
            if (feature.chromosome, feature.strand) != (transcript.chromosome, transcript.strand):
                raise ValueError(context + ': inconsistent feature location for ' + identifier)
            if feature.start < transcript.start or feature.end > transcript.end:
                raise ValueError(context + ': feature outside transcript interval for ' + identifier)
            if feature.kind in EXON_FEATURES and not any(
                    exon.start <= feature.start <= feature.end <= exon.end for exon in exons):
                raise ValueError(context + ': ' + feature.kind + ' lacks a containing exon for ' + identifier)
        model['exons'] = exons
    return models


def terminal_utr(strand, side):
    return 'five_prime_UTR' if (side == 'left') == (strand == '+') else 'three_prime_UTR'


def normalized(feature, original, extended):
    """Remove only additions beyond the raw transcript's two outer boundaries."""
    old, new = original['transcript'], extended['transcript']
    start, end = feature.start, feature.end
    for side in ('left', 'right'):
        left = side == 'left'
        old_outer, new_outer = (old.start, new.start) if left else (old.end, new.end)
        outer = feature.start if left else feature.end
        if not (outer < old_outer if left else outer > old_outer):
            continue
        eligible = feature.kind == 'transcript'
        if feature.kind == 'exon':
            eligible = feature == extended['exons'][0 if left else -1]
        if feature.kind == terminal_utr(old.strand, side) and outer == new_outer:
            eligible = any(f.kind == feature.kind and (f.start if left else f.end) == old_outer
                           for f in original['features'])
        if eligible:
            if left:
                start = old_outer
            else:
                end = old_outer
    return feature._replace(start=start, end=end)


def record_details(feature, count):
    return dict(feature_type=feature.kind, chromosome=feature.chromosome,
                start=feature.start, end=feature.end, strand=feature.strand,
                source=feature.source, score=feature.score, phase=feature.phase,
                attributes=dict(feature.attributes), count=count)


def extension_consistency(original, extended):
    details = []
    old, new = original['transcript'], extended['transcript']
    for label, model in (('current', original), ('extended', extended)):
        transcript = model['transcript']
        if (transcript.start, transcript.end) != (model['exons'][0].start, model['exons'][-1].end):
            details.append(dict(issue='transcript_exon_span_mismatch', annotation=label))
    for side in ('left', 'right'):
        left = side == 'left'
        old_outer, new_outer = (old.start, new.start) if left else (old.end, new.end)
        if old_outer == new_outer:
            continue
        kind = terminal_utr(old.strand, side)
        expected = sum(f.kind == kind and (f.start if left else f.end) == old_outer
                       for f in original['features'])
        actual = sum(f.kind == kind and (f.start if left else f.end) == new_outer
                     for f in extended['features'])
        if expected and actual != expected:
            details.append(dict(issue='terminal_utr_does_not_follow_extension', side=side,
                                feature_type=kind, expected=expected, observed=actual,
                                current_boundary=old_outer, extended_boundary=new_outer))
    return details


def compare_transcripts(current, extended):
    for identifier in sorted(current.keys() | extended.keys()):
        original, new = current.get(identifier), extended.get(identifier)
        row = dict(transcript_id=identifier)
        for label, model in (('current', original), ('extended', new)):
            if model is not None:
                transcript = model['transcript']
                for field in ('chromosome', 'strand', 'start', 'end'):
                    row[label + '_' + field] = getattr(transcript, field)
                row[label + '_feature_counts'] = compact(Counter(f.kind for f in model['features']))
        if original is None:
            row.update(status='missing_from_raw', failed=True)
        elif new is None:
            row.update(status='missing_from_extended', failed=False)
        else:
            old_tx, new_tx = original['transcript'], new['transcript']
            row['location_mismatch'] = (old_tx.chromosome, old_tx.strand) != (new_tx.chromosome, new_tx.strand)
            row['shortened'] = new_tx.start > old_tx.start or new_tx.end < old_tx.end
            row['feature_count_mismatch'] = Counter(f.kind for f in original['features']) != Counter(f.kind for f in new['features'])
            extension_details = extension_consistency(original, new)
            row['inconsistent_extension'] = bool(extension_details)
            before = Counter(original['features'])
            after = Counter(normalized(f, original, new) for f in new['features'])
            row['feature_records_differ'] = before != after
            row['failed'] = any(row[flag] for flag in FLAGS)
            row['status'] = ('changed_features' if row['failed'] else
                             'extended' if (old_tx.start, old_tx.end) != (new_tx.start, new_tx.end)
                             else 'unchanged')
            row['extension_details'] = compact(extension_details)
            if before != after:
                kinds = {f.kind for f in (before - after).keys() | (after - before).keys()}
                row['feature_differences'] = compact(dict(
                    current_only=[record_details(f, n) for f, n in sorted((before - after).items())],
                    normalized_extended_only=[record_details(f, n) for f, n in sorted((after - before).items())],
                    actual_extended_features=[record_details(f, 1) for f in sorted(new['features'])
                                              if f.kind in kinds]))
            else:
                row['feature_differences'] = '{}'
        yield row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('current_annotation', help='Raw pipeline GFF3 (optionally .gz)')
    parser.add_argument('extended_annotation', help='Extended pipeline GFF3 (optionally .gz)')
    parser.add_argument('output', help='Comparison CSV output path')
    args = parser.parse_args()
    try:
        output = Path(args.output)
        for source in (args.current_annotation, args.extended_annotation):
            if output.resolve() == Path(source).resolve() or (output.exists() and output.samefile(source)):
                raise ValueError('output must not overwrite an input annotation')
        current = read_annotation(args.current_annotation)
        extended = read_annotation(args.extended_annotation)
        counts, flags, failed = Counter(), Counter(), 0
        with output.open('w', encoding='utf-8', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=FIELDS)
            writer.writeheader()
            for row in compare_transcripts(current, extended):
                writer.writerow(row)
                counts[row['status']] += 1
                failed += row['failed']
                for flag in FLAGS:
                    flags[flag] += row.get(flag) is True
        print('{}: {} failing transcript(s)'.format('FAIL' if failed else 'PASS', failed))
        for status in STATUSES:
            print('{}: {}'.format(status, counts[status]))
        for flag in FLAGS:
            print('flag_{}: {}'.format(flag, flags[flag]))
        return int(bool(failed))
    except (OSError, EOFError, UnicodeError, ValueError, csv.Error) as error:
        print('ERROR: {}'.format(error), file=sys.stderr)
        return 2


if __name__ == '__main__':
    sys.exit(main())
