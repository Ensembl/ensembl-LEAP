#!/usr/bin/env python3
"""Restore internal terminal-exon extensions using LEAP2GO's CDS-extreme policy.

Standalone LEAP asset, also reused by manual-review modules. Inputs are never modified.
"""
import argparse
from collections import Counter
import csv
import gzip
import json
import logging
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile

STATS_FIELDS = ['gene_id', 'transcript_id', 'transcript_location', 'five_prime_extension',
                'three_prime_extension', 'five_prime_support', 'five_prime_read_count',
                'three_prime_support', 'three_prime_read_count']


def stem(path):
    return re.sub(r'\.gff3?(\.gz)?$', '', path.name)


def records(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', encoding='utf-8') as handle:
        for number, text in enumerate(handle, 1):
            if text.strip() == '##FASTA':
                break
            if not text.strip() or text.startswith('#'):
                continue
            fields = text.rstrip('\r\n').split('\t')
            context = '{}:{}'.format(path, number)
            if len(fields) != 9:
                raise ValueError(context + ': expected 9 GFF3 columns')
            try:
                start, end = int(fields[3]), int(fields[4])
            except ValueError:
                raise ValueError(context + ': coordinates must be integers') from None
            if not 1 <= start <= end:
                raise ValueError(context + ': invalid coordinates')
            attrs = {}
            for item in fields[8].split(';'):
                if not item:
                    continue
                key, sep, value = item.partition('=')
                if not sep or not key or key in attrs:
                    raise ValueError(context + ': malformed or duplicate attribute')
                attrs[key] = value
            yield fields, attrs, context


def exon_number(attrs, context):
    try:
        value = int(attrs['exon_number'])
        if value < 1:
            raise ValueError()
        return value
    except (KeyError, ValueError):
        raise ValueError(context + ': a positive exon_number is required') from None


def scan_extended(path):
    models, extremes = {}, {}
    for fields, attrs, context in records(path):
        identifier = attrs.get('transcript_id')
        if not identifier or fields[2] == 'gene':
            continue
        if not attrs.get('gene_id') or fields[6] not in ('+', '-'):
            raise ValueError(context + ': gene_id and explicit strand required')
        key = (fields[0], identifier)
        location = (attrs['gene_id'], fields[6])
        model = models.setdefault(key, dict(location=location, transcript=None, exons={}, cds=[]))
        if model['location'] != location:
            raise ValueError(context + ': inconsistent gene or strand for ' + identifier)
        interval = (int(fields[3]), int(fields[4]))
        if fields[2] == 'transcript':
            if model['transcript'] is not None:
                raise ValueError(context + ': duplicate transcript ' + identifier)
            model['transcript'] = (fields, attrs)
        elif fields[2] == 'exon':
            number = exon_number(attrs, context)
            if number in model['exons']:
                raise ValueError(context + ': duplicate exon_number for ' + identifier)
            model['exons'][number] = interval
        elif fields[2] == 'CDS':
            model['cds'].append(interval)
            gene_key = (fields[0], attrs['gene_id'])
            # Preserve LEAP2GO's first-seen tie behavior and boundary-OR rule.
            if gene_key not in extremes:
                extremes[gene_key] = [interval, interval]
            else:
                if interval[0] < extremes[gene_key][0][0]:
                    extremes[gene_key][0] = interval
                if interval[1] > extremes[gene_key][1][1]:
                    extremes[gene_key][1] = interval
    if not models:
        raise ValueError(str(path) + ': no transcript-associated records')
    retained = {}
    dropped = []
    for key, model in models.items():
        if not model['cds']:
            dropped.append(key)
            continue
        if model['transcript'] is None or not model['exons']:
            raise ValueError('{}: coding transcript {} lacks transcript/exon records'.format(path, key[1]))
        low, high = extremes[(key[0], model['location'][0])]
        model['restore_left'] = not any(a == low[0] or b == low[1] for a, b in model['cds'])
        model['restore_right'] = not any(a == high[0] or b == high[1] for a, b in model['cds'])
        numbers = sorted(model['exons'])
        if model['location'][1] == '-':
            numbers.reverse()
        model['replace'] = set()
        if model['restore_left']:
            model['replace'].add(numbers[0])
        if model['restore_right']:
            model['replace'].add(numbers[-1])
        classes = [label for flag, label in [('restore_left', 'internal_left'), ('restore_right', 'internal_right')] if model[flag]]
        model['class'] = ','.join(classes) if classes else 'external'
        retained[key] = model
    return retained, dropped


def scan_raw(path, models):
    raw = {key: dict(exons={}, transcript=None) for key in models}
    replacements = []
    for fields, attrs, context in records(path):
        key = (fields[0], attrs.get('transcript_id'))
        if key not in models or fields[2] == 'gene':
            continue
        model = models[key]
        if (attrs.get('gene_id'), fields[6]) != model['location']:
            raise ValueError(context + ': raw/extended gene or strand differs')
        if fields[2] == 'transcript':
            if raw[key]['transcript'] is not None:
                raise ValueError(context + ': duplicate raw transcript')
            raw[key]['transcript'] = fields
        if fields[2] == 'exon':
            number = exon_number(attrs, context)
            if number in raw[key]['exons']:
                raise ValueError(context + ': duplicate raw exon_number')
            raw[key]['exons'][number] = (int(fields[3]), int(fields[4]))
        if fields[2] != 'transcript' and attrs.get('exon_number') and exon_number(attrs, context) in model['replace']:
            replacements.append(fields)
    for key, model in models.items():
        original = raw[key]
        if original['transcript'] is None or original['exons'].keys() != model['exons'].keys():
            raise ValueError('{}: missing raw transcript or different exon numbers for {}'.format(path, key[1]))
        all_exons = sorted(model['exons'].items(), key=lambda item: item[1])
        for number in model['exons']:
            if number in model['replace']:
                continue
            old, new = original['exons'][number], model['exons'][number]
            # Only terminal outer boundaries may differ before processing.
            if ((number != all_exons[0][0] and old[0] != new[0])
                    or (number != all_exons[-1][0] and old[1] != new[1])):
                raise ValueError('{}: protected exon boundary differs for {}'.format(path, key[1]))
    return raw, replacements


def serialize(fields, attrs=None):
    fields = list(fields)
    if attrs is not None:
        fields[8] = ';'.join(key + '=' + value for key, value in attrs.items())
    return '\t'.join(fields) + '\n'


def load_stats(path):
    result = {}
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', encoding='utf-8', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if not set(STATS_FIELDS) <= set(reader.fieldnames or []):
            raise ValueError(str(path) + ': missing extension stats columns')
        for row in reader:
            if any(row.get(key) is None for key in STATS_FIELDS):
                raise ValueError(str(path) + ': incomplete stats row')
            identifier = row['transcript_id']
            if identifier in result:
                raise ValueError(str(path) + ': duplicate stats transcript_id ' + identifier)
            result[identifier] = row
    return result


def prepare_models(models, raw, stats):
    for key, model in models.items():
        identifier = key[1]
        if identifier not in stats:
            raise ValueError('Missing extension stats for ' + identifier)
        row = dict(stats[identifier])
        original = raw[key]
        before, attrs = model['transcript']
        if row['gene_id'] != model['location'][0] or row['transcript_location'] != '{}:{}-{}'.format(key[0], before[3], before[4]):
            raise ValueError('Input stats gene/location differs for ' + identifier)
        final_exons = {number: original['exons'][number] if number in model['replace'] else interval for number, interval in model['exons'].items()}
        raw_start = min(a for a, b in original['exons'].values())
        raw_end = max(b for a, b in original['exons'].values())
        start = min(a for a, b in final_exons.values())
        end = max(b for a, b in final_exons.values())
        if start > raw_start or end < raw_end:
            raise ValueError('Processed transcript is shorter than raw: ' + identifier)
        tags = attrs.get('leap_tag', 'LEAP_none').split(',')
        for side, final, previous, baseline, restored in [
                ('left', start, int(before[3]), raw_start, model['restore_left']),
                ('right', end, int(before[4]), raw_end, model['restore_right'])]:
            prime = 'five_prime' if (side == 'left') == (model['location'][1] == '+') else 'three_prime'
            old_length, length = abs(previous - baseline), abs(final - baseline)
            try:
                if int(row[prime + '_extension']) != old_length:
                    raise ValueError()
            except ValueError:
                raise ValueError('Input stats {} length differs for {}'.format(prime, identifier)) from None
            row[prime + '_extension'] = str(length)
            if restored or final != previous:
                if length:
                    raise ValueError('Cannot infer support for changed extension: ' + identifier)
                # Internal/external classification does not change an existing
                # 5' blacklist declaration, even when that end is restored.
                if not (prime == 'five_prime' and row[prime + '_support'] == 'blacklisted'):
                    row[prime + '_support'], row[prime + '_read_count'] = 'none', '0'
            if not length:
                tag = 'LEAP_5' if prime == 'five_prime' else 'LEAP_3'
                tags = [value for value in tags if value != tag]
        row['transcript_location'] = '{}:{}-{}'.format(key[0], start, end)
        attrs = dict(attrs, leap_tag=','.join(tags) if tags else 'LEAP_none', transcript_class=model['class'])
        fields = list(before)
        fields[3:5] = [str(start), str(end)]
        model['output_transcript'] = serialize(fields, attrs)
        model['output_stats'] = row


def sorted_copy(source, destination):
    with destination.open('w', encoding='utf-8') as handle:
        subprocess.run(['sort', '-s', '-k1,1', '-k4,4n', '-k5,5n', str(source)], stdout=handle,
                       env=dict(os.environ, LC_ALL='C'), check=True)
    source.unlink()


def process_file(extended, raw_path, stats_path, folder):
    models, dropped = scan_extended(extended)
    raw, replacements = scan_raw(raw_path, models)
    prepare_models(models, raw, load_stats(stats_path))
    name = stem(extended)
    temporary = folder / (name + '.unsorted.gff')
    potential = folder / (name + '.potential.unsorted.gff')
    with temporary.open('w') as output, potential.open('w') as internal:
        output.write('##gff-version 3\n')
        internal.write('##gff-version 3\n')
        for fields, attrs, context in records(extended):
            key = (fields[0], attrs.get('transcript_id'))
            if key not in models or fields[2] == 'gene':
                continue
            model = models[key]
            if model['replace']:
                internal.write(serialize(fields))
            if fields[2] == 'transcript':
                output.write(model['output_transcript'])
            elif not attrs.get('exon_number') or exon_number(attrs, context) not in model['replace']:
                if model['class'] == 'external':
                    attrs['transcript_class'] = 'external'
                output.write(serialize(fields, attrs))
        for fields in replacements:
            output.write(serialize(fields))
    sorted_copy(temporary, folder / 'gff' / (name + '.gff'))
    sorted_copy(potential, folder / 'potential_internal' / (name + '.gff'))
    with (folder / 'stats' / re.sub(r'\.gz$', '', stats_path.name)).open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=STATS_FIELDS, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(model['output_stats'] for model in models.values())
    counts = Counter()
    for label in ('external', 'internal_left', 'internal_right', 'internal_both'):
        with (folder / 'ids' / (name + '_' + label + '.txt')).open('w') as handle:
            for key, model in models.items():
                category = 'internal_both' if model['restore_left'] and model['restore_right'] else model['class']
                if category == label:
                    handle.write(key[1].split('.')[0] + '\n')
                    counts[label] += 1
    with (folder / 'ids' / (name + '_classification.tsv')).open('w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['chromosome', 'transcript_id', 'classification'])
        writer.writerows((key[0], key[1], model['class']) for key, model in models.items())
        writer.writerows((key[0], key[1], 'dropped_no_cds') for key in dropped)
    counts['dropped_no_cds'] = len(dropped)
    return dict(input=str(extended), raw=str(raw_path), stats=str(stats_path), counts=dict(counts))


def run(input_dir, output_dir):
    input_dir, output_dir = input_dir.resolve(), output_dir.resolve()
    if not shutil.which('sort'):
        raise ValueError('Required command not found: sort')
    if output_dir.exists():
        raise ValueError('Output already exists; choose a new directory: ' + str(output_dir))
    files = sorted(p for p in (input_dir / 'extended_transcripts').iterdir() if re.search(r'\.gff3?(\.gz)?$', p.name))
    if not files:
        raise ValueError('No extended GFF3 files found')
    if len({stem(path) for path in files}) != len(files):
        raise ValueError('Duplicate extended annotation basenames would overwrite outputs')
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.' + output_dir.name + '-build-', dir=output_dir.parent) as directory:
        staging = Path(directory) / 'bundle'
        folder = staging / 'external_extensions'
        for name in ('gff', 'stats', 'ids', 'potential_internal'):
            (folder / name).mkdir(parents=True, exist_ok=True)
        logger = logging.getLogger('internal_extensions')
        logger.setLevel(logging.INFO)
        handler = logging.FileHandler(staging / 'processing.log')
        logger.addHandler(handler)
        results = []
        try:
            for extended in files:
                name = stem(extended)
                match = re.fullmatch(r'(.+)_(positive|negative)_extended', name)
                if not match:
                    raise ValueError('Expected <chromosome>_<positive|negative>_extended.gff: ' + extended.name)
                chromosome = match[1]
                candidates = [input_dir / 'raw_annotation' / ('current_annotation_' + chromosome + suffix) for suffix in ('.gff3', '.gff', '.gff3.gz', '.gff.gz')]
                raw_files = [path for path in candidates if path.is_file()]
                if len(raw_files) != 1:
                    raise ValueError('Expected exactly one raw annotation for ' + chromosome)
                candidates = [input_dir / 'extended_transcripts_stats' / (chromosome + '_' + match[2] + '_stats.tsv' + suffix) for suffix in ('', '.gz')]
                stat_files = [path for path in candidates if path.is_file()]
                if len(stat_files) != 1:
                    raise ValueError('Expected exactly one stats TSV for ' + name)
                result = process_file(extended, raw_files[0], stat_files[0], folder)
                results.append(result)
                message = '{}: {}'.format(extended.name, result['counts'])
                logger.info(message)
                print(message, flush=True)
            totals = Counter()
            for result in results:
                totals.update(result['counts'])
            (staging / 'processing-summary.json').write_text(json.dumps(dict(input_directory=str(input_dir), totals=dict(totals), files=results), indent=2) + '\n')
        finally:
            logger.removeHandler(handler)
            handler.close()
        if output_dir.exists():
            raise ValueError('Output appeared during processing: ' + str(output_dir))
        staging.rename(output_dir)
    print('READY: ' + str(output_dir / 'external_extensions'))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('leap_output_dir', type=Path)
    parser.add_argument('output_dir', type=Path, help='New output directory; must not exist')
    args = parser.parse_args()
    try:
        run(args.leap_output_dir, args.output_dir)
        return 0
    except (OSError, EOFError, UnicodeError, ValueError, csv.Error, subprocess.SubprocessError) as error:
        print('ERROR: ' + str(error), file=sys.stderr)
        return 2


if __name__ == '__main__':
    sys.exit(main())
