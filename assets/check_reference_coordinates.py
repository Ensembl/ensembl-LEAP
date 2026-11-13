#!/usr/bin/env python3
"""Validate all GFF feature coordinates against explicit reference lengths."""

import argparse
import csv
import gzip
import math
import re
import shutil
import sys
import tempfile
import zlib
from collections import Counter
from pathlib import Path

FIELDS = ['status', 'source_path', 'line', 'feature_type', 'feature_id', 'transcript_id',
          'gene_id', 'chromosome', 'strand', 'start', 'end', 'reference_length', 'details']


def open_text(path):
    return (gzip.open if str(path).endswith('.gz') else open)(path, 'rt', encoding='utf-8', newline='')


def protect_output(output, inputs):
    for source in inputs:
        if output.resolve() == source.resolve() or (output.exists() and source.exists() and output.samefile(source)):
            raise ValueError('report would overwrite input: ' + str(source))


def annotation_paths(source):
    if not source.is_dir():
        if not source.is_file():
            raise ValueError(str(source) + ': expected a GFF file or directory')
        return [source]
    paths, identities = [], {}
    for path in sorted(source.iterdir()):
        if not path.is_file() or not path.name.endswith(('.gff', '.gff3', '.gff.gz', '.gff3.gz')):
            continue
        identity = re.sub(r'\.gff3?(\.gz)?$', '', path.name)
        if identity in identities:
            raise ValueError('ambiguous annotation copies: {} and {}'.format(identities[identity], path))
        identities[identity] = path
        paths.append(path)
    if not paths:
        raise ValueError(str(source) + ': no GFF/GFF3 files found')
    return paths


def reference_lengths(path):
    lengths, layout = {}, None
    with open_text(path) as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith('#'):
                continue
            context = '{}:{}'.format(path, number)
            fields = line.split()
            if len(fields) not in (2, 5) or (layout is not None and len(fields) != layout):
                raise ValueError(context + ': expected consistently two-column lengths or five-column FAI')
            layout = len(fields)
            if fields[0] == '.' or fields[0] in lengths:
                raise ValueError(context + ': invalid or duplicate reference sequence name')
            if not all(re.fullmatch(r'[0-9]+', v) for v in fields[1:]):
                raise ValueError(context + ': reference lengths/index fields must be integers')
            values = list(map(int, fields[1:]))
            if values[0] <= 0 or (layout == 5 and (values[2] <= 0 or values[3] < values[2])):
                raise ValueError(context + ': invalid reference length or FAI line layout')
            lengths[fields[0]] = values[0]
    if not lengths:
        raise ValueError(str(path) + ': no reference sequences found')
    return lengths


def parse_attributes(value, context):
    if value == '.':
        return {}
    if not value:
        raise ValueError(context + ': empty attributes column; use a dot for no attributes')
    result = {}
    for item in value.split(';'):
        if not item:
            continue
        key, sep, text = item.partition('=')
        if not sep or not key or key in result:
            raise ValueError(context + ': malformed or duplicate attribute')
        result[key] = text
    return result


def violations(path, lengths, counts):
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
            if not re.fullmatch(r'-?[0-9]+', f[3]) or not re.fullmatch(r'-?[0-9]+', f[4]):
                raise ValueError(context + ': coordinates must be integers')
            start, end = int(f[3]), int(f[4])
            if f[0] in ('', '.') or not f[1] or f[2] in ('', '.') or f[6] not in ('+', '-', '.', '?') or f[7] not in ('.', '0', '1', '2'):
                raise ValueError(context + ': malformed chromosome, source, feature, strand or phase')
            try:
                if f[5] != '.' and not math.isfinite(float(f[5])):
                    raise ValueError()
            except ValueError:
                raise ValueError(context + ': invalid score') from None
            attrs = parse_attributes(f[8], context)
            counts['features'] += 1
            length = lengths.get(f[0])
            issues = []
            if length is None:
                issues.append('unknown_chromosome')
            if start < 1:
                issues.append('start_before_reference')
            if end < start:
                issues.append('reversed_interval')
            if end < 1:
                issues.append('end_before_reference')
            if length is not None and (start > length or end > length):
                issues.append('beyond_reference_end')
            if issues:
                counts['violating_features'] += 1
                for issue in issues:
                    counts[issue] += 1
                yield dict(status=';'.join(issues), source_path=str(path), line=number,
                           feature_type=f[2], feature_id=attrs.get('ID', ''),
                           transcript_id=attrs.get('transcript_id', ''), gene_id=attrs.get('gene_id', ''),
                           chromosome=f[0], strand=f[6], start=start, end=end,
                           reference_length=length if length is not None else '',
                           details='Unknown exact sequence name' if length is None else
                           'Required: 1 <= start <= end <= {}'.format(length))


def run(args):
    source, reference, output = map(Path, (args.annotation, args.reference_lengths, args.output))
    paths = annotation_paths(source)
    # Guard unselected files too when a report is placed inside an input directory.
    inputs = paths + [reference]
    if source.is_dir():
        inputs += [p for p in source.iterdir() if p.is_file()]
    protect_output(output, inputs)
    lengths = reference_lengths(reference)
    counts = Counter(files=len(paths), features=0, violating_features=0)
    with tempfile.TemporaryFile(mode='w+', newline='', encoding='utf-8') as spool:
        writer = csv.DictWriter(spool, fieldnames=FIELDS)
        writer.writeheader()
        for path in paths:
            writer.writerows(violations(path, lengths, counts))
        spool.seek(0)
        with output.open('w', newline='', encoding='utf-8') as destination:
            shutil.copyfileobj(spool, destination)
    print(('FAIL' if counts['violating_features'] else 'PASS') + ': {} violating feature(s)'.format(counts['violating_features']))
    for name in ('files', 'features', 'violating_features'):
        print('{}: {}'.format(name, counts[name]))
    for name in sorted(counts.keys() - {'files', 'features', 'violating_features'}):
        print('{}: {}'.format(name, counts[name]))
    return int(bool(counts['violating_features']))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('annotation')
    parser.add_argument('reference_lengths')
    parser.add_argument('output')
    args = parser.parse_args()
    try:
        return run(args)
    except (OSError, ValueError, csv.Error, EOFError, zlib.error) as error:
        print('ERROR: ' + str(error), file=sys.stderr)
        return 2


if __name__ == '__main__':
    sys.exit(main())
