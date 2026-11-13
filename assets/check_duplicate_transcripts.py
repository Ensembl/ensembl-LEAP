#!/usr/bin/env python3
"""Find different transcript IDs with identical exon and CDS structures."""

import argparse
import csv
import gzip
import hashlib
import json
import math
import re
import sys
import zlib
from collections import defaultdict
from pathlib import Path


FIELDS = ['group_id', 'group_size', 'chromosome', 'start', 'end', 'strand',
          'transcript_id', 'gene_id', 'gene_name', 'transcript_name', 'gff_id',
          'exon_count', 'exon_coordinates', 'spliced_length', 'cds_count',
          'cds_coordinates', 'cds_phases', 'cds_length', 'jbrowse_label',
          'source_path', 'transcript_line']


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


def protect_output(output, inputs):
    for source in inputs:
        if output.resolve() == source.resolve() or (output.exists() and output.samefile(source)):
            raise ValueError('report would overwrite input: ' + str(source))


def parse_attributes(value, context):
    if value == '.':
        return {}
    attributes = {}
    if not value:
        raise ValueError(context + ': empty attributes column')
    for item in value.split(';'):
        if not item:
            continue
        key, separator, text = item.partition('=')
        if not separator or not key or key in attributes:
            raise ValueError(context + ': malformed or duplicate attribute')
        attributes[key] = text
    return attributes


def read_annotation(path):
    """Load one partition; children may occur before transcript records.

    Chromosome is part of model identity because X/Y PAR copies can share the
    same transcript_id. Repeated records for that ID on one chromosome are errors.
    """
    models = {}
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', encoding='utf-8') as handle:
        for number, line in enumerate(handle, 1):
            if line.strip() == '##FASTA':
                break
            if not line.strip() or line.startswith('#'):
                continue
            context = '{}:{}'.format(path, number)
            f = line.rstrip('\r\n').split('\t')
            if len(f) != 9:
                raise ValueError(context + ': expected 9 GFF3 columns')
            try:
                start, end = int(f[3]), int(f[4])
            except ValueError:
                raise ValueError(context + ': coordinates must be integers') from None
            if not 1 <= start <= end:
                raise ValueError(context + ': invalid coordinates')
            if (f[0] in ('', '.') or not f[1] or f[2] in ('', '.')
                    or f[6] not in ('+', '-', '.', '?') or f[7] not in ('.', '0', '1', '2')):
                raise ValueError(context + ': malformed chromosome, source, feature, strand or phase')
            try:
                if f[5] != '.' and not math.isfinite(float(f[5])):
                    raise ValueError()
            except ValueError:
                raise ValueError(context + ': invalid score') from None
            attrs = parse_attributes(f[8], context)
            if f[2] not in ('transcript', 'exon', 'CDS'):
                continue
            identifier = attrs.get('transcript_id', '')
            if not identifier.strip() or identifier == '.':
                raise ValueError(context + ': missing transcript_id')
            if f[6] not in ('+', '-'):
                raise ValueError(context + ': transcript/exon/CDS strand must be + or -')
            key = (f[0], identifier)
            if key not in models:
                models[key] = dict(chromosome=f[0], strand=f[6], transcript_id=identifier,
                                   exons=[], cds=[], children=[], record=None, context=context)
            model = models[key]
            if model['strand'] != f[6]:
                raise ValueError(context + ': inconsistent strand for ' + identifier)
            if f[2] == 'transcript':
                if model['record'] is not None:
                    raise ValueError(context + ': duplicate transcript record for ' + identifier)
                model.update(record=(start, end), source_path=str(path.resolve()),
                             transcript_line=number, gene_id=attrs.get('gene_id', ''),
                             gene_name=attrs.get('gene_name', ''),
                             transcript_name=attrs.get('transcript_name', ''),
                             gff_id=attrs.get('ID', ''))
            else:
                model['children'].append((start, end, attrs.get('gene_id', ''), context))
                if f[2] == 'exon':
                    model['exons'].append((start, end))
                else:
                    model['cds'].append((start, end, f[7]))
    for model in models.values():
        context, identifier = model['context'], model['transcript_id']
        if model['record'] is None:
            raise ValueError(context + ': orphan exon/CDS for ' + identifier)
        if not model['exons']:
            raise ValueError(context + ': transcript has no exons: ' + identifier)
        model['exons'] = tuple(sorted(model['exons']))
        model['cds'] = tuple(sorted(model['cds']))
        if model['record'] != (model['exons'][0][0], model['exons'][-1][1]):
            raise ValueError(context + ': transcript bounds do not match exon span for ' + identifier)
        for intervals, label in ((model['exons'], 'exon'), (model['cds'], 'CDS')):
            for first, second in zip(intervals, intervals[1:]):
                if second[0] <= first[1]:
                    raise ValueError(context + ': overlapping or duplicate {} intervals for {}'.format(label, identifier))
        for start, end, gene, child_context in model.pop('children'):
            if gene and model['gene_id'] and gene != model['gene_id']:
                raise ValueError(child_context + ': inconsistent gene_id for ' + identifier)
            if not model['record'][0] <= start <= end <= model['record'][1]:
                raise ValueError(child_context + ': child outside transcript interval for ' + identifier)
        for start, end, phase in model['cds']:
            if not any(left <= start <= end <= right for left, right in model['exons']):
                raise ValueError(context + ': CDS outside exons for ' + identifier)
    return models


def structure_key(model):
    return (model['chromosome'], model['strand'], model['exons'], model['cds'])


def duplicate_groups(paths):
    groups, seen = defaultdict(list), {}
    for path in paths:
        for identity, model in read_annotation(path).items():
            if identity in seen:
                raise ValueError('{}: transcript {} on {} also occurs in {}'.format(
                    model['context'], identity[1], identity[0], seen[identity]))
            seen[identity] = model['context']
            groups[structure_key(model)].append(model)
    duplicates = [(key, members) for key, members in sorted(groups.items()) if len(members) > 1]
    return duplicates, len(seen)


def report_rows(groups):
    for key, members in groups:
        # The label stays the same across reordered files and subset reports.
        group_id = 'DUP_' + hashlib.sha256(json.dumps(key, separators=(',', ':')).encode()).hexdigest()[:20]
        for model in sorted(members, key=lambda item: item['transcript_id']):
            exons, cds = model['exons'], model['cds']
            start, end = model['record']
            row = {name: model[name] for name in ('chromosome', 'strand', 'transcript_id',
                'gene_id', 'gene_name', 'transcript_name', 'gff_id', 'source_path', 'transcript_line')}
            row.update(group_id=group_id, group_size=len(members), start=start, end=end,
                       exon_count=len(exons), exon_coordinates=';'.join('{}-{}'.format(*e) for e in exons),
                       spliced_length=sum(right - left + 1 for left, right in exons),
                       cds_count=len(cds), cds_coordinates=';'.join('{}-{}'.format(*c[:2]) for c in cds),
                       cds_phases=';'.join(c[2] for c in cds),
                       cds_length=sum(right - left + 1 for left, right, phase in cds),
                       jbrowse_label='{}:{}-{}'.format(model['chromosome'], start, end))
            yield row


def run(args):
    source, output = Path(args.annotation), Path(args.output)
    paths = annotation_paths(source)
    protect_output(output, paths)
    # Finish validation before opening the report, preserving it on input errors.
    groups, transcripts = duplicate_groups(paths)
    with output.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(report_rows(groups))
    duplicate_transcripts = sum(len(members) for key, members in groups)
    print('{}: {} duplicate group(s), {} transcript(s) in groups'.format(
        'FAIL' if groups else 'PASS', len(groups), duplicate_transcripts))
    for name, count in (('files', len(paths)), ('transcripts', transcripts),
                        ('duplicate_groups', len(groups)), ('duplicate_transcripts', duplicate_transcripts),
                        ('redundant_transcripts', duplicate_transcripts - len(groups)),
                        ('unique_transcripts', transcripts - duplicate_transcripts)):
        print('{}: {}'.format(name, count))
    return int(bool(groups))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('annotation', help='Pipeline GFF/GFF3 file or directory (including gzip)')
    parser.add_argument('output', help='CSV: one row per member of each duplicate group')
    args = parser.parse_args()
    try:
        return run(args)
    except (OSError, ValueError, csv.Error, EOFError, zlib.error) as error:
        print('ERROR: ' + str(error), file=sys.stderr)
        return 2


if __name__ == '__main__':
    sys.exit(main())
