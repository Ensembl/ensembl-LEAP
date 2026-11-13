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
from urllib.parse import unquote


FIELDS = ['group_id', 'group_size', 'chromosome', 'start', 'end', 'strand',
          'transcript_id', 'gene_id', 'gene_name', 'transcript_name', 'gff_id',
          'exon_count', 'exon_coordinates', 'spliced_length', 'cds_count',
          'cds_coordinates', 'cds_phases', 'cds_length', 'jbrowse_label',
          'source_path', 'transcript_line', 'observed_tags', 'observed_leap_tags',
          'classification', 'is_mane', 'is_canonical_extended', 'is_blacklisted',
          'recommendation', 'reason', 'retained_member_ids', 'cross_gene_group']


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
                             gff_id=attrs.get('ID', ''), tags=tag_values(attrs.get('tag', '')),
                             leap_tags=tag_values(attrs.get('leap_tag', '')))
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


def tag_values(value):
    # Split the GFF3 list before decoding escaped commas inside individual values.
    return tuple(sorted({unquote(item).strip() for item in value.split(',') if item.strip()}))


def tag_classification(model):
    normalized = {re.sub(r'[ -]', '_', tag.casefold()) for tag in model['tags'] + model['leap_tags']}
    mane = bool(normalized & {'mane_select', 'mane_plus_clinical'})
    canonical = 'ensembl_canonical_extended' in normalized
    blacklisted = bool(normalized & {'leap_blacklisted', 'upstream_atg', 'cage_supported_tss', 'rp_supported_tis'})
    return dict(classification='mane' if mane else 'canonical_extended' if canonical else 'other',
                is_mane=mane, is_canonical_extended=canonical, is_blacklisted=blacklisted)


def numeric_stable_id(identifier):
    match = re.fullmatch(r'[A-Za-z_]+([0-9]+)(?:\.[0-9]+)?(?:_PAR_Y)?', identifier)
    return int(match[1]) if match else None


def recommendations(members):
    classes = {m['transcript_id']: tag_classification(m) for m in members}
    candidates = {m['transcript_id'] for m in members}
    result = {}
    mane = {tid for tid in candidates if classes[tid]['is_mane']}
    if mane:
        return {tid: ('keep', 'protected_MANE') if tid in mane else
                ('remove_from_gff', 'MANE_member_preferred') for tid in candidates}
    noncanonical = {tid for tid in candidates if not classes[tid]['is_canonical_extended']}
    canonical_only = not noncanonical
    if noncanonical:
        for tid in candidates - noncanonical:
            result[tid] = ('remove_from_gff', 'noncanonical_member_preferred')
        candidates = noncanonical
    blacklisted = {tid for tid in candidates if classes[tid]['is_blacklisted']}
    if blacklisted:
        for tid in candidates - blacklisted:
            result[tid] = ('remove_from_gff', 'blacklisted_member_preferred')
        candidates = blacklisted
    if len(candidates) == 1:
        winner = next(iter(candidates))
        result[winner] = ('keep', 'only_remaining_candidate')
    else:
        numbers = {tid: numeric_stable_id(tid) for tid in candidates}
        lowest = min((v for v in numbers.values() if v is not None), default=None)
        winners = {tid for tid, value in numbers.items() if value == lowest}
        if None in numbers.values() or len(winners) != 1:
            for tid in candidates:
                result[tid] = ('review', 'numeric_stable_id_cannot_resolve_candidates')
        else:
            winner = next(iter(winners))
            for tid in candidates:
                result[tid] = ('keep', 'lowest_numeric_stable_id') if tid == winner else (
                    'remove_from_gff', 'lower_numeric_stable_id_preferred')
    if canonical_only:
        result = {tid: (action, 'canonical_only_group:' + reason) for tid, (action, reason) in result.items()}
    return result


def report_rows(groups, recommend=False):
    for key, members in groups:
        # The label stays the same across reordered files and subset reports.
        group_id = 'DUP_' + hashlib.sha256(json.dumps(key, separators=(',', ':')).encode()).hexdigest()[:20]
        decisions = recommendations(members) if recommend else {}
        retained = ';'.join(sorted(tid for tid, (action, reason) in decisions.items() if action == 'keep'))
        cross_gene = len({m['gene_id'] for m in members}) > 1
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
            action, reason = decisions.get(model['transcript_id'], ('', ''))
            row.update(tag_classification(model), observed_tags=';'.join(model['tags']),
                       observed_leap_tags=';'.join(model['leap_tags']), recommendation=action,
                       reason=reason, retained_member_ids=retained, cross_gene_group=cross_gene)
            yield row


def run(args):
    source, output = Path(args.annotation), Path(args.output)
    paths = annotation_paths(source)
    protect_output(output, paths)
    # Finish validation before opening the report, preserving it on input errors.
    groups, transcripts = duplicate_groups(paths)
    action_counts = defaultdict(int)
    with output.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        for row in report_rows(groups, args.recommend_removals):
            writer.writerow(row)
            action_counts[row['recommendation']] += 1
    duplicate_transcripts = sum(len(members) for key, members in groups)
    print('{}: {} duplicate group(s), {} transcript(s) in groups'.format(
        ('REVIEW' if args.recommend_removals else 'FAIL') if groups else 'PASS', len(groups), duplicate_transcripts))
    for name, count in (('files', len(paths)), ('transcripts', transcripts),
                        ('duplicate_groups', len(groups)), ('duplicate_transcripts', duplicate_transcripts),
                        ('redundant_transcripts', duplicate_transcripts - len(groups)),
                        ('unique_transcripts', transcripts - duplicate_transcripts)):
        print('{}: {}'.format(name, count))
    if args.recommend_removals:
        for action in ('keep', 'remove_from_gff', 'review'):
            print('recommend_{}: {}'.format(action, action_counts[action]))
        print('Recommendations apply only to LEAP GFF models; no annotation or database records were changed.')
    return int(bool(groups))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('annotation', help='Pipeline GFF/GFF3 file or directory (including gzip)')
    parser.add_argument('output', help='CSV: one row per member of each duplicate group')
    parser.add_argument('--recommend-removals', action='store_true',
                        help='Recommend GFF-only retention/removal from tags and numeric stable IDs; never edit inputs')
    args = parser.parse_args()
    try:
        return run(args)
    except (OSError, ValueError, csv.Error, EOFError, zlib.error) as error:
        print('ERROR: ' + str(error), file=sys.stderr)
        return 2


if __name__ == '__main__':
    sys.exit(main())
