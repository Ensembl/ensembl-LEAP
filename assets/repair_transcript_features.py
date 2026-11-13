#!/usr/bin/env python3
"""Rebuild retained LEAP models from raw features, without selecting new extensions."""

import argparse
from collections import Counter, defaultdict
import csv
import gzip
import json
import math
from pathlib import Path
import sys
import tempfile
from typing import NamedTuple

sys.dont_write_bytecode = True
# These are the independent, unchanged datachecks used to validate the result.
import check_transcript_preservation as preservation
import check_stats_gff_consistency as consistency
import check_gff_structural_integrity as structural

IGNORED = {'leap_tag', 'transcript_class'}
FIELDS = ['chromosome', 'transcript_id', 'status', 'original_start', 'original_end',
          'repaired_start', 'repaired_end', 'removed_records', 'added_records',
          'removed_features', 'added_features', 'stats_changes', 'generated_attributes', 'message']


class Record(NamedTuple):
    fields: tuple
    attrs: dict


def compact(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'))


def read_records(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', encoding='utf-8') as handle:
        for number, line in enumerate(handle, 1):
            if line.strip() == '##FASTA': break
            if not line.strip() or line.startswith('#'): continue
            f = line.rstrip('\r\n').split('\t')
            context = '{}:{}'.format(path, number)
            if len(f) != 9: raise ValueError(context + ': expected nine GFF3 columns')
            try:
                start, end = int(f[3]), int(f[4])
                valid_score = f[5] == '.' or math.isfinite(float(f[5]))
            except ValueError:
                raise ValueError(context + ': invalid coordinates or score') from None
            if (not 1 <= start <= end or not valid_score or f[0] in ('', '.') or not f[1]
                    or f[2] in ('', '.') or f[6] not in ('+', '-') or f[7] not in ('.', '0', '1', '2')):
                raise ValueError(context + ': invalid GFF3 record')
            attrs = {}
            for item in f[8].split(';'):
                if not item: continue
                key, sep, value = item.partition('=')
                if not sep or not key or key in attrs:
                    raise ValueError(context + ': malformed or duplicate attribute')
                attrs[key] = value
            if f[2] != 'gene' and not attrs.get('transcript_id', '').strip():
                raise ValueError(context + ': missing transcript_id')
            yield Record(tuple(f[:8]), attrs)


def key(record):
    return record.fields[0], record.attrs['transcript_id']


def feature(record):
    f = record.fields
    return preservation.Feature(f[0], f[1], f[2], int(f[3]), int(f[4]), f[5], f[6], f[7],
        tuple(sorted((k, v) for k, v in record.attrs.items() if k not in IGNORED)))


def text(record):
    return '\t'.join(record.fields) + '\t' + ';'.join(k + '=' + v for k, v in record.attrs.items()) + '\n'


def changed(record, start=None, end=None, attrs=None):
    f = list(record.fields)
    if start is not None: f[3] = str(start)
    if end is not None: f[4] = str(end)
    return Record(tuple(f), dict(record.attrs) if attrs is None else attrs)


def write_gff(path, records):
    with path.open('w', encoding='utf-8') as handle:
        handle.write('##gff-version 3\n')
        for record in sorted(records, key=lambda r: (r.fields[0], int(r.fields[3]), int(r.fields[4]),
                0 if r.fields[2] == 'gene' else 1 if r.fields[2] == 'transcript' else 2, r.fields[2], compact(r.attrs))):
            handle.write(text(record))


def write_csv(path, fields, rows, delimiter=','):
    with path.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter=delimiter, extrasaction='ignore')
        writer.writeheader(); writer.writerows(rows)


def grouped(records):
    models = defaultdict(list)
    for record in records:
        if record.fields[2] != 'gene': models[key(record)].append(record)
    return models


def validate_raw(models, directory):
    """Scope to represented models; unrelated filtered raw models are not repaired."""
    result = {}
    chromosomes = defaultdict(list)
    for (chromosome, identifier), records in models.items(): chromosomes[chromosome].extend(records)
    for n, records in enumerate(chromosomes.values()):
        path = directory / ('raw_' + str(n) + '.gff3')
        write_gff(path, records)
        for identifier, model in preservation.read_annotation(path).items():
            if preservation.extension_consistency(model, model):
                raise ValueError('Raw transcript/exon boundaries disagree: ' + identifier)
            result[(model['transcript'].chromosome, identifier)] = model
    return result


def rebuild(old_records, raw_records, model, entry, candidates, stage):
    chromosome, identifier = key(raw_records[0])
    old_tx = [r for r in old_records if r.fields[2] == 'transcript']
    if len(old_tx) > 1:
        raise ValueError('Duplicate transcript rows make the selected boundaries ambiguous')
    raw_tx = next(r for r in raw_records if r.fields[2] == 'transcript')
    strand = raw_tx.fields[6]
    if old_tx:
        if old_tx[0].fields[6] != strand: raise ValueError('Edited transcript strand disagrees with raw')
        start, end = map(int, old_tx[0].fields[3:5])
    else:
        location_chrom, location = entry['transcript_location'].rsplit(':', 1)
        start, end = map(int, location.split('-'))
        if location_chrom != chromosome: raise ValueError('Stats chromosome differs from raw')
        # Without a transcript row, require two agreeing descriptions of its bounds.
        left, right = int(raw_tx.fields[3]) - start, end - int(raw_tx.fields[4])
        lengths = (left, right) if strand == '+' else (right, left)
        if list(lengths) != [entry[s + '_extension'] for s in consistency.SIDES]:
            raise ValueError('Missing transcript row and inconsistent stats: cannot infer selected boundaries')
    raw_start, raw_end = int(raw_tx.fields[3]), int(raw_tx.fields[4])
    if start > raw_start or end < raw_end:
        raise ValueError('Selected transcript is shortened; repair will not choose new outer boundaries')
    if entry['gene_id'] != model['gene']:
        raise ValueError('Stats gene membership disagrees with raw')
    classes = {r.attrs['transcript_class'] for r in old_records if r.attrs.get('transcript_class')}
    if len(classes) > 1: raise ValueError('Conflicting transcript_class values')
    tags = {t for r in old_tx for t in r.attrs.get('leap_tag', '').split(',')}
    blacklisted = 'LEAP_blacklisted' in tags or entry['five_prime_support'] == 'blacklisted'
    row = dict(entry, transcript_location='{}:{}-{}'.format(chromosome, start, end))
    expected_tags = set()
    baseline = dict(chromosome=chromosome, strand=strand, gene_id=model['gene'],
        exons=[dict(start=f.start, end=f.end) for f in model['exons']])
    for side in consistency.SIDES:
        left = (side == 'five_prime') == (strand == '+')
        length = raw_start - start if left else end - raw_end
        row[side + '_extension'] = length
        if side == 'five_prime' and blacklisted:
            if length: raise ValueError('A blacklisted 5-prime end cannot remain extended')
            row[side + '_support'] = row[side + '_read_count'] = 'blacklisted'
            expected_tags.add('LEAP_blacklisted')
        elif not length:
            row[side + '_support'], row[side + '_read_count'] = 'none', '0'
        else:
            if candidates is None: raise ValueError('Extended ends require --leap-output-dir to verify candidate provenance')
            candidate = candidates.candidate(baseline, side)
            endpoint = start if left else end
            if (candidate is None or model['gene'] not in candidate['genes'] or candidate['read_count'] is None
                    or (endpoint > candidate['endpoint'] if not left else endpoint < candidate['endpoint'])
                    or (stage == 'leap' and endpoint != candidate['endpoint'])
                    or (candidate['support'] == 'reads' and candidate['read_count'] <= 0)):
                raise ValueError('Selected endpoint is incompatible with the last candidate for ' + side)
            support, count = candidate['support'], str(candidate['read_count'])
            # Preserve equivalent repeated support entries from the existing writer.
            if not (all(s == support for s in row[side + '_support'].split(',')) and
                    all(n == count for n in row[side + '_read_count'].split(','))):
                row[side + '_support'], row[side + '_read_count'] = support, count
            expected_tags.add('LEAP_5' if side == 'five_prime' else 'LEAP_3')
    output = []
    for record in raw_records:
        f = feature(record)
        a, b = f.start, f.end
        if f.kind == 'transcript': a, b = start, end
        if f == model['exons'][0]: a = start
        if f == model['exons'][-1]: b = end
        if f.kind == preservation.terminal_utr(strand, 'left') and f.start == raw_start: a = start
        if f.kind == preservation.terminal_utr(strand, 'right') and f.end == raw_end: b = end
        attrs = {k: v for k, v in record.attrs.items() if k not in IGNORED}
        if f.kind == 'transcript':
            expected = expected_tags or {'LEAP_none'}
            previous = old_tx[0].attrs.get('leap_tag', '') if old_tx else ''
            attrs['leap_tag'] = previous if set(previous.split(',')) == expected else ','.join(
                tag for tag in ('LEAP_5', 'LEAP_3', 'LEAP_blacklisted', 'LEAP_none') if tag in expected)
            if classes: attrs['transcript_class'] = next(iter(classes))
        output.append(changed(record, a, b, attrs))
    return output, row


def gene_records(repaired, raw_genes):
    genes = {}
    for record in repaired:
        if record.fields[2] != 'transcript': continue
        parent = record.attrs.get('Parent', '')
        if not parent or ',' in parent: raise ValueError('Raw transcript must have exactly one gene Parent')
        k = record.fields[0], parent
        if k not in genes:
            template = raw_genes.get(k)
            if template is None:
                attrs = {a: v for a, v in record.attrs.items() if a.startswith('gene_')}
                attrs['ID'] = parent
                f = list(record.fields); f[2] = 'gene'; f[7] = '.'
                template = Record(tuple(f), attrs)
            genes[k] = changed(template, int(record.fields[3]), int(record.fields[4]))
        else:
            gene = genes[k]
            if gene.fields[6] != record.fields[6]: raise ValueError('Gene has conflicting transcript strands')
            genes[k] = changed(gene, min(int(gene.fields[3]), int(record.fields[3])),
                               max(int(gene.fields[4]), int(record.fields[4])))
    return list(genes.values())


def audit_row(k, before, after, old_stats, new_stats):
    def full(r):
        return feature(r)._replace(attributes=tuple(sorted(r.attrs.items())))
    a, b = Counter(full(r) for r in before), Counter(full(r) for r in after)
    removed = [preservation.record_details(f, n) for f, n in sorted((a - b).items())]
    added = [preservation.record_details(f, n) for f, n in sorted((b - a).items())]
    old_tx = next((r for r in before if r.fields[2] == 'transcript'), None)
    new_tx = next(r for r in after if r.fields[2] == 'transcript')
    generated = dict(before=old_tx.attrs if old_tx else {}, after=new_tx.attrs)
    generated = {label: {a: v for a, v in attrs.items() if a in IGNORED} for label, attrs in generated.items()}
    stats_changes = {field: dict(before=old_stats[field], after=new_stats[field])
                     for field in consistency.STATS_FIELDS if old_stats[field] != new_stats[field]}
    different = bool(removed or added or stats_changes or generated['before'] != generated['after'])
    return dict(chromosome=k[0], transcript_id=k[1], status='repaired' if different else 'unchanged',
        original_start=old_tx.fields[3] if old_tx else '', original_end=old_tx.fields[4] if old_tx else '',
        repaired_start=new_tx.fields[3], repaired_end=new_tx.fields[4], removed_records=sum((a - b).values()),
        added_records=sum((b - a).values()), removed_features=compact(removed), added_features=compact(added),
        stats_changes=compact(stats_changes), generated_attributes=compact(generated))


def run(args):
    output = Path(args.output_dir).absolute()
    if output.exists() or output.is_symlink(): raise ValueError('Output must be a new directory: ' + str(output))
    output.parent.mkdir(parents=True, exist_ok=True)
    rows, validations, status, message, repaired, corrected, genes = [], [], 'INPUT_ERROR', '', [], {}, []
    with tempfile.TemporaryDirectory(prefix='.' + output.name + '-', dir=output.parent) as temporary:
        staging = Path(temporary) / 'report'; staging.mkdir()
        try:
            raw_all = list(read_records(args.current_annotation))
            edited_all = list(read_records(args.extended_annotation))
            edited = grouped(edited_all)
            if not edited: raise ValueError('No represented transcript models in edited GFF')
            stats = consistency.read_stats(args.extension_stats)
            if set(stats) != set(edited):
                raise ValueError('GFF/stats membership differs; repair does not infer retention or add wholly missing models')
            raw_groups = grouped(raw_all)
            if not set(edited) <= set(raw_groups): raise ValueError('Edited transcript ID/chromosome absent from raw; exact versions are required')
            raw_groups = {k: raw_groups[k] for k in edited}
            models = validate_raw(raw_groups, Path(temporary))
            raw_genes = {}
            for record in raw_all:
                if record.fields[2] != 'gene': continue
                k = record.fields[0], record.attrs.get('ID', '')
                if k in raw_genes: raise ValueError('Duplicate raw gene ID: ' + str(k))
                raw_genes[k] = record
            context = consistency.CandidateContext(args.leap_output_dir,
                {k: dict(chromosome=k[0], strand=m['transcript'].strand) for k, m in models.items()}) if args.leap_output_dir else None
            for k, before in sorted(edited.items()):
                try:
                    after, row = rebuild(before, raw_groups[k], models[k], stats[k], context, args.stage)
                    rows.append(audit_row(k, before, after, stats[k], row))
                    repaired.extend(after); corrected[k] = row
                except ValueError as error:
                    rows.append(dict(chromosome=k[0], transcript_id=k[1], status='unrepairable', message=str(error)))
            if any(r['status'] == 'unrepairable' for r in rows):
                status, message = 'UNREPAIRABLE', 'At least one model is ambiguous or has unverified selected boundaries'
            else:
                if args.stage == 'leap': genes = gene_records(repaired, raw_genes)
                write_gff(staging / 'repaired.gff3', repaired + genes)
                write_csv(staging / 'repaired_stats.tsv', consistency.STATS_FIELDS,
                          [corrected[k] for k in sorted(corrected)], '\t')
                repaired_groups = grouped(repaired)
                checked_models = validate_raw(repaired_groups, Path(temporary))
                for k, records in repaired_groups.items():
                    checked = checked_models[k]
                    comparison = next(preservation.compare_transcripts({k[1]: models[k]}, {k[1]: checked}))
                    if comparison['failed']: raise ValueError('Post-repair preservation validation failed: ' + k[1])
                validations.append(dict(check='transcript_preservation', status='PASS', details='All retained models preserve raw protected features'))
                raw_stats = {(k[0], k[1]): dict(chromosome=k[0], strand=m['transcript'].strand,
                    gene_id=m['gene'], start=m['transcript'].start, end=m['transcript'].end) for k, m in models.items()}
                final_stats = consistency.read_annotation(staging / 'repaired.gff3', include_exons=True)
                failures = [r for r in consistency.compare(raw_stats, final_stats, corrected) if r['failed']]
                if failures: raise ValueError('Post-repair stats/GFF validation failed: ' + compact(failures))
                validations.append(dict(check='stats_gff', status='PASS', details='Lengths, support, counts and tags agree; candidate provenance verified during repair'))
                issues, _ = structural.check_annotation(structural.parse_annotation(staging / 'repaired.gff3'), allow_missing_genes=args.stage == 'processed')
                if issues: raise ValueError('Post-repair structural validation failed: ' + compact(issues))
                validations.append(dict(check='structural_integrity', status='PASS', details='IDs, parents and containment validated'))
                status = 'PASS'
        except (OSError, EOFError, UnicodeError, ValueError, KeyError, csv.Error) as error:
            message = str(error)
        if status != 'PASS':
            for name in ('repaired.gff3', 'repaired_stats.tsv'):
                (staging / name).unlink(missing_ok=True)
        gene_changes = dict(before=[text(r).rstrip() for r in edited_all if r.fields[2] == 'gene'],
                            after=[text(r).rstrip() for r in genes]) if 'edited_all' in locals() else {}
        write_csv(staging / 'changes.csv', FIELDS, rows)
        write_csv(staging / 'validation.csv', ['check', 'status', 'details'], validations)
        (staging / 'gene_changes.json').write_text(json.dumps(gene_changes, indent=2) + '\n')
        summary = dict(status=status, message=message, stage=args.stage, counts=dict(Counter(r['status'] for r in rows)),
            gene_records_changed=(Counter(feature(r) for r in edited_all if r.fields[2] == 'gene') != Counter(feature(r) for r in genes)) if 'edited_all' in locals() else False,
            inputs=dict(raw=str(Path(args.current_annotation).resolve()), edited=str(Path(args.extended_annotation).resolve()),
                        stats=str(Path(args.extension_stats).resolve()), leap_output_dir=args.leap_output_dir),
            validation='Preservation, stats/GFF and structural checks only; rerun the full datacheck suite before release.')
        (staging / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
        staging.rename(output)
    print('{}: {}; {}'.format(status, summary['counts'], message or 'repaired copies and audit saved'))
    return {'PASS': 0, 'UNREPAIRABLE': 1, 'INPUT_ERROR': 2}[status]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('current_annotation'); parser.add_argument('extended_annotation'); parser.add_argument('output_dir')
    parser.add_argument('--extension-stats', required=True)
    parser.add_argument('--leap-output-dir', help='Original run context; required to retain any extended end')
    parser.add_argument('--stage', choices=('leap', 'processed'), default='leap')
    args = parser.parse_args()
    try:
        return run(args)
    except (OSError, ValueError, csv.Error) as error:
        print('ERROR: ' + str(error), file=sys.stderr); return 2


if __name__ == '__main__':
    sys.exit(main())
