#!/usr/bin/env python3
"""Build offline LEAP impact statistics from declared GFF/stats partitions."""

import argparse
from bisect import bisect_left
from collections import Counter, defaultdict
import csv
import json
from pathlib import Path
import re
import statistics
import sys
import tempfile

sys.dont_write_bytecode = True
import check_stats_gff_consistency as consistency
import process_internal_extensions as processor

STAGES = ('leap', 'external', 'internal')
STRANDS = {'positive': '+', 'negative': '-'}
SIDES = consistency.SIDES
BINS = [(1, 100), (101, 500), (501, 1000), (1001, 5000), (5001, 10000), (10001, None)]
COVERAGE_FIELDS = ['unique_genomic_extension_bases', 'newly_covered_bases', 'newly_exonic_bases',
                   'new_same_strand_feature_bases', 'unique_genomic_read_supported_bases', 'new_read_supported_bases']
COVERAGE_DEFINITIONS = dict(
    newly_covered='All extension intervals merged by chromosome across transcripts and strands, minus every raw GFF feature interval on either strand. Gene, transcript and unfamiliar feature spans are included. Each genomic base counts once per result view.',
    newly_exonic='Extension union minus all raw exons on either strand; includes new exonic coverage inside existing gene or transcript spans. Each genomic base counts once per result view.',
    new_same_strand='Extensions minus raw features on their own strand, with unstranded features protecting both; resulting positions are merged across strands and counted once.',
    new_read_supported='Read-supported extension union across transcripts and strands, minus every raw feature on either strand. Each genomic base counts once per result view.')
TRANSCRIPT_FIELDS = ['chromosome', 'strand', 'gene_id', 'transcript_id', 'jbrowse_label',
    'raw_start', 'raw_end', 'leap_start', 'leap_end', 'external_start', 'external_end',
    'classification', 'dropped_no_cds'] + [side + '_' + field for side in SIDES
    for field in ('extension', 'support', 'read_count', 'external_extension', 'internal_extension')]


def compact(value):
    return json.dumps(value, separators=(',', ':'), ensure_ascii=False)


def render_dashboard(summary):
    payload = compact(summary).replace('&', '\\u0026').replace('<', '\\u003c').replace('>', '\\u003e')
    template = Path(__file__).with_name('leap_impact_dashboard.html').read_text(encoding='utf-8')
    return template.replace('__LEAP_IMPACT_DATA__', payload)


def written_summary(summary):
    text = [summary['label'], 'Statistics: {}; partitions: {}/{}; datachecks: {}'.format(summary['status'],
        summary['completed_partitions'], summary['expected_partitions'], summary['datachecks']['status']), '']
    for stage, values in summary['totals'].items():
        text += [stage.upper(), '{} / {} models extended; {} genes extended'.format(
            values.get('extended_transcripts', 0), values.get('transcripts', 0), values.get('extended_genes', 0)),
            '{} unique genomic extension bp, merged across transcripts and strands'.format(values.get('unique_genomic_extension_bases', 0)),
            '{} added bp summed across transcripts (not a unique base count)'.format(values.get('added_bases', 0)),
            '{} new read-supported genomic bp outside every raw feature on either strand'.format(
                values.get('new_read_supported_bases', 0)),
            '{} newly covered genomic bp outside every raw feature on either strand'.format(values.get('newly_covered_bases', 0)),
            '{} new exonic bp outside raw exons on either strand'.format(values.get('newly_exonic_bases', 0)), '']
    text += ['Internal and external models/unions can overlap; do not add their counts or unique base totals.',
             'Each end has one support source: reads or current annotation. Datachecks describe original LEAP output.',
             'See index.html, CSV files and summary.json for complete definitions, examples and diagnostics.']
    return '\n'.join(text) + '\n'


def merged(intervals):
    result = []
    for start, end in sorted(intervals):
        if result and start <= result[-1][1] + 1:
            result[-1] = (result[-1][0], max(end, result[-1][1]))
        else:
            result.append((start, end))
    return result


def bases(intervals):
    return sum(b - a + 1 for a, b in intervals)


def coverage_accumulator():
    return {name: [] for name in ('all_features', 'all_exons', 'positive_features', 'negative_features')}


def add_raw_coverage(coverage, fields):
    interval = int(fields[3]), int(fields[4])
    coverage['all_features'].append(interval)
    if fields[2] == 'exon':
        coverage['all_exons'].append(interval)
    for label, sign in STRANDS.items():
        if fields[6] in (sign, '.', '?'):
            coverage[label + '_features'].append(interval)


def extension_interval(row, side, length):
    left = (side == 'five_prime') == (row['strand'] == 'positive')
    return ((row['raw_start'] - length, row['raw_start'] - 1) if left else
            (row['raw_end'] + 1, row['raw_end'] + length))


def genomic_coverage(rows, raw_coverage):
    """Union by chromosome, including across strands; never sum isoform coverage."""
    additions, reads = defaultdict(list), defaultdict(list)
    for row in rows:
        for stage in STAGES:
            if stage != 'leap' and row['dropped_no_cds']:
                continue
            selected = stage_row(row, stage)
            for side in SIDES:
                length = selected[side + '_extension']
                if length:
                    key = row['chromosome'], stage, row['strand']
                    interval = extension_interval(row, side, length)
                    additions[key].append(interval)
                    if support_category(selected[side + '_support']) == 'reads':
                        reads[key].append(interval)
    results, intervals = [], []
    for chromosome, protected in sorted(raw_coverage.items()):
        for stage in STAGES:
            by_strand = {s: merged(additions[(chromosome, stage, s)]) for s in STRANDS}
            same_strand = {s: outside(by_strand[s], protected[s + '_features']) for s in STRANDS}
            for view in ('positive', 'negative', 'both'):
                selected = list(STRANDS) if view == 'both' else [view]
                union = merged(i for s in selected for i in by_strand[s])
                read_union = merged(i for s in selected for i in reads[(chromosome, stage, s)])
                new = outside(union, protected['all_features'])
                results.append(dict(stage=stage, chromosome=chromosome, strand=view,
                    unique_genomic_extension_bases=bases(union), newly_covered_bases=bases(new),
                    unique_genomic_read_supported_bases=bases(read_union),
                    new_read_supported_bases=bases(outside(read_union, protected['all_features'])),
                    newly_exonic_bases=bases(outside(union, protected['all_exons'])),
                    new_same_strand_feature_bases=bases(merged(i for s in selected for i in same_strand[s]))))
                if view == 'both':
                    intervals.extend(dict(stage=stage, chromosome=chromosome, start=a, end=b, bases=b-a+1)
                                     for a, b in new)
    return results, intervals


def add_coverage_totals(summary, coverage):
    summary['schema_version'] = 2
    summary['genomic_coverage'] = coverage
    summary['definitions'].update(COVERAGE_DEFINITIONS)
    for stage in STAGES:
        summary['totals'][stage].update({field: sum(r[field] for r in coverage
            if r['stage'] == stage and r['strand'] == 'both') for field in COVERAGE_FIELDS})


def outside(intervals, protected, protected_ends=None):
    """Both collections are inclusive, sorted, disjoint interval unions."""
    if not intervals:
        return []
    if protected_ends is None:
        protected_ends = [b for a, b in protected]
    result, index = [], bisect_left(protected_ends, intervals[0][0])
    for start, end in intervals:
        while index < len(protected) and protected[index][1] < start:
            index += 1
        cursor, other = start, index
        while other < len(protected) and protected[other][0] <= end:
            a, b = protected[other]
            if a > cursor:
                result.append((cursor, min(a - 1, end)))
            cursor = max(cursor, b + 1)
            if cursor > end:
                break
            other += 1
        if cursor <= end:
            result.append((cursor, end))
    return result


def support_category(value):
    labels = set(value.split(','))
    if len(labels) != 1:
        raise ValueError('Each transcript end must have one support source; conflicting labels: ' + value)
    return next(iter(labels))


def table(path, fields, rows):
    with path.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction='ignore')
        writer.writeheader()
        for row in rows:
            # Keep exact IDs; neutralise spreadsheet formulas in free-text fields.
            writer.writerow({k: "'" + v if isinstance(v, str) and v.startswith(('=', '+', '-', '@')) else v
                             for k, v in row.items()})


def annotations(path, chromosome, strand, raw=False, coverage=None):
    transcripts, exons, kinds = {}, defaultdict(list), defaultdict(Counter)
    for f, attrs, context in processor.records(path):
        if f[0] != chromosome or (not raw and f[6] != strand):
            raise ValueError(context + ': annotation is outside selected chromosome/strand')
        if coverage is not None:
            add_raw_coverage(coverage, f)
        if f[6] != strand or f[2] == 'gene':
            continue
        identifier = attrs.get('transcript_id')
        if raw and not identifier and f[2] not in ('transcript', 'exon', 'CDS', 'start_codon',
                                                  'stop_codon', 'five_prime_UTR', 'three_prime_UTR'):
            # Other gene-level or reference features still protect coverage.
            continue
        if not identifier or identifier == '.' or not attrs.get('gene_id'):
            raise ValueError(context + ': missing transcript_id or gene_id')
        k = (chromosome, identifier)
        kinds[k][f[2]] += 1
        if f[2] == 'transcript':
            if k in transcripts:
                raise ValueError(context + ': duplicate transcript row')
            transcripts[k] = dict(chromosome=chromosome, strand=strand, gene_id=attrs['gene_id'],
                start=int(f[3]), end=int(f[4]), leap_tag=attrs.get('leap_tag'),
                transcript_type=attrs.get('transcript_type'), context=context)
        elif f[2] == 'exon':
            exons[k].append(dict(start=int(f[3]), end=int(f[4]), strand=strand, gene_id=attrs['gene_id']))
    for k, model in transcripts.items():
        model['exons'] = sorted(exons[k], key=lambda r: (r['start'], r['end']))
    if not raw and set(kinds) != set(transcripts):
        raise ValueError(str(path) + ': child records without transcript rows')
    return transcripts, exons, kinds


def stats(path):
    with consistency.open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        header = reader.fieldnames or []
        if len(header) != len(set(header)) or not set(consistency.STATS_FIELDS) <= set(header):
            raise ValueError(str(path) + ': missing or duplicate stats columns')
        has_rows = next(reader, None) is not None
    return consistency.read_stats(path) if has_rows else {}


def analyse_partition(job):
    chromosome, label = job['chromosome'], job['strand']
    strand = STRANDS[label]
    paths = {}
    for field in ('raw', 'gff', 'stats'):
        if len(job[field]) != 1:
            raise ValueError('Expected exactly one {} input; received {}'.format(field, len(job[field])))
        paths[field] = Path(job[field][0])
    raw_coverage = coverage_accumulator()
    raw, raw_exons, kinds = annotations(paths['raw'], chromosome, strand, raw=True, coverage=raw_coverage)
    raw_coverage = {k: merged(v) for k, v in raw_coverage.items()}
    final, _, _ = annotations(paths['gff'], chromosome, strand)
    entries = stats(paths['stats'])
    if set(entries) != set(final):
        raise ValueError('Final annotation/stats membership differs')
    findings = [r for r in consistency.compare(raw, final, entries) if r['failed']]
    if findings:
        raise ValueError('GFF/stats inconsistency: ' + compact([
            dict(transcript_id=r['transcript_id'], details=r['details']) for r in findings[:10]]))
    protected = merged((r['start'], r['end']) for records in raw_exons.values() for r in records)
    protected_ends = [b for a, b in protected]
    eligible = sum(r['transcript_type'] == 'protein_coding' and kinds[k]['exon'] > 1
                   and kinds[k]['five_prime_UTR'] > 0 and kinds[k]['three_prime_UTR'] > 0
                   for k, r in raw.items())
    models = {}
    if final:
        models, _ = processor.scan_extended(paths['gff'])
        baseline, _ = processor.scan_raw(paths['raw'], models)
        processor.prepare_models(models, baseline, processor.load_stats(paths['stats']))
    rows = []
    for k, model in final.items():
        entry, original, processed = entries[k], raw[k], models.get(k)
        ext = processed['output_stats'] if processed else None
        coordinates = ext['transcript_location'].rsplit(':', 1)[1].split('-') if ext else ('', '')
        row = dict(chromosome=chromosome, strand=label, gene_id=model['gene_id'], transcript_id=k[1],
            jbrowse_label='{}:{}-{}'.format(chromosome, model['start'], model['end']),
            raw_start=original['start'], raw_end=original['end'], leap_start=model['start'], leap_end=model['end'],
            external_start=int(coordinates[0]) if ext else '', external_end=int(coordinates[1]) if ext else '',
            classification=processed['class'] if processed else 'dropped_no_cds', dropped_no_cds=not bool(processed))
        for side in SIDES:
            for field in ('extension', 'support', 'read_count'):
                row[side + '_' + field] = entry[side + '_' + field]
            length = int(entry[side + '_extension'])
            retained = int(ext[side + '_extension']) if ext else 0
            row[side + '_external_extension'] = retained
            row[side + '_internal_extension'] = length - retained if processed else 0
        rows.append(row)
    groups, genes = [], []
    baseline = dict(raw_transcripts=len(raw), eligible_raw_transcripts=eligible,
        raw_genes=len({r['gene_id'] for r in raw.values()}), raw_exon_union_bases=bases(protected))
    for stage in STAGES:
        selected = [stage_row(r, stage) for r in rows
                    if stage == 'leap' or (not r['dropped_no_cds'] and
                        (stage == 'external' or any(r[s + '_internal_extension'] for s in SIDES)))]
        group = measure(selected, protected, protected_ends)
        group.update(chromosome=chromosome, strand=label, stage=stage, **baseline)
        group['top_transcripts'] = sorted(selected,
            key=lambda r: (-sum(r[s + '_extension'] for s in SIDES), r['transcript_id']))[:20]
        groups.append(group)
        by_gene = defaultdict(list)
        for row in selected:
            by_gene[row['gene_id']].append(row)
        for identifier, members in sorted(by_gene.items()):
            measured = measure(members, protected, protected_ends)
            genes.append(dict(chromosome=chromosome, strand=label, stage=stage, gene_id=identifier,
                **{key: value for key, value in measured.items() if isinstance(value, (int, float))}))
    return dict(groups=groups, genes=genes, transcripts=rows, baseline=baseline,
        raw_coverage=raw_coverage, sources={k: str(v) for k, v in paths.items()})


def stage_row(row, stage):
    result = dict(row)
    for side in SIDES:
        length = row[side + ('_external_extension' if stage == 'external' else
                             '_internal_extension' if stage == 'internal' else '_extension')]
        result[side + '_extension'] = length
        if not length and not (side == 'five_prime' and row[side + '_support'] == 'blacklisted'):
            result[side + '_support'], result[side + '_read_count'] = 'none', '0'
    if stage == 'external':
        result['jbrowse_label'] = '{}:{}-{}'.format(row['chromosome'], row['external_start'], row['external_end'])
    return result


def measure(rows, raw_exons, raw_exon_ends=None):
    lengths = {side: [] for side in SIDES}
    supports = {side: Counter() for side in SIDES}
    intervals, reads, genes, read_counts = [], [], set(), []
    histogram = {side: [0] * len(BINS) for side in SIDES}
    metrics = Counter(transcripts=len(rows), extended_transcripts=0, unchanged_transcripts=0,
        extended_genes=0, both_ends_extended=0, added_bases=0, reads_added_bases=0,
        annotation_added_bases=0, blacklisted_transcripts=0, dropped_no_cds=0,
        five_prime_extended=0, three_prime_extended=0, five_prime_added_bases=0, three_prime_added_bases=0)
    for row in rows:
        extended = any(row[s + '_extension'] > 0 for s in SIDES)
        metrics['extended_transcripts' if extended else 'unchanged_transcripts'] += 1
        if extended:
            genes.add(row['gene_id'])
        metrics['both_ends_extended'] += all(row[s + '_extension'] > 0 for s in SIDES)
        metrics['blacklisted_transcripts'] += row['five_prime_support'] == 'blacklisted'
        metrics['dropped_no_cds'] += row['dropped_no_cds']
        for side in SIDES:
            length = row[side + '_extension']
            category = support_category(row[side + '_support'])
            supports[side][category] += 1
            if not length:
                continue
            metrics[side + '_extended'] += 1
            metrics[side + '_added_bases'] += length
            metrics['added_bases'] += length
            metrics[category + '_added_bases' if category != 'current_annotation' else 'annotation_added_bases'] += length
            lengths[side].append(length)
            for i, (a, b) in enumerate(BINS):
                if a <= length and (b is None or length <= b):
                    histogram[side][i] += 1
                    break
            interval = extension_interval(row, side, length)
            intervals.append(interval)
            if category == 'reads':
                reads.append(interval)
                read_counts.extend(int(n) for s, n in zip(row[side + '_support'].split(','),
                    row[side + '_read_count'].split(',')) if s == 'reads')
    metrics['extended_genes'] = len(genes)
    metrics['genes'] = len({r['gene_id'] for r in rows})
    metrics['unique_added_bases'] = bases(merged(intervals))
    metrics['unique_read_supported_bases'] = bases(merged(reads))
    metrics['novel_read_supported_bases'] = bases(outside(merged(reads), raw_exons, raw_exon_ends))
    metrics['recorded_read_count_min'] = min(read_counts, default=0)
    metrics['recorded_read_count_max'] = max(read_counts, default=0)
    result = dict(metrics, support={s: dict(v) for s, v in supports.items()}, histogram=histogram)
    for side, values in lengths.items():
        result[side + '_lengths'] = describe(values)
    return result


def describe(values):
    values = sorted(values)
    def percentile(p):
        if not values:
            return 0
        position = (len(values) - 1) * p
        a = int(position)
        return values[a] + (values[min(a + 1, len(values) - 1)] - values[a]) * (position - a)
    return dict(count=len(values), mean=statistics.mean(values) if values else 0,
                median=percentile(.5), p90=percentile(.9), p95=percentile(.95), maximum=max(values, default=0))


def one(folder, name, suffixes):
    return [str(p) for suffix in suffixes if (p := folder / (name + suffix)).is_file()]


def discover(root, chromosomes, strand):
    if not chromosomes or len(set(chromosomes)) != len(chromosomes):
        raise ValueError('Require an explicit unique --chromosomes list')
    jobs = []
    for chromosome in chromosomes:
        for label in STRANDS if strand == 'both' else [strand]:
            jobs.append(dict(chromosome=chromosome, strand=label,
                raw=one(root / 'raw_annotation', 'current_annotation_' + chromosome, ('.gff3', '.gff3.gz')),
                gff=one(root / 'extended_transcripts', chromosome + '_' + label + '_extended', ('.gff', '.gff.gz', '.gff3', '.gff3.gz')),
                stats=one(root / 'extended_transcripts_stats', chromosome + '_' + label + '_stats', ('.tsv', '.tsv.gz'))))
    return dict(partitions=jobs, issues=[])


def build(manifest, output, check_summary=None, label='LEAP impact'):
    output = Path(output).absolute()
    if output.exists() or output.is_symlink():
        raise ValueError('Output must be a new directory: ' + str(output))
    jobs = manifest['partitions']
    scopes = [(j['chromosome'], j['strand']) for j in jobs]
    if not scopes or len(scopes) != len(set(scopes)):
        raise ValueError('Expected nonempty unique chromosome/strand scopes')
    for chromosome, strand in scopes:
        if not re.fullmatch(r'[A-Za-z0-9_.-]+', chromosome) or chromosome in ('.', '..') or strand not in STRANDS:
            raise ValueError('Invalid chromosome/strand scope')
    errors = [dict(chromosome='', strand='', message=str(issue)) for issue in manifest.get('issues', [])]
    groups, genes, rows, provenance, raw_coverage = [], [], [], [], {}
    for job in jobs:
        try:
            result = analyse_partition(job)
            chromosome = job['chromosome']
            if chromosome in raw_coverage and raw_coverage[chromosome] != result['raw_coverage']:
                raise ValueError('Different raw feature coverage supplied for strands of chromosome ' + chromosome)
            raw_coverage[chromosome] = result['raw_coverage']
            groups += result['groups']; genes += result['genes']; rows += result['transcripts']
            provenance.append(dict(chromosome=job['chromosome'], strand=job['strand'], **result['sources']))
            print('Analysed {} {}: {} output transcripts'.format(job['chromosome'], job['strand'],
                len(result['transcripts'])), flush=True)
        except (OSError, EOFError, UnicodeError, ValueError, KeyError, csv.Error) as error:
            errors.append(dict(chromosome=job['chromosome'], strand=job['strand'], message=str(error)))
            print('INPUT_ERROR {} {}: {}'.format(job['chromosome'], job['strand'], error), file=sys.stderr, flush=True)
    checks = dict(status='NOT_RUN', counts={}, results=[])
    if check_summary:
        path = Path(check_summary)
        if path.is_dir():
            path /= 'summary.json'
        try:
            checks = json.loads(path.read_text())
            if (checks.get('status') not in ('PASS', 'REVIEW', 'FAIL', 'INPUT_ERROR')
                    or not isinstance(checks['results'], list) or not isinstance(checks['counts'], dict)
                    or any(not isinstance(r, dict) or r.get('status') not in
                        ('PASS', 'REVIEW', 'FAIL', 'INPUT_ERROR', 'NOT_APPLICABLE') for r in checks['results'])):
                raise ValueError('Invalid datacheck summary schema')
        except (OSError, ValueError, KeyError) as error:
            errors.append(dict(chromosome='', strand='', message='Datacheck summary: ' + str(error)))
            checks = dict(status='INPUT_ERROR', counts={}, results=[])
    coverage, new_intervals = genomic_coverage(rows, raw_coverage)
    summary = dict(schema_version=2, label=label, status='INPUT_ERROR' if errors else 'COMPLETE',
        exit_code=2 if errors else 0, expected_partitions=len(scopes), completed_partitions=len(provenance),
        groups=groups, genomic_coverage=coverage, datachecks=checks, diagnostics=errors, sources=provenance, expected_sources=jobs,
        definitions=dict(coordinates='Inclusive genomic bases, exact chromosome names; strand-aware transcript ends.',
            unique_bases='Legacy unique_added_bases, unique_read_supported_bases and novel_read_supported_bases columns merge within each strand. New genomic_coverage metrics also merge across strands, counting each chromosome position once.',
            read_supported='The selected end is read-supported. Each end has one support source; annotation and reads may support opposite ends of one transcript.',
            novel='Read-supported extension union outside every raw exon on the same chromosome and strand, including partial and noncoding exons.',
            internal='Only the ends restored by the existing internal-extension processor; a model may also retain an external end.',
            counts='Transcript totals count exact versioned IDs separately on each chromosome/strand. Gene totals use chromosome/strand/gene_id.',
            exclusions='External processing excludes models without CDS; those are reported separately from internal additions.',
            limitations='Added sequence is genomic span growth, not measured read coverage, new molecules, new transcripts, coding growth or assembly validation. Quantiles describe positive extension ends only.'))
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.' + output.name + '-', dir=output.parent) as temp:
        staging = Path(temp) / 'impact'; staging.mkdir()
        table(staging / 'transcripts.csv', TRANSCRIPT_FIELDS, rows)
        metric_fields = sorted({k for g in groups for k, v in g.items() if isinstance(v, (int, float))})
        table(staging / 'partition_statistics.csv', ['stage', 'chromosome', 'strand'] + metric_fields, groups)
        table(staging / 'genes.csv', ['stage', 'chromosome', 'strand', 'gene_id'] + sorted(
            {k for g in genes for k, v in g.items() if isinstance(v, (int, float))}), genes)
        support_rows, lengths = [], []
        for group in groups:
            identity = {k: group[k] for k in ('stage', 'chromosome', 'strand')}
            for side in SIDES:
                support_rows.extend(dict(identity, end=side, support=k, transcript_ends=v)
                                    for k, v in group['support'][side].items())
        for stage in STAGES:
            selected = [stage_row(r, stage) for r in rows if stage == 'leap' or
                (not r['dropped_no_cds'] and (stage == 'external' or any(r[s + '_internal_extension'] for s in SIDES)))]
            for side in SIDES:
                lengths.append(dict(stage=stage, end=side, **describe([r[side + '_extension'] for r in selected if r[side + '_extension'] > 0])))
        summary['length_statistics'] = lengths
        # Gene/transcript counts aggregate safely because scopes are disjoint.
        summary['totals'] = {stage: {field: sum(g.get(field, 0) for g in groups if g['stage'] == stage)
                                     for field in metric_fields if not field.startswith('recorded_read_count_')}
                             for stage in STAGES}
        add_coverage_totals(summary, coverage)
        table(staging / 'genomic_coverage.csv', ['stage', 'chromosome', 'strand'] + COVERAGE_FIELDS, coverage)
        table(staging / 'newly_covered_intervals.csv', ['stage', 'chromosome', 'start', 'end', 'bases'], new_intervals)
        table(staging / 'overall_statistics.csv', ['stage'] + sorted(
            {k for values in summary['totals'].values() for k in values}),
            [dict(stage=stage, **values) for stage, values in summary['totals'].items()])
        (staging / 'summary.txt').write_text(written_summary(summary), encoding='utf-8')
        table(staging / 'support.csv', ['stage', 'chromosome', 'strand', 'end', 'support', 'transcript_ends'], support_rows)
        table(staging / 'extension_lengths.csv', ['stage', 'end', 'count', 'mean', 'median', 'p90', 'p95', 'maximum'], lengths)
        table(staging / 'diagnostics.csv', ['chromosome', 'strand', 'message'], errors)
        (staging / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
        (staging / 'index.html').write_text(render_dashboard(summary), encoding='utf-8')
        (staging / 'README.txt').write_text('Open index.html in a browser. All charts work offline. CSVs contain complete details.\n'
            'An incomplete report is not evidence of a successful run. Read summary.json and diagnostics.csv.\n'
            'These statistics do not repair or modify annotations. See summary.json definitions for counting rules.\n')
        staging.rename(output)
    print('{}: {}/{} partitions; {} LEAP models; {} external extended models'.format(summary['status'],
        len(provenance), len(scopes), summary['totals']['leap'].get('transcripts', 0), summary['totals']['external'].get('extended_transcripts', 0)))
    return summary['exit_code']


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument('--leap-output-dir', type=Path)
    source.add_argument('--manifest', type=Path, help='Declared partition/source manifest for workflow execution')
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--chromosomes', nargs='+')
    parser.add_argument('--strand', choices=('positive', 'negative', 'both'), default='both')
    parser.add_argument('--datacheck-summary', type=Path)
    parser.add_argument('--label', default='LEAP impact')
    parser.add_argument('--defer-errors', action='store_true', help='Publish diagnostics before a separate workflow gate fails')
    args = parser.parse_args()
    try:
        manifest = json.loads(args.manifest.read_text()) if args.manifest else discover(args.leap_output_dir, args.chromosomes, args.strand)
        summary = args.datacheck_summary
        if not summary and args.leap_output_dir:
            candidate = args.leap_output_dir / 'datachecks/summary/summary.json'
            if candidate.is_file():
                summary = candidate
        code = build(manifest, args.output, summary, args.label)
        return 0 if args.defer_errors else code
    except (OSError, EOFError, UnicodeError, ValueError, KeyError, TypeError, csv.Error) as error:
        print('ERROR: ' + str(error), file=sys.stderr)
        return 2


if __name__ == '__main__':
    sys.exit(main())
