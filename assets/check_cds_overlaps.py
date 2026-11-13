#!/usr/bin/env python3
"""Check LEAP overlap rules, plus the intentional within-gene partial-CDS check."""
import argparse
from bisect import bisect_right
from collections import Counter, defaultdict
import csv
import gzip
from pathlib import Path
import sys

UTRS = {'five_prime_UTR', 'three_prime_UTR'}
FIELDS = [
    'source_transcript',
    'source_gene',
    'target_transcript',
    'target_gene',
    'chromosome',
    'strand',
    'target_feature',
    'status',
    'added_overlap_intervals',
    'added_overlap_bases',
    'allowed_intervals',
    'allowed_bases',
    'failing_intervals',
    'failing_bases',
    'shared_splice_sites',
    'annotation_allowed_intervals',
    'splice_site_allowed_intervals',
    'annotation_extension_ends',
    'leap_allowed_intervals',
    'grouped_annotation_allowed_intervals',
    'grouped_annotation_boundaries',
    'grouped_annotation_donors',
    'check_rule',
    'gene_relationship',
    'target_is_partial',
    'target_missing_utr',
    'intentional_partial_cds_check',
]


def read_extension_stats(path):
    """Read LEAP's per-end provenance; never infer support from read counts."""
    path = Path(path)
    paths = (
        sorted(p for p in path.iterdir() if p.name.endswith(('.tsv', '.tsv.gz')))
        if path.is_dir()
        else [path]
    )
    if not paths:
        raise ValueError('{}: no extension stats TSV files found'.format(path))
    entries = {}
    required = {
        'gene_id',
        'transcript_id',
        'transcript_location',
        'five_prime_extension',
        'three_prime_extension',
        'five_prime_support',
        'three_prime_support',
    }
    for source in paths:
        opener = gzip.open if str(source).endswith('.gz') else open
        with opener(source, 'rt', encoding='utf-8', newline='') as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            if not required <= set(reader.fieldnames or []):
                raise ValueError('{}: missing extension stats columns'.format(source))
            for row in reader:
                context = '{}:{}'.format(source, reader.line_num)
                if any(not row.get(key) for key in required):
                    raise ValueError(context + ': incomplete extension stats row')
                identifier = row['transcript_id']
                try:
                    chromosome, interval = row['transcript_location'].rsplit(':', 1)
                    start, end = map(int, interval.split('-'))
                    if not chromosome or not 1 <= start <= end:
                        raise ValueError()
                except ValueError:
                    raise ValueError(
                        context + ': invalid transcript_location'
                    ) from None
                key = (chromosome, identifier)
                if key in entries:
                    raise ValueError(
                        context + ': duplicate stats transcript_id ' + identifier
                    )
                for end in ('five_prime', 'three_prime'):
                    try:
                        length = int(row[end + '_extension'])
                    except ValueError:
                        raise ValueError(
                            context + ': extension length must be an integer'
                        ) from None
                    if length < 0:
                        raise ValueError(
                            context + ': extension length must not be negative'
                        )
                    row[end + '_extension'] = length
                    support = row[end + '_support'].split(',')
                    if any(
                        token
                        not in ('current_annotation', 'reads', 'none', 'blacklisted')
                        for token in support
                    ):
                        raise ValueError(context + ': unknown extension support')
                # PAR transcript IDs can occur on both X and Y.
                entries[key] = row
    return entries


def validate_extension_stats(raw, extended, stats):
    """Reject stale provenance instead of exempting the wrong extension."""
    for identifier, source in extended.items():
        gene, chromosome, strand = source['location']
        entry = stats.get((chromosome, identifier))
        if entry is None:
            continue
        intervals = source['transcript'] or source['exon']
        start = min(interval[0] for interval in intervals)
        end = max(interval[1] for interval in intervals)
        if entry['gene_id'] != gene or entry[
            'transcript_location'
        ] != '{}:{}-{}'.format(chromosome, start, end):
            raise ValueError('extension stats gene/location differs for ' + identifier)
        if identifier not in raw or not source['exon'] or not raw[identifier]['exon']:
            continue
        original = raw[identifier]
        deltas = {
            'five_prime': max(0, original['exon'][0][0] - source['exon'][0][0]),
            'three_prime': max(0, source['exon'][-1][1] - original['exon'][-1][1]),
        }
        if strand == '-':
            deltas['five_prime'], deltas['three_prime'] = (
                deltas['three_prime'],
                deltas['five_prime'],
            )
        for side, length in deltas.items():
            if entry[side + '_extension'] != length:
                raise ValueError(
                    'extension stats {} length differs for {}'.format(side, identifier)
                )


def annotation_extensions(source, original, entry):
    """Return only the outward bases of ends supported purely by annotation."""
    if not entry or not source['exon'] or not original['exon']:
        return []
    strand = source['location'][2]
    result = []
    for side, interval in (
        (
            'five_prime' if strand == '+' else 'three_prime',
            (source['exon'][0][0], original['exon'][0][0] - 1),
        ),
        (
            'three_prime' if strand == '+' else 'five_prime',
            (original['exon'][-1][1] + 1, source['exon'][-1][1]),
        ),
    ):
        tokens = entry[side + '_support'].split(',')
        if interval[0] <= interval[1] and set(tokens) == {'current_annotation'}:
            result.append((interval, side))
    return result


def merge(intervals):
    result = []
    for start, end in sorted(intervals):
        if result and start <= result[-1][1] + 1:
            result[-1] = (result[-1][0], max(end, result[-1][1]))
        else:
            result.append((start, end))
    return result


def subtract(intervals, removed):
    """Subtract inclusive, merged removed intervals from intervals."""
    removed = merge(removed)
    for start, end in merge(intervals):
        cursor = start
        for left, right in removed:
            if right < cursor:
                continue
            if left > end:
                break
            if left > cursor:
                yield cursor, left - 1
            cursor = max(cursor, right + 1)
            if cursor > end:
                break
        if cursor <= end:
            yield cursor, end


def read_annotation(path):
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
            kind = fields[2]
            if kind not in UTRS | {'transcript', 'exon', 'CDS'}:
                continue
            attrs = {}
            for item in fields[8].split(';'):
                if not item:
                    continue
                key, sep, value = item.partition('=')
                if not sep or not key or key in attrs:
                    raise ValueError(context + ': malformed or duplicate attribute')
                attrs[key] = value
            for key in ('transcript_id', 'gene_id'):
                if not attrs.get(key, '').strip() or attrs[key] == '.':
                    raise ValueError(context + ': missing ' + key)
            try:
                start, end = int(fields[3]), int(fields[4])
            except ValueError:
                raise ValueError(context + ': coordinates must be integers') from None
            if start < 1 or end < start:
                raise ValueError(context + ': invalid coordinates')
            if fields[0] in ('', '.') or fields[6] not in ('+', '-'):
                raise ValueError(
                    context + ': chromosome and explicit + or - strand required'
                )
            identifier = attrs['transcript_id']
            location = (attrs['gene_id'], fields[0], fields[6])
            model = models.setdefault(
                identifier,
                dict(
                    location=location,
                    exon=set(),
                    CDS=set(),
                    utr=set(),
                    transcript=set(),
                    utr_types=set(),
                    transcript_type=None,
                    exon_numbers={},
                ),
            )
            biotype = attrs.get('transcript_type')
            if biotype:
                if model['transcript_type'] not in (None, biotype):
                    raise ValueError(context + ': inconsistent transcript_type')
                model['transcript_type'] = biotype
            if kind in UTRS:
                model['utr_types'].add(kind)
            if kind == 'exon' and attrs.get('exon_number'):
                try:
                    model['exon_numbers'][(start, end)] = int(attrs['exon_number'])
                except ValueError:
                    raise ValueError(context + ': invalid exon_number') from None
            if model['location'] != location:
                raise ValueError(
                    context + ': inconsistent gene/chromosome/strand for ' + identifier
                )
            if kind == 'transcript' and model['transcript']:
                raise ValueError(
                    context + ': duplicate transcript record for ' + identifier
                )
            model['utr' if kind in UTRS else kind].add((start, end))
    if not models:
        raise ValueError(str(path) + ': no transcript-associated records found')
    for identifier, model in models.items():
        model['exon'] = sorted(model['exon'])
        for previous, following in zip(model['exon'], model['exon'][1:]):
            if previous[1] >= following[0]:
                raise ValueError(
                    '{}: overlapping exons for {}'.format(path, identifier)
                )
        for interval in model['CDS'] | model['utr']:
            if not any(
                start <= interval[0] <= interval[1] <= end
                for start, end in model['exon']
            ):
                raise ValueError(
                    '{}: CDS/UTR lacks a containing exon for {}'.format(
                        path, identifier
                    )
                )
    return models


class IntervalIndex:
    """Sorted starts and prefix maximum ends prune non-overlapping records."""

    def __init__(self, records):
        self.records = sorted(records, key=lambda record: record[0])
        self.starts = [record[0] for record in self.records]
        self.ends = []
        maximum = 0
        for record in self.records:
            maximum = max(maximum, record[1])
            self.ends.append(maximum)

    def query(self, start, end):
        index = bisect_right(self.starts, end) - 1
        while index >= 0 and self.ends[index] >= start:
            record = self.records[index]
            if record[1] >= start:
                yield record
            index -= 1


def shared_site(source, raw_source, source_exon, target, target_exon, added):
    """Return evidence only for actual outward terminal-exon extension."""
    exons, old, targets = source['exon'], raw_source['exon'], target['exon']
    if len(exons) < 2 or len(old) < 2 or len(targets) < 2:
        return ''
    if (
        source_exon == exons[0]
        and source_exon[1] == old[0][1]
        and added[1] < old[0][0]
        and target_exon != targets[-1]
        and source_exon[1] == target_exon[1]
    ):
        return 'end:{}'.format(source_exon[1])
    if (
        source_exon == exons[-1]
        and source_exon[0] == old[-1][0]
        and added[0] > old[-1][1]
        and target_exon != targets[0]
        and source_exon[0] == target_exon[0]
    ):
        return 'start:{}'.format(source_exon[0])
    return ''


def format_intervals(intervals):
    return ';'.join('{}-{}'.format(*interval) for interval in intervals)


def full_length(model):
    return model['transcript_type'] == 'protein_coding' and model['utr_types'] >= UTRS


def terminal_kind(strand, side):
    return 'FE' if (side == 'left') == (strand == '+') else 'LE'


def terminal_safe(fe_genes, le_genes, natural):
    """Keep overlap_safety_checker.py's ordered zip and shared-group semantics."""
    if set(fe_genes) & set(le_genes):
        return 'shared_gene_group'
    if any(pair in natural for pair in zip(fe_genes, le_genes)):
        return 'recorded_natural_pair'
    return ''


class LeapContext:
    """Run-level provenance and the filtered guard regions used by LEAP."""

    def __init__(self, raw, extended, directory):
        import re

        self.raw, self.extended = raw, extended
        root = Path(directory)
        self.natural, self.candidates = {}, {}
        self.candidate_locations, self.candidate_files = {}, {}
        self.blockers, self.flanks, self.special = {}, {}, {}
        partitions = {model['location'][1:] for model in extended.values()}
        for chromosome, strand in partitions:
            label = 'positive' if strand == '+' else 'negative'
            path = root / 'natural_overlaps' / (chromosome + '_' + label + '_NO.txt')
            pairs = set()
            with path.open() as handle:
                for line in handle:
                    if not line.strip():
                        continue
                    fields = line.rstrip('\r\n').split('\t')
                    if len(fields) != 2 or not all(fields):
                        raise ValueError(str(path) + ': expected two ordered gene IDs')
                    pairs.add(tuple(fields))
            self.natural[(chromosome, strand)] = pairs
            for kind in ('FE', 'LE'):
                pattern = re.compile(
                    r'(?:^|_)' + re.escape(chromosome + '_' + label + '_' + kind) + r'_'
                )
                files = [
                    p
                    for p in (root / 'fixed_extensions').glob('*.bed')
                    if pattern.search(p.name)
                ]
                if len(files) != 1:
                    raise ValueError(
                        'Expected one fixed {} BED for {} {}'.format(
                            kind, chromosome, strand
                        )
                    )
                self.candidate_files[(chromosome, strand, kind)] = files[0]
                with files[0].open() as handle:
                    for number, line in enumerate(handle, 1):
                        fields = line.rstrip('\r\n').split('\t')
                        if (
                            len(fields) != 10
                            or fields[0] != chromosome
                            or fields[5] != strand
                        ):
                            raise ValueError(
                                str(files[0]) + ': invalid fixed-extension record'
                            )
                        start, end = int(fields[1]), int(fields[2])
                        if not 0 <= start <= end:
                            raise ValueError(str(files[0]) + ': invalid BED interval')
                        inner = end if (kind == 'FE') == (strand == '+') else start
                        genes = tuple(
                            (
                                fields[6] if fields[3] == 'single_exon' else fields[9]
                            ).split(',')
                        )
                        if not all(genes):
                            raise ValueError(str(files[0]) + ': missing gene IDs')
                        baseline = (
                            None
                            if fields[3] in ('.', 'single_exon')
                            else int(fields[3].split('/')[-1])
                        )
                        if baseline is not None and baseline < 1:
                            raise ValueError(
                                str(files[0]) + ': invalid original extension boundary'
                            )
                        # The writer overwrites its splice-site dictionary in file
                        # order. Support and group metadata must come from that same
                        # final row, even when earlier rows have identical geometry.
                        key = (chromosome, strand, kind, inner)
                        self.candidates[key] = (start, end, genes, fields[3] == '.', baseline)
                        self.candidate_locations[key] = (files[0], number)
        self.raw_terminals = defaultdict(list)
        for tid, model in raw.items():
            if model['transcript_type'] is None:
                raise ValueError(
                    'Raw transcript_type required for LEAP eligibility: ' + tid
                )
            gene, chromosome, strand = model['location']
            if (chromosome, strand) not in partitions or not full_length(model) or len(model['exon']) < 2:
                continue
            for side, exon in [('left', model['exon'][0]), ('right', model['exon'][-1])]:
                kind = terminal_kind(strand, side)
                inner, outer = (exon[1], exon[0]) if side == 'left' else (exon[0], exon[1])
                self.raw_terminals[(chromosome, strand, kind, inner)].append((outer, gene, tid))
        for partition in partitions:
            records, multi, single_cds = [], {}, defaultdict(set)
            for tid, model in raw.items():
                gene, chromosome, strand = model['location']
                if (chromosome, strand) != partition:
                    continue
                if model['transcript_type'] != 'protein_coding':
                    records.extend(
                        (a, b, tid, 'exon', 'noncoding_exon') for a, b in model['exon']
                    )
                elif full_length(model):
                    rule = (
                        'single_exon_cds'
                        if len(model['exon']) == 1
                        else 'full_length_cds'
                    )
                    records.extend((a, b, tid, 'CDS', rule) for a, b in model['CDS'])
                    if len(model['exon']) > 1:
                        multi[tid] = model
                        single_cds[gene].update(model['CDS'])
            natural_genes = {
                gene
                for pair in self.natural[partition]
                if pair[0] != pair[1]
                for gene in pair
            }
            special = {
                gene
                for gene, cds in single_cds.items()
                if len(cds) == 1 and gene in natural_genes
            }
            self.special[partition] = special
            self.blockers[partition] = IntervalIndex(records)
            for kind in ('FE', 'LE'):
                records = []
                for tid, model in multi.items():
                    if model['location'][0] not in special:
                        continue
                    exons = model['exon']
                    coding = [
                        exon
                        for exon in exons
                        if any(exon[0] <= a <= b <= exon[1] for a, b in model['CDS'])
                    ]
                    if not coding:
                        continue  # LEAP skips special-gene isoforms with no CDS row.
                    numbers = model['exon_numbers']
                    if any(exon not in numbers for exon in exons) or len(
                        {numbers[exon] for exon in exons}
                    ) != len(exons):
                        raise ValueError(
                            'Unambiguous exon_number required for LEAP single-CDS-gene flanks: '
                            + tid
                        )
                    coding_number = numbers[coding[0]]
                    selected = [
                        exon
                        for exon in exons
                        if (
                            numbers[exon] < coding_number
                            if kind == 'FE'
                            else numbers[exon] > coding_number
                        )
                    ] or coding
                    records.extend(
                        (a, b, tid, 'exon', 'single_cds_gene_flank')
                        for a, b in selected
                    )
                self.flanks[partition + (kind,)] = IntervalIndex(records)
        # The later overlap stage checks candidate FE/LE contacts and raw single exons.
        terminals = defaultdict(list)
        for tid, model in extended.items():
            if tid not in raw or not full_length(raw[tid]) or len(model['exon']) < 2:
                continue
            gene, chromosome, strand = model['location']
            for side, exon in [
                ('left', model['exon'][0]),
                ('right', model['exon'][-1]),
            ]:
                kind = terminal_kind(strand, side)
                candidate = self.candidate(model, exon, kind, required=False)
                if candidate:
                    terminals[(chromosome, strand, kind)].append(
                        (*exon, tid, candidate[0], kind)
                    )
        for tid, model in raw.items():
            gene, chromosome, strand = model['location']
            if (
                (chromosome, strand) in partitions
                and full_length(model)
                and len(model['exon']) == 1
            ):
                terminals[(chromosome, strand, 'SE')].append(
                    (*model['exon'][0], tid, (gene,), 'SE')
                )
        self.terminals = {key: IntervalIndex(value) for key, value in terminals.items()}

    def candidate(self, source, exon, kind, required=True):
        """Select the writer's last row first, then validate its applicability."""
        gene, chromosome, strand = source['location']
        inner = exon[1] if (kind == 'FE') == (strand == '+') else exon[0]
        key = (chromosome, strand, kind, inner)
        row = self.candidates.get(key)
        path = self.candidate_files[(chromosome, strand, kind)]
        location = self.candidate_locations.get(key)
        context = '{}:{}'.format(*location) if location else str(path)
        description = '{} {} {} {} splice_site={}'.format(gene, chromosome, strand, kind, inner)
        if row is None:
            if not required:
                return None
            raise ValueError('{}: missing fixed-extension context for {}: no record for splice site'.format(
                context, description))
        if gene not in row[2]:
            raise ValueError('{}: incompatible final fixed-extension record for {}: selected gene group {} omits source gene'.format(
                context, description, ','.join(row[2])))
        if required and not (row[0] <= exon[0] and exon[1] <= row[1]):
            raise ValueError('{}: incompatible final fixed-extension record for {}: selected interval {}-{} does not contain output exon {}-{}'.format(
                context, description, row[0], row[1], exon[0], exon[1]))
        # An unchanged terminal exon may exceed a shorter selected candidate:
        # LEAP never shrinks it. Processed extensions may lie inside the candidate.
        return row[2], row[3], row[4]

    def grouped_annotation_parts(self, source, original, ends):
        """Added isoform bases already covered by its original eligible splice group.

        Read-support labels do not turn pre-existing annotation into new read
        sequence. Recorded baselines cap the allowance, and raw full-length
        multi-exon donors verify coverage; partial/internal exons cannot supply it.
        """
        gene, chromosome, strand = source['location']
        for exon, kind, (genes, annotation, recorded), guard, added in ends:
            left = (kind == 'FE') == (strand == '+')
            inner = exon[1] if left else exon[0]
            old = original['exon'][0 if left else -1]
            if inner != (old[1] if left else old[0]):
                continue  # A changed splice site cannot acquire another raw group's allowance.
            donors = [row for row in self.raw_terminals.get((chromosome, strand, kind, inner), [])
                      if row[1] in genes]
            if not donors:
                continue
            boundary = (min if left else max)(row[0] for row in donors)
            if recorded is not None:
                boundary = max(boundary, recorded) if left else min(boundary, recorded)
            interval = (max(added[0], boundary), added[1]) if left else (added[0], min(added[1], boundary))
            if interval[0] > interval[1]:
                continue
            supporting = tuple(sorted(row[2] for row in donors if
                                      (row[0] <= boundary if left else row[0] >= boundary)))
            yield interval, ('five_prime' if kind == 'FE' else 'three_prime'), boundary, supporting

    def source_ends(self, source, original):
        if not source['exon'] or not original['exon']:
            return
        for side, exon, old in [
            ('left', source['exon'][0], original['exon'][0]),
            ('right', source['exon'][-1], original['exon'][-1]),
        ]:
            grew = exon[0] < old[0] if side == 'left' else exon[1] > old[1]
            if not grew:
                continue
            kind = terminal_kind(source['location'][2], side)
            candidate = self.candidate(source, exon, kind)
            # The CDS clipping helper includes the old outer boundary; bedtools does not count touching endpoints.
            baseline = (
                candidate[2]
                if candidate[2] is not None
                else (old[0] if side == 'left' else old[1])
            )
            guard = (exon[0], baseline) if side == 'left' else (baseline, exon[1])
            if guard[0] >= guard[1]:
                guard = (
                    None  # This isoform grew only to an existing annotation boundary.
                )
            added = (exon[0], old[0] - 1) if side == 'left' else (old[1] + 1, exon[1])
            yield exon, kind, candidate, guard, added


def compare(raw, extended, counts, extension_stats=None, leap_context=None):
    """Yield rule-specific findings. Intervals are inclusive report coordinates."""
    if leap_context is None:
        raise ValueError('LEAP run context is required for cross-gene checks')
    extension_stats = {} if extension_stats is None else extension_stats
    validate_extension_stats(raw, extended, extension_stats)
    records = defaultdict(list)
    for identifier, model in raw.items():
        for kind in ('exon', 'CDS'):
            for interval in model[kind]:
                exon = next(
                    exon
                    for exon in model['exon']
                    if exon[0] <= interval[0] <= interval[1] <= exon[1]
                )
                records[model['location'][1:]].append(
                    (*interval, identifier, kind, exon)
                )
    indexes = {key: IntervalIndex(value) for key, value in records.items()}
    for identifier, source in sorted(extended.items()):
        if not source['utr']:
            continue
        gene, chromosome, strand = source['location']
        if identifier not in raw:
            counts['missing_from_raw'] += 1
            yield dict(
                source_transcript=identifier,
                source_gene=gene,
                chromosome=chromosome,
                strand=strand,
                status='missing_from_raw',
            )
            continue
        original = raw[identifier]
        if original['location'] != source['location']:
            raise ValueError('raw/extended gene or location differs for ' + identifier)
        counts['evaluated_transcripts'] += 1
        partition = (chromosome, strand)
        annotation_parts = annotation_extensions(
            source, original, extension_stats.get((chromosome, identifier))
        )
        source_ends = list(leap_context.source_ends(source, original))
        grouped_parts = list(leap_context.grouped_annotation_parts(source, original, source_ends))
        added_utrs = list(subtract(source['utr'], original['utr']))
        groups, audit, classified = {}, defaultdict(list), defaultdict(list)

        def add(target_id, kind, rule, overlap, outcome, site='', end='', boundary='', donors=()):
            group = groups.setdefault(
                (target_id, kind, rule),
                dict(
                    annotation=[],
                    splice=[],
                    failing=[],
                    leap=[],
                    grouped=[],
                    boundaries=set(),
                    donors=set(),
                    sites=set(),
                    ends=set(),
                ),
            )
            group[outcome].append(overlap)
            if boundary:
                group['boundaries'].add(boundary)
            group['donors'].update(donors)
            if site:
                group['sites'].add(site)
            if end:
                group['ends'].add(end)
            classified[target_id].append(overlap)

        # The intentional within-gene check includes raw partial-model CDS.
        for added in added_utrs:
            for source_exon in source['exon']:
                piece = (max(added[0], source_exon[0]), min(added[1], source_exon[1]))
                if piece[0] > piece[1]:
                    continue
                for start, end, target_id, kind, target_exon in indexes.get(
                    partition, IntervalIndex([])
                ).query(*piece):
                    if target_id == identifier:
                        continue
                    target = raw[target_id]
                    same_gene = gene == target['location'][0]
                    overlap = (max(piece[0], start), min(piece[1], end))
                    if not same_gene:
                        if kind == 'exon':
                            audit[(target_id, 'exon')].append(overlap)
                        continue
                    if kind != 'CDS':
                        continue
                    if target['transcript_type'] != 'protein_coding' or full_length(
                        target
                    ):
                        audit[(target_id, 'CDS')].append(overlap)
                        continue
                    exempt = []
                    for interval, side, boundary, donors in grouped_parts:
                        part = (max(overlap[0], interval[0]), min(overlap[1], interval[1]))
                        if part[0] <= part[1]:
                            exempt.append(part)
                            add(target_id, kind, 'within_gene_partial_cds', part, 'grouped',
                                boundary='{}:{}'.format(side, boundary), donors=donors)
                    for interval, side in annotation_parts:
                        part = (
                            max(overlap[0], interval[0]),
                            min(overlap[1], interval[1]),
                        )
                        if part[0] <= part[1]:
                            exempt.append(part)
                            add(
                                target_id,
                                kind,
                                'within_gene_partial_cds',
                                part,
                                'annotation',
                                end=side,
                            )
                    site = shared_site(
                        source, original, source_exon, target, target_exon, piece
                    )
                    for part in subtract([overlap], exempt):
                        add(
                            target_id,
                            kind,
                            'within_gene_partial_cds',
                            part,
                            'splice' if site else 'failing',
                            site=site,
                        )

        # LEAP clips against eligible features from every gene, including the source gene.
        # Its grouped candidate metadata, rather than per-transcript stats, controls bypasses.
        for (
            exon,
            kind,
            (genes, annotation, baseline),
            guard,
            added,
        ) in source_ends:
            blockers = [leap_context.blockers[partition]]
            if len(genes) == 1 and genes[0] in leap_context.special[partition]:
                blockers.append(leap_context.flanks[partition + (kind,)])
            for index in blockers if guard else []:
                for start, end, target_id, feature, rule in index.query(*guard):
                    if target_id == identifier:
                        continue
                    overlap = (max(start, guard[0]), min(end, guard[1]))
                    side = 'five_prime' if kind == 'FE' else 'three_prime'
                    add(
                        target_id,
                        feature,
                        rule,
                        overlap,
                        'annotation' if annotation else 'failing',
                        end=side if annotation else '',
                    )
            opposite = 'LE' if kind == 'FE' else 'FE'
            for target_kind in (opposite, 'SE'):
                index = leap_context.terminals.get(partition + (target_kind,))
                if index is None:
                    continue
                for start, end, target_id, target_genes, _ in index.query(*exon):
                    if raw[target_id]['location'][0] == gene or max(
                        start, exon[0]
                    ) >= min(end, exon[1]):
                        continue  # bedtools requires positive-length BED overlap.
                    # LEAP tests the whole candidate exon, not just its newly added bases.
                    overlap = (max(start, exon[0]), min(end, exon[1]))
                    fe_genes, le_genes = (
                        (genes, target_genes) if kind == 'FE' else (target_genes, genes)
                    )
                    reason = terminal_safe(
                        fe_genes, le_genes, leap_context.natural[partition]
                    )
                    rule = (
                        'terminal_single_exon'
                        if target_kind == 'SE'
                        else 'terminal_fe_le'
                    )
                    add(
                        target_id,
                        'exon',
                        rule + (':' + reason if reason else ':unknown_pair'),
                        overlap,
                        'leap' if reason else 'failing',
                    )
        # Retain non-failing raw exon contacts as an audit trail, without duplicating classified bases.
        for (target_id, audit_feature), intervals in audit.items():
            target = raw[target_id]
            reason = (
                'ineligible_partial_model'
                if target['transcript_type'] == 'protein_coding'
                and not full_length(target)
                else (
                    'outside_leap_extension_guard'
                    if audit_feature == 'CDS'
                    else 'outside_leap_protected_regions'
                )
            )
            for part in subtract(intervals, classified[target_id]):
                add(target_id, audit_feature, reason, part, 'leap')
        for (target_id, feature, rule), group in sorted(groups.items()):
            failing = merge(group['failing'])
            grouped = list(subtract(group['grouped'], failing))
            annotation = list(subtract(group['annotation'], failing + grouped))
            splice = list(subtract(group['splice'], failing + grouped + annotation))
            leap = list(subtract(group['leap'], failing + grouped + annotation + splice))
            allowed = merge(grouped + annotation + splice + leap)
            same_gene = gene == raw[target_id]['location'][0]
            if failing:
                status = (
                    (
                        'new_cds_overlap'
                        if feature == 'CDS'
                        else 'new_within_gene_exon_overlap'
                    )
                    if same_gene
                    else (
                        'unsafe_terminal_overlap'
                        if rule.startswith('terminal_')
                        else (
                            'new_cross_gene_cds_overlap'
                            if feature == 'CDS'
                            else 'new_cross_gene_exon_overlap'
                        )
                    )
                )
            elif leap:
                status = 'allowed_leap_overlap'
            elif sum(bool(parts) for parts in (grouped, annotation, splice)) > 1:
                status = 'allowed_multiple_reasons'
            elif grouped:
                status = 'allowed_grouped_annotation'
            else:
                status = (
                    'allowed_current_annotation'
                    if annotation
                    else 'allowed_shared_splice_site'
                )
            counts[status] += 1
            if failing:
                counts[
                    (
                        'partial_cds_failures'
                        if rule == 'within_gene_partial_cds'
                        else 'leap_rule_failures'
                    )
                ] += 1
            target = raw[target_id]
            partial = target['transcript_type'] == 'protein_coding' and not full_length(
                target
            )
            row = dict(
                source_transcript=identifier,
                source_gene=gene,
                target_transcript=target_id,
                target_gene=target['location'][0],
                chromosome=chromosome,
                strand=strand,
                target_feature=feature,
                status=status,
                shared_splice_sites=';'.join(sorted(group['sites'])),
                annotation_allowed_intervals=format_intervals(annotation),
                splice_site_allowed_intervals=format_intervals(splice),
                leap_allowed_intervals=format_intervals(leap),
                grouped_annotation_allowed_intervals=format_intervals(grouped),
                grouped_annotation_boundaries=';'.join(sorted(group['boundaries'])),
                grouped_annotation_donors=';'.join(sorted(group['donors'])),
                check_rule=rule,
                gene_relationship='within_gene' if same_gene else 'cross_gene',
                target_is_partial=str(partial).lower(),
                intentional_partial_cds_check=str(
                    rule == 'within_gene_partial_cds'
                ).lower(),
                target_missing_utr=(
                    ';'.join(sorted(UTRS - target['utr_types'])) if partial else ''
                ),
                annotation_extension_ends=';'.join(sorted(group['ends'])),
            )
            for prefix, intervals in (
                ('added_overlap', merge(allowed + failing)),
                ('allowed', allowed),
                ('failing', failing),
            ):
                row[prefix + '_intervals'] = format_intervals(intervals)
                row[prefix + '_bases'] = sum(
                    end - start + 1 for start, end in intervals
                )
            yield row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('current_annotation')
    parser.add_argument('extended_annotation')
    parser.add_argument('output', help='CSV report file (not a directory)')
    parser.add_argument(
        '--leap-output-dir',
        required=True,
        help='Original LEAP run directory containing natural_overlaps/ and fixed_extensions/',
    )
    parser.add_argument(
        '--extension-stats',
        help='Validated support TSV file or directory; enables the partial-CDS annotation allowance',
    )
    args = parser.parse_args()
    try:
        output = Path(args.output)
        sources = [args.current_annotation, args.extended_annotation]
        if args.extension_stats:
            stats_path = Path(args.extension_stats)
            if stats_path.is_dir():
                sources.extend(str(p) for p in stats_path.iterdir() if p.is_file())
            else:
                sources.append(args.extension_stats)
        context_root = Path(args.leap_output_dir)
        sources.extend(
            str(p)
            for folder in ('natural_overlaps', 'fixed_extensions')
            for p in (context_root / folder).glob('*')
            if p.is_file()
        )
        for source in sources:
            if output.resolve() == Path(source).resolve() or (
                output.exists() and output.samefile(source)
            ):
                raise ValueError('output must not overwrite an input annotation')
        raw, extended = read_annotation(args.current_annotation), read_annotation(
            args.extended_annotation
        )
        extension_stats = (
            read_extension_stats(args.extension_stats) if args.extension_stats else {}
        )
        validate_extension_stats(raw, extended, extension_stats)
        # Check matched locations before opening the report.
        for identifier in raw.keys() & extended.keys():
            if raw[identifier]['location'] != extended[identifier]['location']:
                raise ValueError(
                    'raw/extended gene or location differs for ' + identifier
                )
        context = LeapContext(raw, extended, args.leap_output_dir)
        for identifier in raw.keys() & extended.keys():
            list(context.source_ends(extended[identifier], raw[identifier]))
        counts = Counter()
        with output.open('w', newline='', encoding='utf-8') as handle:
            writer = csv.DictWriter(handle, fieldnames=FIELDS)
            writer.writeheader()
            writer.writerows(compare(raw, extended, counts, extension_stats, context))
        failures = counts['leap_rule_failures'] + counts['partial_cds_failures']
        print(
            '{}: {} failing overlap comparison(s)'.format(
                'FAIL' if failures else 'PASS', failures
            )
        )
        print(
            'partial_cds_current_annotation_exception: '
            + ('enabled' if args.extension_stats else 'disabled (no --extension-stats)')
        )
        print('leap_rule_context: ' + str(Path(args.leap_output_dir).resolve()))
        for key in (
            'evaluated_transcripts',
            'leap_rule_failures',
            'partial_cds_failures',
            'allowed_shared_splice_site',
            'allowed_current_annotation',
            'allowed_grouped_annotation',
            'allowed_multiple_reasons',
            'new_cds_overlap',
            'new_within_gene_exon_overlap',
            'new_cross_gene_cds_overlap',
            'new_cross_gene_exon_overlap',
            'unsafe_terminal_overlap',
            'allowed_leap_overlap',
            'missing_from_raw',
        ):
            print('{}: {}'.format(key, counts[key]))
        return 1 if failures else 0
    except (OSError, EOFError, UnicodeError, ValueError, csv.Error) as error:
        print('ERROR: {}'.format(error), file=sys.stderr)
        return 2


if __name__ == '__main__':
    sys.exit(main())
