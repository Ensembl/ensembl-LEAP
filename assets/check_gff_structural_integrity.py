#!/usr/bin/env python3
"""Validate LEAP GFF3 hierarchy, identity, exon structure and feature containment."""

import argparse
import csv
import gzip
import hashlib
import json
import math
import re
import sys
from collections import Counter, defaultdict, deque
from pathlib import Path
from typing import NamedTuple
from urllib.parse import unquote


TRANSCRIPTS = {'transcript', 'mRNA'}
EXON_FEATURES = {'CDS', 'start_codon', 'stop_codon', 'five_prime_UTR', 'three_prime_UTR'}
FIELDS = ['issue', 'line', 'feature_type', 'feature_id', 'transcript_id', 'gene_id',
          'chromosome', 'start', 'end', 'strand', 'parent_id', 'details']


class Feature(NamedTuple):
    line: int
    chromosome: str
    kind: str
    start: int
    end: int
    strand: str
    phase: str
    identifier: str
    parents: tuple
    transcripts: tuple
    genes: tuple


def identifiers(value, context):
    """Split before percent decoding, so an encoded comma remains part of an ID."""
    if value is None:
        return ()
    values = value.split(',')
    if any(not v or v == '.' or re.search(r'%(?![0-9A-Fa-f]{2})', v) for v in values):
        raise ValueError(context + ': empty or malformed identifier list')
    return tuple(sys.intern(unquote(v, errors='strict')) for v in values)


def parse_annotation(path):
    features, groups, duplicates, seen = [], defaultdict(list), [], {}
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
            if not 1 <= start <= end or f[0] in ('', '.') or not f[1] or f[2] in ('', '.'):
                raise ValueError(context + ': invalid coordinates, chromosome, source or feature type')
            if f[6] not in ('+', '-', '.', '?') or f[7] not in ('.', '0', '1', '2'):
                raise ValueError(context + ': invalid strand or phase')
            if f[5] != '.':
                try:
                    if not math.isfinite(float(f[5])):
                        raise ValueError()
                except ValueError:
                    raise ValueError(context + ': invalid score') from None
            attrs = {}
            if f[8] != '.':
                for item in f[8].split(';'):
                    if not item:
                        continue
                    key, sep, value = item.partition('=')
                    if not sep or not key or key in attrs:
                        raise ValueError(context + ': malformed or duplicate attribute')
                    attrs[key] = value
            ids = identifiers(attrs.get('ID'), context)
            if len(ids) > 1:
                raise ValueError(context + ': ID must have one value; encode literal commas')
            parents = identifiers(attrs.get('Parent'), context)
            tids = identifiers(attrs.get('transcript_id'), context)
            gids = identifiers(attrs.get('gene_id'), context)
            feature = Feature(number, sys.intern(f[0]), sys.intern(f[2]), start, end, f[6], f[7],
                              ids[0] if ids else '', parents, tids, gids)
            # Canonicalize attribute order and Parent order for duplicate detection.
            canonical = dict(attrs)
            if parents:
                canonical['Parent'] = tuple(sorted(parents))
            signature = hashlib.sha256(repr((tuple(f[:3]), start, end, tuple(f[5:8]), sorted(canonical.items()))).encode()).digest()
            if signature in seen:
                duplicates.append((feature, seen[signature]))
            else:
                seen[signature] = number
            features.append(feature)
            if feature.identifier:
                groups[feature.identifier].append(feature)
    if not features:
        raise ValueError(str(path) + ': no feature records found')
    if not any(feature.kind in TRANSCRIPTS for feature in features):
        raise ValueError(str(path) + ': no transcript/mRNA records found')
    return features, dict(groups), duplicates


def merged(intervals):
    result = []
    for start, end in sorted(intervals):
        if result and start <= result[-1][1] + 1:
            result[-1] = (result[-1][0], max(end, result[-1][1]))
        else:
            result.append((start, end))
    return result


def check_annotation(parsed, allow_missing_genes=False):
    features, groups, duplicates = parsed
    issues = []

    def issue(code, feature, parent='', **details):
        issues.append(dict(issue=code, line=feature.line, feature_type=feature.kind,
                           feature_id=feature.identifier, transcript_id=','.join(feature.transcripts),
                           gene_id=','.join(feature.genes), chromosome=feature.chromosome,
                           start=feature.start, end=feature.end, strand=feature.strand,
                           parent_id=parent, details=json.dumps(details, sort_keys=True, separators=(',', ':'))))

    genes = [f for f in features if f.kind == 'gene']
    transcripts = [f for f in features if f.kind in TRANSCRIPTS]
    omit_genes = allow_missing_genes and not genes
    bad_ids = set()
    for feature, first in duplicates:
        issue('duplicate_record', feature, first_line=first)
    for identifier, parts in groups.items():
        first = parts[0]
        metadata = {(p.chromosome, p.strand, p.kind, tuple(sorted(p.parents)),
                     tuple(sorted(p.genes)), tuple(sorted(p.transcripts))) for p in parts}
        if len(parts) > 1 and (first.kind in TRANSCRIPTS | {'gene'} or len(metadata) > 1):
            bad_ids.add(identifier)
            issue('conflicting_id', parts[1], first_line=first.line,
                  reason='gene/transcript IDs must be unique; multipart feature identity must agree',
                  record_lines=[p.line for p in parts])
    for feature in features:
        if feature.kind in TRANSCRIPTS | {'gene'} and not feature.identifier:
            issue('missing_id', feature)
        if len(set(feature.parents)) != len(feature.parents):
            issue('duplicate_parent_reference', feature)
        if feature.identifier in feature.parents:
            issue('self_parent', feature, feature.identifier)
        if feature.kind == 'gene' and feature.parents:
            issue('invalid_parent_type', feature, reason='gene must be a root feature')
        if feature.kind in TRANSCRIPTS and len(feature.parents) != 1:
            issue('invalid_transcript_parent_count', feature, observed=len(feature.parents), expected=1)
        if (feature.kind in TRANSCRIPTS | EXON_FEATURES | {'exon'} or feature.transcripts or feature.genes) and feature.kind != 'gene' and not feature.parents:
            issue('missing_parent', feature)
        if feature.kind in TRANSCRIPTS and len(feature.transcripts) > 1:
            issue('transcript_id_mismatch', feature, reason='transcript row must describe one transcript')
        if feature.kind == 'gene' and len(feature.genes) > 1:
            issue('gene_id_mismatch', feature, reason='gene row must describe one gene')
        if feature.kind == 'CDS' and feature.phase == '.':
            issue('missing_cds_phase', feature)
    for records, attribute, code in [(transcripts, 'transcripts', 'duplicate_transcript_id'),
                                     (genes, 'genes', 'duplicate_gene_id')]:
        identities = {}
        for feature in records:
            for identifier in getattr(feature, attribute):
                key = (feature.chromosome, identifier)
                if key in identities:
                    issue(code, feature, first_line=identities[key])
                else:
                    identities[key] = feature.line

    # Resolve ID dependencies without recursion. Bad identities and cycles never
    # silently select one of several incompatible parent definitions.
    usable = {key:parts for key,parts in groups.items() if key not in bad_ids}
    dependencies, children = {}, defaultdict(set)
    for key, parts in usable.items():
        dependencies[key] = {p for p in parts[0].parents if p in usable}
        for parent in dependencies[key]:
            children[parent].add(key)
    pending = {key:len(value) for key,value in dependencies.items()}
    ready = deque(key for key,n in pending.items() if n == 0)
    ancestors = {}
    while ready:
        key = ready.popleft();feature = usable[key][0]
        tx, gn = set(), set()
        for parent in dependencies[key]:
            pt, pg = ancestors[parent];tx.update(pt);gn.update(pg)
        if feature.kind in TRANSCRIPTS:
            tx.add(key)
        if feature.kind == 'gene':
            gn.add(key)
        ancestors[key] = (frozenset(tx), frozenset(gn))
        for child in children[key]:
            pending[child] -= 1
            if pending[child] == 0:
                ready.append(child)
    blocked = set(usable) - ancestors.keys()
    for key in sorted(blocked):
        issue('cyclic_parent_graph', usable[key][0], reason='cycle or dependency on a cycle')

    intervals = {key:merged((p.start,p.end) for p in parts) for key,parts in usable.items()}
    tx_exons, tx_features, gene_transcripts = defaultdict(list), defaultdict(list), defaultdict(set)
    skipped_gene_parents = 0
    # Custom Ensembl IDs are checked against resolved ancestry, not substituted
    # for GFF IDs when resolving Parent links.
    for feature in features:
        tx_ids, gene_ids = set(), set()
        for parent in set(feature.parents):
            if parent in bad_ids:
                issue('ambiguous_parent', feature, parent)
                continue
            if parent not in usable:
                if omit_genes and feature.kind in TRANSCRIPTS and len(feature.parents) == 1:
                    skipped_gene_parents += 1
                    continue
                issue('unresolved_parent', feature, parent)
                continue
            parts = usable[parent];p = parts[0]
            if feature.kind in TRANSCRIPTS and p.kind != 'gene':
                issue('invalid_parent_type', feature, parent, expected='gene', observed=p.kind)
            if feature.kind == 'exon' and p.kind not in TRANSCRIPTS:
                issue('invalid_parent_type', feature, parent, expected='transcript/mRNA', observed=p.kind)
            if feature.chromosome != p.chromosome:
                issue('chromosome_mismatch', feature, parent, parent_chromosome=p.chromosome)
            if feature.strand != p.strand:
                issue('strand_mismatch', feature, parent, parent_strand=p.strand)
            if not any(a <= feature.start <= feature.end <= b for a,b in intervals[parent]):
                issue('outside_parent', feature, parent, parent_intervals=intervals[parent])
            if parent in ancestors:
                pt, pg = ancestors[parent];tx_ids.update(pt);gene_ids.update(pg)
        if feature.kind in TRANSCRIPTS and feature.identifier:
            tx_ids.add(feature.identifier)
        if feature.kind == 'gene' and feature.identifier:
            gene_ids.add(feature.identifier)
        if feature.kind in EXON_FEATURES | {'exon'} and not tx_ids and not blocked.intersection(feature.parents):
            issue('missing_transcript_ancestor', feature)
        expected_tids = {t for key in tx_ids if key in usable for t in (usable[key][0].transcripts or (key,))}
        expected_gids = {g for key in gene_ids if key in usable for g in (usable[key][0].genes or (key,))}
        if omit_genes:
            expected_gids.update(g for key in tx_ids if key in usable for g in usable[key][0].genes)
        if feature.transcripts and set(feature.transcripts) != expected_tids:
            issue('transcript_id_mismatch', feature, observed=feature.transcripts, expected=sorted(expected_tids))
        if feature.genes and set(feature.genes) != expected_gids:
            issue('gene_id_mismatch', feature, observed=feature.genes, expected=sorted(expected_gids))
        # Explicit ancestor checks identify every outlying descendant, including
        # those attached through another (possibly malformed) child feature.
        for ancestor in (tx_ids | gene_ids) - {feature.identifier}:
            if ancestor not in usable:
                continue
            parent = usable[ancestor][0]
            if feature.chromosome != parent.chromosome or feature.strand != parent.strand:
                # Direct-parent mismatch already has a precise diagnostic.
                if ancestor not in feature.parents:
                    issue('ancestor_location_mismatch', feature, ancestor)
            if ancestor not in feature.parents and not any(a <= feature.start <= feature.end <= b for a,b in intervals[ancestor]):
                issue('outside_ancestor', feature, ancestor, ancestor_type=parent.kind,
                      ancestor_intervals=intervals[ancestor])
        if feature.kind == 'exon':
            for tx in tx_ids:
                tx_exons[tx].append(feature)
        elif feature.kind in EXON_FEATURES:
            for tx in tx_ids:
                tx_features[tx].append(feature)
        if feature.kind in TRANSCRIPTS:
            for gene in gene_ids:
                gene_transcripts[gene].add(feature.identifier)

    for tx in transcripts:
        exons = sorted(tx_exons[tx.identifier], key=lambda f:(f.start,f.end,f.line))
        if not exons:
            issue('transcript_without_exons', tx)
        for a,b in zip(exons,exons[1:]):
            if a.end >= b.start:
                issue('overlapping_exons', b, tx.identifier, previous_line=a.line)
        for feature in tx_features[tx.identifier]:
            if not any(exon.chromosome == feature.chromosome and exon.strand == feature.strand
                       and exon.start <= feature.start <= feature.end <= exon.end for exon in exons):
                issue('feature_outside_exon', feature, tx.identifier)
    for gene in genes:
        if gene.identifier and not gene_transcripts[gene.identifier]:
            issue('gene_without_transcripts', gene)
    if omit_genes:
        aliases = defaultdict(list)
        for tx in transcripts:
            if len(tx.parents) == 1 and tx.parents[0] not in groups and tx.genes:
                aliases[tx.parents[0]].append(tx)
        for parent, members in aliases.items():
            if len({g for tx in members for g in tx.genes}) > 1:
                for tx in members:
                    issue('gene_id_mismatch', tx, parent,
                          reason='gene-free transcripts sharing a gene Parent declare different gene IDs')
    counts = dict(features=len(features), genes=len(genes), transcripts=len(transcripts),
                  exons=sum(f.kind=='exon' for f in features),
                  omitted_gene_parent_checks=skipped_gene_parents,
                  gene_boundaries_checked=bool(genes))
    return issues, counts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('annotation', help='One GFF3 annotation, optionally gzipped')
    parser.add_argument('output', help='CSV issue report filename')
    parser.add_argument('--allow-missing-genes', action='store_true',
                        help='Allow intentionally gene-free files; all other parent links remain required')
    args = parser.parse_args()
    try:
        output, source = Path(args.output), Path(args.annotation)
        if output.resolve() == source.resolve() or (output.exists() and output.samefile(source)):
            raise ValueError('output must not overwrite the input annotation')
        parsed = parse_annotation(source)
        issues, counts = check_annotation(parsed, args.allow_missing_genes)
        with output.open('w', newline='', encoding='utf-8') as handle:
            writer = csv.DictWriter(handle, fieldnames=FIELDS);writer.writeheader();writer.writerows(issues)
        print('{}: {} structural issue(s)'.format('FAIL' if issues else 'PASS', len(issues)))
        for key,value in counts.items():
            print('{}: {}'.format(key,value))
        for code,value in sorted(Counter(r['issue'] for r in issues).items()):
            print('{}: {}'.format(code,value))
        return int(bool(issues))
    except (OSError, EOFError, UnicodeError, ValueError, csv.Error) as error:
        print('ERROR: ' + str(error), file=sys.stderr)
        return 2


if __name__ == '__main__':
    sys.exit(main())
