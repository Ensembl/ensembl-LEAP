"""Preserve the Perl exporter's chromosome-Y PAR annotation projection."""

from collections import defaultdict
from dataclasses import replace
from pathlib import Path
import re
import sys
import tempfile

from .annotation import read_features
from .regions import fetch_regions


def _par_ids(features, fmt):
    result = []
    for feature in features:
        attrs = dict(feature.attrs)
        gene, transcript = attrs.get('gene_id'), attrs.get('transcript_id')
        def suffixed(value):
            return value if value.endswith('_PAR_Y') else value + '_PAR_Y'
        if fmt == 'gff3':
            for key in ('ID', 'Parent'):
                value = attrs.get(key)
                if not value:
                    continue
                # Replace the embedded transcript ID, not the child exon rank.
                for identifier in (transcript, gene):
                    if identifier:
                        value = re.sub(re.escape(identifier) + r'(?=:|,|$)', lambda match: suffixed(match[0]), value)
                attrs[key] = value
        else:
            for key in ('gene_id', 'transcript_id'):
                if attrs.get(key): attrs[key] = suffixed(attrs[key])
        if feature.feature_type == 'gene':
            tags = attrs.get('tag', '').split(',') if attrs.get('tag') else []
            if 'PAR' not in tags: tags.append('PAR')
            attrs['tag'] = ','.join(tags)
        result.append(replace(feature, attrs=attrs))
    return result


def add_y_par_features(secrets, args, regions, features, export):
    hosts = [r for r in regions if r['name'] == 'Y']
    if not hosts:
        return features
    from sqlalchemy import bindparam, text
    from .exporter import sql_connection
    mappings = []
    gene_order = {}
    with sql_connection(secrets) as connection:
        for host in hosts:
            rows = connection.execute(text('''
                SELECT a.seq_region_start AS host_start, a.seq_region_end AS host_end,
                    a.exc_seq_region_id AS source_id, a.exc_seq_region_start AS source_start,
                    a.exc_seq_region_end AS source_end, a.ori, s.name AS source_name
                FROM assembly_exception a JOIN seq_region s ON s.seq_region_id=a.exc_seq_region_id
                WHERE a.seq_region_id=:host_id AND a.exc_type='PAR'
                ORDER BY a.seq_region_start
            '''), dict(host_id=host['seq_region_id']))
            for row in rows.mappings():
                mapping = dict(row)
                resolved = fetch_regions(connection, chromosome=mapping['source_name'])
                if len(resolved) != 1 or resolved[0]['seq_region_id'] != mapping['source_id']:
                    raise ValueError('Y PAR counterpart is not an unambiguous default-assembly region')
                if mapping['ori'] not in (1, -1) or mapping['host_end']-mapping['host_start'] != mapping['source_end']-mapping['source_start']:
                    raise ValueError('Invalid Y PAR projection interval')
                mappings.append(mapping)
        if not mappings:
            return features
        ids = tuple({r['seq_region_id'] for r in regions} | {m['source_id'] for m in mappings})
        query = text('SELECT stable_id, version, gene_id FROM gene WHERE is_current=1 AND seq_region_id IN :ids').bindparams(bindparam('ids', expanding=True))
        for stable, version, db_id in connection.execute(query, dict(ids=ids)):
            identifier = stable if args.no_version or version in (None, '') else '{}.{}'.format(stable, version)
            if identifier in gene_order and gene_order[identifier] != db_id:
                raise ValueError('Ambiguous current gene identity in PAR annotation context')
            gene_order[identifier] = db_id
    host_name = 'chrY' if args.add_chr else 'Y'
    groups = defaultdict(list)
    for feature in features:
        groups[(feature.chrom, feature.attrs['gene_id'])].append(feature)
    with tempfile.TemporaryDirectory(prefix='gxf-par-') as directory:
        for index, mapping in enumerate(mappings):
            destination = Path(directory) / ('source_' + str(index) + '.' + args.format)
            command = ['--db-secrets', args.db_secrets, '--chromosome', mapping['source_name'],
                       '--start', str(mapping['source_start']), '--end', str(mapping['source_end']),
                       '--format', args.format, '--output', str(destination)]
            for name in ('vega', 'artifacts', 'no_version', 'include_hidden_remarks'):
                if getattr(args, name): command.append('--' + name.replace('_', '-'))
            for name in ('source', 'analysis', 'gene_list', 'transcript_list', 'mane'):
                value = getattr(args, name)
                if value: command += ['--' + name.replace('_', '-'), value]
            for biotype in args.biotype or []:
                command += ['--biotype', biotype]
            export(command)
            source_groups = defaultdict(list)
            for feature in read_features(destination, args.format):
                source_groups[feature.attrs['gene_id']].append(feature)
            for gene_id, records in source_groups.items():
                gene = next(f for f in records if f.feature_type == 'gene')
                if not mapping['source_start'] <= gene.start <= gene.end <= mapping['source_end']:
                    # Ensembl's PAR transfer requires an intact gene projection.
                    # XG/TMLHE/SPRY3 in the supplied legacy exports demonstrate
                    # that overlapping the interval alone does not create a Y copy.
                    print('Excluded non-projectable Y PAR gene {}: crosses the projection boundary'.format(gene_id), file=sys.stderr)
                    continue
                projected = []
                for feature in records:
                    if mapping['ori'] == 1:
                        start = mapping['host_start'] + feature.start - mapping['source_start']
                        end = mapping['host_start'] + feature.end - mapping['source_start']
                        strand = feature.strand
                    else:
                        start = mapping['host_end'] - (feature.end - mapping['source_start'])
                        end = mapping['host_end'] - (feature.start - mapping['source_start'])
                        strand = {'+':'-', '-':'+'}.get(feature.strand, feature.strand)
                    projected.append(replace(feature, chrom=host_name, start=start, end=end, strand=strand))
                projected_gene = next(f for f in projected if f.feature_type == 'gene')
                if args.start is not None and projected_gene.end < args.start: continue
                if args.end is not None and projected_gene.start > args.end: continue
                key = (host_name, gene_id)
                if key in groups:
                    existing = next(f for f in groups[key] if f.feature_type == 'gene')
                    if (existing.start, existing.end) != (projected_gene.start, projected_gene.end):
                        raise ValueError('Conflicting stored and projected Y PAR gene')
                else:
                    groups[key] = projected
    result = []
    for (chromosome, gene_id), records in sorted(groups.items(), key=lambda item: (item[0][0], gene_order[item[0][1]])):
        gene = next(f for f in records if f.feature_type == 'gene')
        if chromosome == host_name and any(m['host_start'] <= gene.start <= gene.end <= m['host_end'] for m in mappings):
            records = _par_ids(records, args.format)
        result.extend(records)
    return result
