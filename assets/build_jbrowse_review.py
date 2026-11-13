#!/usr/bin/env python3
"""Build a JBrowse Desktop review bundle from LEAP annotations and datachecks."""
import argparse
from collections import Counter
from contextlib import ExitStack
import csv
import gzip
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from urllib.parse import quote, unquote

ASSEMBLY_NAME = 'GRCh38'
REFERENCE = 'https://jbrowse.org/genomes/GRCh38/fasta/hg38.prefix.fa.gz'
LENGTHS = dict(
    zip(
        [str(i) for i in range(1, 23)] + ['X', 'Y'],
        [
            248956422,
            242193529,
            198295559,
            190214555,
            181538259,
            170805979,
            159345973,
            145138636,
            138394717,
            133797422,
            135086622,
            133275309,
            114364328,
            107043718,
            101991189,
            90338345,
            83257441,
            80373285,
            58617616,
            64444167,
            46709983,
            50818468,
            156040895,
            57227415,
        ],
    )
)
HUMAN_LENGTHS = dict(LENGTHS)
# UCSC mm39 primary chromosome sizes; annotation records keep their original unprefixed names.
MOUSE_LENGTHS = dict(
    zip(
        [str(i) for i in range(1, 20)] + ['X', 'Y'],
        [
            195154279,
            181755017,
            159745316,
            156860686,
            151758149,
            149588044,
            144995196,
            130127694,
            124359700,
            130530862,
            121973369,
            120092757,
            120883175,
            125139656,
            104073951,
            98008968,
            95294699,
            90720763,
            61420004,
            169476592,
            91455967,
        ],
    )
)


def configure_assembly(name):
    global ASSEMBLY_NAME, LENGTHS
    if name not in ('GRCh38', 'GRCm39'):
        raise ValueError('Unsupported reference assembly: ' + name)
    ASSEMBLY_NAME = name
    LENGTHS = dict(HUMAN_LENGTHS if name == 'GRCh38' else MOUSE_LENGTHS)


TRACKS = [
    ('raw', 'Raw annotation — all transcripts', '#476a96'),
    ('extended', 'Extended annotation — both strands', '#7257a3'),
    ('raw_cds', 'Raw CDS — coding segments', '#494949'),
    ('raw_utr', 'Raw UTR — original UTR segments', '#7897b8'),
    ('same_gene_failures', 'FAIL — same-gene CDS/exon overlap', '#c0392b'),
    ('cross_gene_failures', 'FAIL — cross-gene LEAP rule violation', '#c47415'),
    ('shared_site_allowed', 'PASS — shared internal splice site', '#27804a'),
    ('annotation_allowed', 'PASS — current-annotation extension', '#207e91'),
    ('grouped_annotation_allowed', 'PASS — existing grouped annotation', '#368c9b'),
    ('leap_allowed', 'AUDIT — exon contacts allowed by LEAP', '#537ca6'),
]
REPORT_COLUMNS = {
    'source_transcript',
    'source_gene',
    'target_transcript',
    'target_gene',
    'chromosome',
    'strand',
    'target_feature',
    'status',
    'failing_intervals',
    'allowed_intervals',
    'shared_splice_sites',
}


def open_text(path):
    return (
        gzip.open(path, 'rt', encoding='utf-8')
        if str(path).endswith('.gz')
        else open(path, encoding='utf-8')
    )


def discover(path, suffixes):
    path = path.resolve()
    if path.is_file():
        return [path]
    if path.is_dir():
        files = sorted(
            p
            for p in path.iterdir()
            if p.is_file() and any(p.name.endswith(s) for s in suffixes)
        )
        if files:
            return files
    raise ValueError('No matching input files: ' + str(path))


def coordinates(chromosome, start, end, strand, context):
    try:
        start, end = int(start), int(end)
    except (TypeError, ValueError):
        raise ValueError(context + ': coordinates must be integers') from None
    if chromosome not in LENGTHS:
        raise ValueError(
            context
            + ': expected '
            + ASSEMBLY_NAME
            + ' primary chromosome names; got '
            + str(chromosome)
        )
    if not 1 <= start <= end <= LENGTHS[chromosome] or strand not in (
        '+',
        '-',
        '.',
        '?',
    ):
        raise ValueError(context + ': invalid coordinates or strand')
    return start, end


def records(path):
    with open_text(path) as handle:
        for number, text in enumerate(handle, 1):
            if text.strip() == '##FASTA':
                break
            if not text.strip() or text.startswith('#'):
                continue
            fields = text.rstrip('\r\n').split('\t')
            context = '{}:{}'.format(path, number)
            if len(fields) != 9:
                raise ValueError(context + ': expected 9 GFF3 columns')
            coordinates(fields[0], fields[3], fields[4], fields[6], context)
            yield fields, context


def attributes(values):
    return ';'.join(
        key + '=' + quote(str(value), safe=':._-')
        for key, value in values.items()
        if value != ''
    )


def gff(chromosome, kind, start, end, strand, values, phase='.'):
    return (
        '\t'.join(
            map(
                str,
                [
                    chromosome,
                    'LEAP_review',
                    kind,
                    start,
                    end,
                    '.',
                    strand,
                    phase,
                    attributes(values),
                ],
            )
        )
        + '\n'
    )


def annotation_files(paths, label, outputs, counts, first_locus):
    for path in paths:
        for fields, context in records(path):
            kind = fields[2]
            if kind not in (
                'transcript',
                'exon',
                'CDS',
                'five_prime_UTR',
                'three_prime_UTR',
            ):
                continue
            attrs = {}
            for item in fields[8].split(';'):
                if not item:
                    continue
                key, sep, value = item.partition('=')
                if not sep or not key or key in attrs:
                    raise ValueError(context + ': malformed or duplicate attribute')
                attrs[key] = unquote(value)
            identifier = attrs.get('transcript_id')
            if not identifier or identifier == '.':
                raise ValueError(context + ': missing transcript_id')
            counts[label] += 1
            first_locus.setdefault(
                label,
                dict(chromosome=fields[0], start=int(fields[3]), end=int(fields[4])),
            )
            values = {
                key: attrs[key]
                for key in (
                    'transcript_id',
                    'gene_id',
                    'gene_name',
                    'exon_id',
                    'exon_number',
                    'transcript_type',
                    'leap_tag',
                    'transcript_class',
                )
                if key in attrs
            }
            if kind == 'transcript':
                values.update(
                    ID=identifier,
                    Name=(attrs.get('gene_name', '') + ' ' + identifier).strip(),
                )
            else:
                values.update(
                    ID='{}_{}'.format(label, counts[label]), Parent=identifier
                )
            outputs[label].write(
                gff(fields[0], kind, fields[3], fields[4], fields[6], values, fields[7])
            )
            feature_track = (
                'raw_cds'
                if kind == 'CDS'
                else 'raw_utr' if kind.endswith('_UTR') else None
            )
            if label == 'raw' and feature_track:
                values.pop('Parent', None)
                values['Name'] = identifier
                outputs[feature_track].write(
                    gff(
                        fields[0],
                        kind,
                        fields[3],
                        fields[4],
                        fields[6],
                        values,
                        fields[7],
                    )
                )
                counts[feature_track] += 1
        print('Prepared ' + str(path), flush=True)
    if not counts[label]:
        raise ValueError('No displayable annotation features in ' + label)


def report_files(paths, outputs, counts, examples, loci):
    serial = 0
    for path in paths:
        with path.open(newline='', encoding='utf-8') as handle:
            reader = csv.DictReader(handle)
            missing = REPORT_COLUMNS - set(reader.fieldnames or [])
            if missing:
                raise ValueError(
                    '{}: missing report columns: {}'.format(
                        path, ', '.join(sorted(missing))
                    )
                )
            for row in reader:
                context = '{}:{}'.format(path, reader.line_num)
                if any(row.get(key) is None for key in REPORT_COLUMNS):
                    raise ValueError(context + ': incomplete report row')
                if row['status'] == 'missing_from_raw':
                    counts['missing_from_raw'] += 1
                    continue
                if row['status'] not in (
                    'new_cds_overlap',
                    'new_within_gene_exon_overlap',
                    'new_cross_gene_exon_overlap',
                    'allowed_shared_splice_site',
                    'allowed_current_annotation',
                    'allowed_grouped_annotation',
                    'allowed_multiple_reasons',
                    'allowed_leap_overlap',
                    'new_cross_gene_cds_overlap',
                    'unsafe_terminal_overlap',
                ):
                    raise ValueError(context + ': unrecognized overlap status')
                if row['target_feature'] not in ('CDS', 'exon'):
                    raise ValueError(context + ': target_feature must be CDS or exon')
                if not all(
                    row[key]
                    for key in (
                        'source_transcript',
                        'target_transcript',
                        'source_gene',
                        'target_gene',
                    )
                ):
                    raise ValueError(context + ': missing transcript or gene ID')
                all_intervals = []
                marker_columns = [
                    (
                        'failing_intervals',
                        (
                            'same_gene_failures'
                            if row['source_gene'] == row['target_gene']
                            else 'cross_gene_failures'
                        ),
                    )
                ]
                if (
                    'annotation_allowed_intervals' in row
                    and 'splice_site_allowed_intervals' in row
                ):
                    marker_columns += [
                        ('annotation_allowed_intervals', 'annotation_allowed'),
                        ('splice_site_allowed_intervals', 'shared_site_allowed'),
                    ]
                elif row['status'] in (
                    'allowed_current_annotation',
                    'allowed_grouped_annotation',
                    'allowed_multiple_reasons',
                    'allowed_leap_overlap',
                    'new_cross_gene_cds_overlap',
                    'unsafe_terminal_overlap',
                ):
                    raise ValueError(
                        context
                        + ': provenance-aware report requires reason-specific interval columns'
                    )
                else:
                    # Reports from before the provenance exception only allowed shared sites.
                    marker_columns.append(('allowed_intervals', 'shared_site_allowed'))
                if row['status'] == 'allowed_grouped_annotation' and 'grouped_annotation_allowed_intervals' not in row:
                    raise ValueError(context + ': grouped-annotation report requires its interval column')
                if 'grouped_annotation_allowed_intervals' in row:
                    marker_columns.append(('grouped_annotation_allowed_intervals', 'grouped_annotation_allowed'))
                if 'leap_allowed_intervals' in row:
                    marker_columns.append(('leap_allowed_intervals', 'leap_allowed'))
                for column, name in marker_columns:
                    if row.get(column) is None:
                        raise ValueError(context + ': incomplete report row')
                    for interval in filter(None, row[column].split(';')):
                        if interval.count('-') != 1:
                            raise ValueError(
                                context + ': expected start-end overlap interval'
                            )
                        start, end = coordinates(
                            row['chromosome'],
                            *interval.split('-'),
                            row['strand'],
                            context
                        )
                        serial += 1
                        values = {
                            key: row[key]
                            for key in (
                                'source_transcript',
                                'target_transcript',
                                'source_gene',
                                'target_gene',
                                'shared_splice_sites',
                            )
                        }
                        values['annotation_extension_ends'] = row.get(
                            'annotation_extension_ends', ''
                        )
                        for key in (
                            'check_rule',
                            'gene_relationship',
                            'target_is_partial',
                            'target_missing_utr',
                            'intentional_partial_cds_check',
                            'grouped_annotation_boundaries',
                            'grouped_annotation_donors',
                        ):
                            values[key] = row.get(key, '')
                        values.update(
                            ID='overlap_' + str(serial),
                            Name=row['source_transcript']
                            + ' vs '
                            + row['target_transcript'],
                            comparison_status=row['status'],
                            interval_outcome=name,
                        )
                        outputs[name].write(
                            gff(
                                row['chromosome'],
                                'sequence_feature',
                                start,
                                end,
                                row['strand'],
                                values,
                            )
                        )
                        counts[name] += 1
                        examples.setdefault(name, dict(row, start=start, end=end))
                        all_intervals.append((start, end))
                if not all_intervals:
                    raise ValueError(context + ': overlap row contains no intervals')
                if all_intervals:
                    start = max(1, min(x[0] for x in all_intervals) - 2000)
                    end = min(
                        LENGTHS[row['chromosome']],
                        max(x[1] for x in all_intervals) + 2000,
                    )
                    loci.writerow(
                        ['{}:{}..{}'.format(row['chromosome'], start, end)]
                        + [
                            row[k]
                            for k in (
                                'status',
                                'source_transcript',
                                'target_transcript',
                                'source_gene',
                                'target_gene',
                                'shared_splice_sites',
                            )
                        ]
                    )
        print('Prepared markers ' + str(path), flush=True)


def index_track(data, name, threads):
    source, output = data / (name + '.unsorted.gff3'), data / (name + '.gff3.gz')
    with output.open('wb') as handle:
        sort = subprocess.Popen(
            [
                'sort',
                '-S',
                '256M',
                '-T',
                str(data),
                '-k1,1',
                '-k4,4n',
                '-k5,5n',
                str(source),
            ],
            stdout=subprocess.PIPE,
            env=dict(os.environ, LC_ALL='C'),
        )
        try:
            compress = subprocess.Popen(
                ['bgzip', '-@', str(threads), '-c'], stdin=sort.stdout, stdout=handle
            )
        except Exception:
            sort.stdout.close()
            sort.terminate()
            sort.wait()
            raise
        sort.stdout.close()
        compressed, sorted_ok = compress.wait(), sort.wait()
        if compressed or sorted_ok:
            raise RuntimeError('sort/bgzip failed for ' + name)
    subprocess.run(['tabix', '-f', '-p', 'gff', str(output)], check=True)
    source.unlink()
    print('Indexed ' + name, flush=True)


def local(path):
    return dict(locationType='LocalPathLocation', localPath=str(path))


def session(name, locus, names):
    chromosome = locus['chromosome']
    start, end = max(0, locus['start'] - 1 - 2000), min(
        LENGTHS[chromosome], locus['end'] + 2000
    )
    scale = (end - start) / 1000
    return dict(
        id='leap-' + names[0] + '-' + name.split(' — ')[0].lower().replace(' ', '-'),
        name=name,
        focusedViewId='leap-review-view',
        views=[
            dict(
                id='leap-review-view',
                type='LinearGenomeView',
                displayedRegions=[
                    dict(
                        refName=chromosome,
                        start=0,
                        end=LENGTHS[chromosome],
                        assemblyName=ASSEMBLY_NAME,
                    )
                ],
                bpPerPx=scale,
                offsetPx=start / scale,
                tracks=[
                    dict(
                        type='FeatureTrack',
                        configuration=track,
                        displays=[
                            dict(
                                type='LinearBasicDisplay',
                                configuration=track + '-LinearBasicDisplay',
                                height=220 if track in ('raw', 'extended') else 95,
                            )
                        ],
                    )
                    for track in names
                ],
            )
        ],
    )


def config_document(root, track_specs, examples, first_locus):
    tracks = []
    for name, label, color in track_specs:
        path = root / 'data' / (name + '.gff3.gz')
        tracks.append(
            dict(
                type='FeatureTrack',
                trackId=name,
                name=label,
                assemblyNames=[ASSEMBLY_NAME],
                category=['LEAP review'],
                adapter=dict(
                    type='Gff3TabixAdapter',
                    gffGzLocation=local(path),
                    index=dict(
                        location=local(Path(str(path) + '.tbi')), indexType='TBI'
                    ),
                ),
                displays=[
                    dict(
                        type='LinearBasicDisplay',
                        displayId=name + '-LinearBasicDisplay',
                        renderer=dict(type='SvgFeatureRenderer', color1=color),
                    )
                ],
            )
        )
    uri = lambda value: dict(uri=value, locationType='UriLocation')
    assembly = dict(
        name='GRCh38',
        aliases=['hg38'],
        sequence=dict(
            type='ReferenceSequenceTrack',
            trackId='GRCh38-ReferenceSequenceTrack',
            adapter=dict(
                type='BgzipFastaAdapter',
                fastaLocation=uri(REFERENCE),
                faiLocation=uri(REFERENCE + '.fai'),
                gziLocation=uri(REFERENCE + '.gzi'),
            ),
        ),
    )
    if ASSEMBLY_NAME == 'GRCm39':
        assembly = dict(
            name='GRCm39',
            aliases=['mm39'],
            sequence=dict(
                type='ReferenceSequenceTrack',
                trackId='GRCm39-ReferenceSequenceTrack',
                adapter=dict(
                    type='TwoBitAdapter',
                    twoBitLocation=uri(
                        'https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.2bit'
                    ),
                    chromSizesLocation=uri(
                        'https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.chrom.sizes'
                    ),
                ),
            ),
            refNameAliases=dict(
                adapter=dict(
                    type='RefNameAliasAdapter',
                    location=local(root / 'data/refname-aliases.tsv'),
                )
            ),
        )
    names = [spec[0] for spec in track_specs]
    sessions = []
    for key, label, _ in track_specs:
        if key in examples:
            sessions.append(
                session(
                    label + ' — ' + examples[key].get('source_transcript', 'example'),
                    examples[key],
                    names,
                )
            )
    if not sessions:
        sessions.append(
            session(
                'LEAP annotation review — no overlap markers',
                first_locus['extended'],
                names,
            )
        )
    return dict(
        assemblies=[assembly],
        tracks=tracks,
        defaultSession=sessions[0],
        preConfiguredSessions=sessions,
    )


def build(args):
    configure_assembly(getattr(args, 'assembly', 'GRCh38'))
    for tool in ('sort', 'bgzip', 'tabix'):
        if not shutil.which(tool):
            raise ValueError('Required command not found: ' + tool)
    if args.threads < 1:
        raise ValueError('--threads must be at least 1')
    root = args.output.resolve()
    if root.exists():
        raise ValueError('Output already exists; choose a new directory: ' + str(root))
    suffixes = ('.gff', '.gff3', '.gff.gz', '.gff3.gz')
    raw, extended = discover(args.raw_annotation, suffixes), discover(
        args.extended_annotation, suffixes
    )
    reports = discover(args.reports, ('.csv',)) if args.reports else []
    extras = []
    for index, argument in enumerate(args.extra_track, 1):
        label, separator, path = argument.partition('=')
        if not separator or not label or not path:
            raise ValueError('--extra-track must be NAME=PATH')
        extras.append(('extra_' + str(index), label, discover(Path(path), suffixes)))
    root.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix='.' + root.name + '-build-', dir=str(root.parent)
    ) as temporary:
        staging = Path(temporary) / 'bundle'
        data = staging / 'data'
        data.mkdir(parents=True)
        specs = (
            TRACKS[:4]
            + (TRACKS[4:] if reports else [])
            + [(key, label, '#397d83') for key, label, _ in extras]
        )
        counts, examples, first_locus = Counter(), {}, {}
        with ExitStack() as stack:
            outputs = {
                name: stack.enter_context(
                    (data / (name + '.unsorted.gff3')).open('w', encoding='utf-8')
                )
                for name, _, _ in specs
            }
            for handle in outputs.values():
                handle.write('##gff-version 3\n')
            annotation_files(raw, 'raw', outputs, counts, first_locus)
            annotation_files(extended, 'extended', outputs, counts, first_locus)
            loci_file = stack.enter_context(
                (staging / 'review-loci.tsv').open('w', newline='', encoding='utf-8')
            )
            loci = csv.writer(loci_file, delimiter='\t')
            loci.writerow(
                [
                    'location',
                    'status',
                    'source_transcript',
                    'target_transcript',
                    'source_gene',
                    'target_gene',
                    'shared_splice_sites',
                ]
            )
            report_files(reports, outputs, counts, examples, loci)
            for key, label, paths in extras:
                for path in paths:
                    for fields, context in records(path):
                        outputs[key].write('\t'.join(fields) + '\n')
                        counts[key] += 1
                        examples.setdefault(
                            key,
                            dict(
                                chromosome=fields[0],
                                start=int(fields[3]),
                                end=int(fields[4]),
                            ),
                        )
        for name, _, _ in specs:
            index_track(data, name, args.threads)
        if ASSEMBLY_NAME == 'GRCm39':
            (data / 'refname-aliases.tsv').write_text(
                ''.join('chr' + chrom + '\t' + chrom + '\n' for chrom in LENGTHS),
                encoding='utf-8',
            )
        # Config locations refer to the final directory, never the staging directory.
        config = config_document(root, specs, examples, first_locus)
        for filename in ('config.json', 'LEAP-overlap-review.jbrowse'):
            (staging / filename).write_text(
                json.dumps(config, indent=2) + '\n', encoding='utf-8'
            )
        summary = dict(
            assembly=ASSEMBLY_NAME,
            raw_annotation=[str(p) for p in raw],
            extended_annotation=[str(p) for p in extended],
            reports=[str(p) for p in reports],
            extra_tracks=[
                dict(name=label, files=[str(p) for p in paths])
                for _, label, paths in extras
            ],
            feature_counts=dict(counts),
            marker_count_unit='individual interval features, not unique transcripts',
        )
        (staging / 'build-summary.json').write_text(
            json.dumps(summary, indent=2) + '\n', encoding='utf-8'
        )
        (staging / 'examples.json').write_text(
            json.dumps(examples, indent=2) + '\n', encoding='utf-8'
        )
        (staging / 'README.txt').write_text(
            'Open LEAP-overlap-review.jbrowse using Open saved session in JBrowse Desktop.\n'
            'Raw/extended annotations are followed by raw CDS/UTR and datacheck tracks.\n'
            'Red: same-gene failure. Orange: cross-gene failure. Green: shared-site allowance. Teal: current-annotation allowance. Blue-green: existing grouped annotation. Blue: contacts allowed by LEAP.\n'
            'Click a marker for source/target IDs, splice evidence and annotation-supported ends. Mixed comparisons may have allowed and failing markers.\n'
            'File > Pre-configured sessions offers example locations. review-loci.tsv lists overlap results.\n'
            'Markers follow rule-specific report intervals. CDS clipping includes the original outer boundary; terminal tests require positive-length BED contact.\n'
            'All intervals use inclusive GFF3 coordinates. Markers follow the input reports without re-evaluating biology.\n'
            'Keep this data directory in place: config paths are absolute and desktop-specific.\n'
            'Internet is required for the public reference sequence. Annotation and marker tracks are local.\n'
            'See build-summary.json for input paths and feature counts. Rebuild after regenerating annotations AND reports.\n',
            encoding='utf-8',
        )
        if root.exists():
            raise ValueError(
                'Output appeared during build; refusing to overwrite: ' + str(root)
            )
        staging.rename(root)
    print('READY: ' + str(root / 'LEAP-overlap-review.jbrowse'))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--assembly',
        choices=('GRCh38', 'GRCm39'),
        default='GRCh38',
        help='Reference assembly; default GRCh38',
    )
    parser.add_argument(
        '--raw-annotation', type=Path, required=True, help='Raw GFF3 file or directory'
    )
    parser.add_argument(
        '--extended-annotation',
        type=Path,
        required=True,
        help='Extended GFF3 file or directory',
    )
    parser.add_argument(
        '--reports', type=Path, help='CDS/UTR overlap CSV file or directory; optional'
    )
    parser.add_argument(
        '--output',
        type=Path,
        required=True,
        help='New review bundle directory; must not exist',
    )
    parser.add_argument(
        '--threads', type=int, default=2, help='bgzip threads (default: 2)'
    )
    parser.add_argument(
        '--extra-track',
        action='append',
        default=[],
        metavar='NAME=PATH',
        help='Additional datacheck GFF3 file or directory; repeatable',
    )
    args = parser.parse_args()
    try:
        build(args)
        return 0
    except (
        OSError,
        EOFError,
        UnicodeError,
        ValueError,
        RuntimeError,
        csv.Error,
        subprocess.SubprocessError,
    ) as error:
        print('ERROR: ' + str(error), file=sys.stderr)
        return 2


if __name__ == '__main__':
    sys.exit(main())
