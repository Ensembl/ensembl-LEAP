#!/usr/bin/env python3
"""Prepare annotation tracks and an editable HAVANA JBrowser launch bundle."""

import argparse
from collections import Counter
import csv
import gzip
import importlib.util
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import tempfile

sys.dont_write_bytecode = True

GROUPS = ('raw', 'leaped', 'external', 'internal')
STRANDS = {'positive': '+', 'negative': '-'}
SHEET_FIELDS = ['sample', 'type', 'file', 'reference_file', 'bin_size', 'normalization',
                'region_format', 'chrom_sizes', 'genome_id', 'track_category',
                'ref_name_aliases', 'species', 'data_category', 'ftp_project', 'publication']
DESCRIPTIONS = {
    'raw': 'Published raw annotation before LEAP extension.',
    'leaped': 'Original LEAP annotation, including any internal extensions.',
    'external': 'Processed annotation: internal terminal ends restored to raw; retained external extensions preserved. Models without CDS omitted.',
    'internal': 'Original LEAP models with actual growth at an internal terminal end: the structures that would be present if internal extensions were retained. May also retain an external extension at the other end.',
}


def processor():
    candidates = [Path(__file__).with_name('process_internal_extensions.py'),
                  Path(__file__).resolve().parents[1] / 'assets/process_internal_extensions.py']
    path = next((p for p in candidates if p.is_file()), None)
    if path is None:
        raise ValueError('Missing process_internal_extensions.py helper')
    spec = importlib.util.spec_from_file_location('review_internal_processor', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def read_processor():
    path = Path(__file__).with_name('prepare_review_reads.py')
    if not path.is_file():
        raise ValueError('Missing prepare_review_reads.py helper')
    spec = importlib.util.spec_from_file_location('review_alignment_processor', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def records(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', encoding='utf-8') as handle:
        for number, line in enumerate(handle, 1):
            if line.strip() == '##FASTA':
                break
            if not line.strip() or line.startswith('#'):
                continue
            fields = line.rstrip('\r\n').split('\t')
            context = '{}:{}'.format(path, number)
            if len(fields) != 9:
                raise ValueError(context + ': expected nine GFF3 columns')
            try:
                start, end = int(fields[3]), int(fields[4])
            except ValueError:
                raise ValueError(context + ': noninteger coordinates') from None
            if not 1 <= start <= end or fields[6] not in ('+', '-', '.', '?'):
                raise ValueError(context + ': invalid interval or strand')
            yield fields


def write_csv(path, fields, rows, delimiter=','):
    with Path(path).open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter=delimiter, lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def new_directory(path):
    path = Path(path)
    if path.exists() or path.is_symlink():
        raise ValueError('Output must be a new directory: ' + str(path))
    path.mkdir(parents=True)
    return path


def validate_scope(chromosome, strand=None):
    if not re.fullmatch(r'[A-Za-z0-9_.-]+', chromosome) or chromosome in ('.', '..'):
        raise ValueError('Invalid chromosome scope')
    if strand is not None and strand not in STRANDS:
        raise ValueError('Invalid strand scope')


def raw_partition(chromosome, annotation, output):
    validate_scope(chromosome)
    output = new_directory(output)
    counts = Counter()
    handles = {s: (output / (s + '.gff3')).open('w') for s in ('positive', 'negative', 'unsplit')}
    try:
        for handle in handles.values():
            handle.write('##gff-version 3\n')
        for fields in records(annotation):
            if fields[0] != chromosome:
                raise ValueError('Raw annotation contains records outside chromosome ' + chromosome)
            line = '\t'.join(fields) + '\n'
            handles['unsplit'].write(line)
            counts['unsplit'] += 1
            for strand, sign in STRANDS.items():
                if fields[6] == sign:
                    handles[strand].write(line)
                    counts[strand] += 1
    finally:
        for handle in handles.values():
            handle.close()
    (output / 'partition.json').write_text(json.dumps(dict(kind='raw', chromosome=chromosome,
        source=Path(annotation).name, counts=dict(counts)), indent=2) + '\n')


def extension_partition(chromosome, strand, raw, extended, stats, output):
    validate_scope(chromosome, strand)
    output = new_directory(output)
    module = processor()
    name = chromosome + '_' + strand + '_extended'
    # Keep the pipeline basename, even when the standalone caller supplied gzip.
    local_extended = output / (name + '.gff')
    count = 0
    with local_extended.open('w') as handle:
        handle.write('##gff-version 3\n')
        for fields in records(extended):
            if fields[0] != chromosome or fields[6] != STRANDS[strand]:
                raise ValueError('Extended annotation contains records outside ' + chromosome + ' ' + strand)
            handle.write('\t'.join(fields) + '\n')
            count += 1
    folder = output / 'processed'
    for sub in ('gff', 'stats', 'ids', 'potential_internal'):
        (folder / sub).mkdir(parents=True)
    if count:
        result = module.process_file(local_extended, Path(raw), Path(stats), folder)
    else:
        if module.load_stats(stats):
            raise ValueError('Empty annotation has nonempty extension stats')
        for sub in ('gff', 'potential_internal'):
            (folder / sub / (name + '.gff')).write_text('##gff-version 3\n')
        write_csv(folder / 'stats' / (chromosome + '_' + strand + '_stats.tsv'), module.STATS_FIELDS, [], '\t')
        write_csv(folder / 'ids' / (name + '_classification.tsv'),
                  ['chromosome', 'transcript_id', 'classification'], [], '\t')
        result = dict(counts={})
    classification = {}
    with (folder / 'ids' / (name + '_classification.tsv')).open() as handle:
        classification = {row['transcript_id']: row['classification'] for row in csv.DictReader(handle, delimiter='\t')}
    original_stats = module.load_stats(stats)
    internal_ids = set()
    for identifier, category in classification.items():
        if not category.startswith('internal'):
            continue
        for side, prime in (('left', 'five_prime' if strand == 'positive' else 'three_prime'),
                            ('right', 'three_prime' if strand == 'positive' else 'five_prime')):
            if 'internal_' + side in category.split(',') and int(original_stats[identifier][prime + '_extension']) > 0:
                internal_ids.add(identifier)
    shutil.copyfile(local_extended, output / 'leaped.gff3')
    shutil.copyfile(folder / 'gff' / (name + '.gff'), output / 'external.gff3')
    with (output / 'internal.gff3').open('w') as handle:
        handle.write('##gff-version 3\n')
        for fields in records(folder / 'potential_internal' / (name + '.gff')):
            attrs = dict(item.split('=', 1) for item in fields[8].split(';') if item)
            if attrs.get('transcript_id') not in internal_ids:
                continue
            if fields[2] == 'transcript':
                attrs['transcript_class'] = classification[attrs['transcript_id']]
                fields[8] = ';'.join(k + '=' + v for k, v in attrs.items())
            handle.write('\t'.join(fields) + '\n')
    local_extended.unlink()
    (output / 'partition.json').write_text(json.dumps(dict(kind='extended', chromosome=chromosome,
        strand=strand, source=Path(extended).name, stats_source=Path(stats).name,
        counts=dict(result['counts'], internal_extension_models=len(internal_ids))), indent=2) + '\n')


def sorted_annotation(paths, destination):
    counts = Counter()
    with tempfile.NamedTemporaryFile(mode='w', dir=destination.parent, delete=False) as temporary:
        source = Path(temporary.name)
        for path in paths:
            for fields in records(path):
                temporary.write('\t'.join(fields) + '\n')
                counts['features'] += 1
                counts['transcripts'] += fields[2] == 'transcript'
    try:
        with destination.open('w') as handle:
            handle.write('##gff-version 3\n')
            subprocess.run(['sort', '-s', '-t', '\t', '-k1,1', '-k4,4n', '-k5,5n', str(source)],
                stdout=handle, env=dict(os.environ, LC_ALL='C'), check=True)
    finally:
        source.unlink()
    return dict(features=counts['features'], transcripts=counts['transcripts'])


def launch_script(settings, offline=False):
    # Values are shell-quoted data; no database values or guessed cluster paths.
    values = dict(JBROWSER_PIPELINE=settings.get('jbrowser_pipeline') or 'EDIT_ME/havana-jbrowser',
                  PROFILE=settings.get('profile') or 'slurm,singularity')
    if not offline:
        values.update(PUBLIC_BASE_URL=settings.get('public_base_url', ''),
                      FTP_TARGET_DIR=settings.get('ftp_target_dir', ''), HUB_EMAIL=settings.get('hub_email', ''))
    assignments = '\n'.join(k + '=' + shlex.quote(str(v)) for k, v in values.items())
    sheet = 'samplesheet_offline.csv' if offline else 'samplesheet.csv'
    suffix = '_offline' if offline else ''
    prepare_sheet = '''# Refresh the offline sheet from the shared settings, forcing every row offline.
python3 - <<'PY_OFFLINE'
import csv
import os
import tempfile
from pathlib import Path
with open('samplesheet.csv', newline='') as handle:
    reader=csv.DictReader(handle)
    fields=list(reader.fieldnames or [])
    rows=list(reader)
if not fields:
    raise SystemExit('samplesheet.csv must contain a CSV header.')
if 'publication' not in fields:
    fields.append('publication')
with tempfile.NamedTemporaryFile(mode='w', newline='', dir='.', delete=False) as handle:
    temporary=Path(handle.name)
    writer=csv.DictWriter(handle, fieldnames=fields)
    writer.writeheader()
    writer.writerows(dict(row, publication='offline') for row in rows)
os.replace(temporary, 'samplesheet_offline.csv')
PY_OFFLINE
''' if offline else ''
    return '''#!/usr/bin/env bash
#SBATCH --job-name=LEAP_manual_review
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=48:00:00
# Add your site's --partition and --account settings if needed.
set -euo pipefail

''' + assignments + '''

# Edit settings above; load your site's Nextflow and container modules as needed.
# REVIEW_DIR is recorded because Slurm spools scripts outside their source folder.
REVIEW_DIR=''' + shlex.quote(settings['published_root']) + '''
cd -- "$REVIEW_DIR"
[[ -f "$JBROWSER_PIPELINE/main.nf" ]] || { echo 'Set JBROWSER_PIPELINE to your server checkout.' >&2; exit 1; }
command -v nextflow >/dev/null || { echo 'Load Nextflow before submitting this script.' >&2; exit 1; }
''' + prepare_sheet + '''
python3 - <<'PY_VALIDATE'
import csv
from pathlib import Path
with open(''' + repr(sheet) + ''', newline='') as handle:
    rows=list(csv.DictReader(handle))
if not rows:
    raise SystemExit('No nonempty review tracks are available.')
for row in rows:
    if any('EDIT_ME' in value for value in row.values()):
        raise SystemExit('Complete reference/assembly/species settings in samplesheet.csv before launching.')
    if not Path(row['file']).is_file():
        raise SystemExit('Samplesheet input is missing: '+row['file'])
    if row['type'] == 'bam' and (not row['reference_file'] or not Path(row['reference_file']).is_file()):
        raise SystemExit('Alignment reference is missing: '+row['reference_file'])
PY_VALIDATE
command=(nextflow run "$JBROWSER_PIPELINE/main.nf"
    --input "$REVIEW_DIR/''' + sheet + '''"
    --sample_metadata "$REVIEW_DIR/sample_metadata.csv"
    --outdir "$REVIEW_DIR/jbrowser_output''' + suffix + '''"
    -profile "$PROFILE" -work-dir "$REVIEW_DIR/jbrowser_work''' + suffix + '''"
    -with-report "$REVIEW_DIR/jbrowser_report''' + suffix + '''.html"
    -with-trace "$REVIEW_DIR/jbrowser_trace''' + suffix + '''.txt")
''' + ('' if offline else '''
if python3 -c "import csv,sys; sys.exit(not any(r['publication']=='public' for r in csv.DictReader(open('samplesheet.csv'))))"; then
    [[ -n "$PUBLIC_BASE_URL" && -n "$FTP_TARGET_DIR" && -n "$HUB_EMAIL" ]] || {
        echo 'Set PUBLIC_BASE_URL, FTP_TARGET_DIR and HUB_EMAIL for public FTP preparation.' >&2; exit 1;
    }
    command+=(--public_base_url "$PUBLIC_BASE_URL" --ftp_target_dir "$FTP_TARGET_DIR" --hub_email "$HUB_EMAIL")
fi
''') + '''
"${command[@]}" "$@"
echo ''' + shlex.quote('Offline JBrowser preparation finished: jbrowser_output_offline.' if offline else
    'JBrowser preparation finished. Review its output and generated FTP copy scripts before publishing.') + '''
'''


def build_bundle(partitions, scopes, settings, output, reference_lengths=None, read_evidence=None):
    output = Path(output).absolute()
    if output.exists() or output.is_symlink():
        raise ValueError('Output must be a new directory: ' + str(output))
    settings = dict(settings)
    if settings.get('publication', 'public') not in ('public', 'offline'):
        raise ValueError('Review publication must be public or offline')
    reference_name = Path(settings.get('reference_genome') or '').name.lower()
    for names, genome, species in ((('grch38', 'hg38'), 'hg38', 'homo_sapiens'),
                                    (('grcm39', 'mm39'), 'mm39', 'mus_musculus')):
        if any(name in reference_name for name in names):
            if not settings.get('genome_id'):
                settings['genome_id'] = genome
            if not settings.get('species'):
                settings['species'] = species
    expected = {(str(c), s) for c, s in scopes}
    if not expected or len(expected) != len(scopes):
        raise ValueError('Expected review scopes must be nonempty and unique')
    for c, s in expected:
        validate_scope(c, s)
    raw, extended = {}, {}
    for folder in partitions:
        folder = Path(folder)
        data = json.loads((folder / 'partition.json').read_text())
        if data.get('kind') not in ('raw', 'extended'):
            raise ValueError('Unknown review partition type: ' + str(folder))
        key = data['chromosome'] if data['kind'] == 'raw' else (data['chromosome'], data['strand'])
        mapping = raw if data['kind'] == 'raw' else extended
        if key in mapping:
            raise ValueError('Duplicate review partition: ' + str(key))
        mapping[key] = (folder, data)
    if set(raw) != {c for c, s in expected} or set(extended) != expected:
        raise ValueError('Missing or unexpected raw/extended review partitions')
    read_summary = dict(status='not_supplied', description='Read evidence was not supplied to this standalone annotation-only build.')
    if read_evidence:
        read_evidence = Path(read_evidence)
        read_summary = json.loads((read_evidence / 'read-summary.json').read_text())
        counts = read_summary.get('counts', {})
        if ({tuple(scope) for scope in read_summary.get('scopes', [])} != expected or
                set(counts) != {'positive', 'negative', 'unsplit'} or
                any(type(n) is not int or n < 0 for n in counts.values()) or
                counts['positive'] + counts['negative'] != counts['unsplit']):
            raise ValueError('Read evidence has incompatible scopes or alignment counts')
        for view in ('positive', 'negative', 'unsplit'):
            for suffix in ('.bam', '.bam.csi'):
                if not (read_evidence / (view + suffix)).is_file():
                    raise ValueError('Missing prepared read-evidence BAM/index: ' + view + suffix)
        read_summary = dict(read_summary, status='available')
        if settings.get('read_input_bam'):
            read_summary['original_input'] = dict(bam=settings['read_input_bam'], index=settings['read_input_index'])
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.' + output.name + '-', dir=output.parent) as temporary:
        staging = Path(temporary) / 'bundle'
        staging.mkdir()
        published = Path(settings.get('published_root') or output).absolute()
        settings = dict(settings, published_root=str(published))
        manifest, samples, metadata = [], [], []
        chrom_sizes = ''
        if reference_lengths:
            spec = importlib.util.spec_from_file_location('review_reference_check',
                next(p for p in (Path(__file__).with_name('check_reference_coordinates.py'),
                    Path(__file__).resolve().parents[1] / 'assets/check_reference_coordinates.py') if p.is_file()))
            checker = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(checker)
            lengths = checker.reference_lengths(Path(reference_lengths))
            (staging / 'reference.chrom.sizes').write_text(''.join('{}\t{}\n'.format(k, v) for k, v in lengths.items()))
            chrom_sizes = str(published / 'reference.chrom.sizes')
        for group in GROUPS:
            folder = staging / 'annotations' / group
            folder.mkdir(parents=True)
            for strand in ('positive', 'negative', 'unsplit'):
                if group == 'raw':
                    paths = [raw[c][0] / (strand + '.gff3') for c in sorted(raw)]
                else:
                    paths = [extended[k][0] / (group + '.gff3') for k in sorted(extended)
                             if strand == 'unsplit' or k[1] == strand]
                target = folder / (strand + '.gff3')
                counts = sorted_annotation(paths, target)
                relative = target.relative_to(staging)
                included = counts['transcripts'] > 0
                manifest.append(dict(group=group, strand=strand, file=str(relative), **counts,
                    included_in_samplesheet=included))
                if included:
                    sample = 'LEAP_' + group + '_' + strand
                    samples.append(dict(sample=sample, type='gff3', file=str(published / relative),
                        chrom_sizes=chrom_sizes, genome_id=settings.get('genome_id') or 'EDIT_ME_ASSEMBLY',
                        species=settings.get('species') or 'EDIT_ME_SPECIES', track_category='LEAP ' + group,
                        data_category='annotation', ftp_project=settings.get('ftp_project') or 'LEAP_review',
                        publication=settings.get('publication', 'public')))
                    metadata.append(dict(sample=sample, tags='LEAP;manual_review;' + group + ';' + strand,
                        description=DESCRIPTIONS[group] + ' Strand view: ' + strand + '.'))
        if read_evidence:
            reads = staging / 'reads'
            reads.mkdir()
            shutil.copyfile(read_evidence / 'read-summary.json', reads / 'read-summary.json')
            for view in ('positive', 'negative', 'unsplit'):
                for suffix in ('.bam', '.bam.csi'):
                    shutil.copyfile(read_evidence / (view + suffix), reads / (view + suffix))
                count = read_summary['counts'][view]
                relative = 'reads/' + view + '.bam'
                manifest.append(dict(group='reads', strand=view, file=relative, track_type='bam',
                    alignment_records=count, index_file=relative + '.csi', included_in_samplesheet=count > 0))
                if count:
                    sample = 'LEAP_reads_' + view
                    samples.append(dict(sample=sample, type='bam', file=str(published / relative),
                        reference_file=settings.get('reference_genome') or 'EDIT_ME/reference.fa',
                        genome_id=settings.get('genome_id') or 'EDIT_ME_ASSEMBLY',
                        species=settings.get('species') or 'EDIT_ME_SPECIES', track_category='LEAP read evidence',
                        data_category='alignment', ftp_project=settings.get('ftp_project') or 'LEAP_review',
                        publication=settings.get('publication', 'public')))
                    metadata.append(dict(sample=sample, tags='LEAP;manual_review;read_evidence;' + view,
                        description=read_summary['description'] + ' Strand view: ' + view + '.'))
        if samples:
            reference = settings.get('reference_genome') or 'EDIT_ME/reference.fa'
            samples.insert(0, dict(sample='LEAP_reference', type='reference', file=reference,
                chrom_sizes=chrom_sizes, genome_id=settings.get('genome_id') or 'EDIT_ME_ASSEMBLY',
                species=settings.get('species') or 'EDIT_ME_SPECIES', data_category='reference',
                ftp_project=settings.get('ftp_project') or 'LEAP_review', publication=settings.get('publication', 'public')))
            metadata.insert(0, dict(sample='LEAP_reference', tags='LEAP;reference', description='Matching reference supplied for manual review.'))
        for key, (folder, data) in extended.items():
            target = staging / 'processing' / (key[0] + '_' + key[1])
            target.parent.mkdir(exist_ok=True)
            shutil.copytree(folder / 'processed' / 'stats', target / 'stats')
            shutil.copytree(folder / 'processed' / 'ids', target / 'ids')
        write_csv(staging / 'samplesheet.csv', SHEET_FIELDS, samples)
        write_csv(staging / 'samplesheet_offline.csv', SHEET_FIELDS,
                  [dict(row, publication='offline') for row in samples])
        write_csv(staging / 'sample_metadata.csv', ['sample', 'tags', 'description'], metadata)
        write_csv(staging / 'track_manifest.tsv', ['group', 'strand', 'file', 'features', 'transcripts',
            'included_in_samplesheet', 'track_type', 'alignment_records', 'index_file'], manifest, '\t')
        summary = dict(settings=settings, scopes=sorted(expected), tracks=manifest,
            partitions=[data for folder, data in extended.values()], raw_partitions=[data for folder, data in raw.values()],
            read_evidence=read_summary)
        (staging / 'review-summary.json').write_text(json.dumps(summary, indent=2) + '\n')
        script = staging / 'run_jbrowser.sbatch'
        script.write_text(launch_script(settings))
        script.chmod(0o755)
        offline_script = staging / 'run_jbrowser_offline.sbatch'
        offline_script.write_text(launch_script(settings, offline=True))
        offline_script.chmod(0o755)
        (staging / 'README.txt').write_text(
            'LEAP manual review bundle\n\n' + '\n'.join(k + ': ' + v for k, v in DESCRIPTIONS.items()) +
            '\n\nPositive/negative tracks are combined across the selected chromosomes; unsplit contains both strands.\n'
            'Coordinates, IDs and metadata are retained; only internal-track transcript_class is added for review.\n'
            'External and internal views can contain the same transcript with different terminal boundaries.\n'
            'Internal tracks require positive extension length at an internal end; classification lists also include unchanged internal models.\n'
            'Raw/leaped views retain gene rows; processed external/internal views are transcript-oriented.\n'
            'Header-only files are retained but omitted from the samplesheet. See track_manifest.tsv.\n'
            'processing/ contains corrected external stats and classification/dropped-no-CDS tables.\n\n'
            'When supplied, reads/ contains positive, negative and unsplit mapped-input BAMs with CSI indexes.\n'
            'BAM flag 0x10 defines orientation. Alignment counts are records, not molecules or selected-extension support.\n'
            'Empty alignment views remain on disk but are omitted from JBrowser sheets.\n'
            'Alignment rows use reference_file; update those entries too when changing the reference after generation.\n'
            'Before launch: complete any EDIT_ME cells in samplesheet.csv; check the server reference and assembly.\n'
            'Edit JBrowser checkout, FTP roots, hub email and Slurm settings in run_jbrowser.sbatch.\n'
            'Submit with: sbatch run_jbrowser.sbatch\n'
            'The launch script prepares JBrowser output and FTP copy scripts. Publishing is a later manual step.\n'
            'For offline review, submit: sbatch run_jbrowser_offline.sbatch (no FTP settings required).\n'
            'The offline script refreshes samplesheet_offline.csv from samplesheet.csv, forcing every row offline.\n'
            'Edit reference paths in the shared samplesheet.csv; offline output/work/reports use separate _offline names.\n'
            'Paths refer to this published bundle. If you move it, update samplesheet.csv and REVIEW_DIR in both scripts.\n')
        staging.rename(output)


def one_file(folder, stem, suffixes):
    paths = [Path(folder) / (stem + suffix) for suffix in suffixes]
    paths = [p for p in paths if p.is_file()]
    if len(paths) != 1:
        raise ValueError('Missing or ambiguous input: ' + str(Path(folder) / stem))
    return paths[0]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    raw = sub.add_parser('raw')
    raw.add_argument('chromosome'); raw.add_argument('annotation'); raw.add_argument('output')
    part = sub.add_parser('partition')
    for name in ('chromosome', 'strand', 'raw', 'extended', 'stats', 'output'):
        part.add_argument(name)
    bundle = sub.add_parser('bundle')
    for name in ('partitions', 'manifest', 'settings', 'output'):
        bundle.add_argument(name)
    bundle.add_argument('--reference-lengths')
    bundle.add_argument('--read-evidence')
    run = sub.add_parser('run')
    run.add_argument('leap_output_dir', type=Path); run.add_argument('output', type=Path)
    run.add_argument('--chromosomes', nargs='+', required=True)
    run.add_argument('--strand', choices=('positive', 'negative', 'both'), default='both')
    run.add_argument('--reference-lengths'); run.add_argument('--reference-genome')
    run.add_argument('--reads', type=Path); run.add_argument('--reads-index', type=Path)
    run.add_argument('--read-threads', type=int, default=1)
    for name in ('genome-id', 'species', 'jbrowser-pipeline', 'public-base-url', 'ftp-target-dir', 'hub-email'):
        run.add_argument('--' + name, default='')
    run.add_argument('--ftp-project', default='LEAP_review')
    run.add_argument('--publication', choices=('public', 'offline'), default='public')
    args = parser.parse_args()
    try:
        if args.command == 'raw':
            raw_partition(args.chromosome, args.annotation, args.output)
        elif args.command == 'partition':
            extension_partition(args.chromosome, args.strand, args.raw, args.extended, args.stats, args.output)
        elif args.command == 'bundle':
            build_bundle(sorted(Path(args.partitions).iterdir()), json.loads(Path(args.manifest).read_text()),
                         json.loads(Path(args.settings).read_text()), args.output, args.reference_lengths, args.read_evidence)
        else:
            if bool(args.reads) != bool(args.reads_index):
                raise ValueError('--reads and --reads-index must be supplied together')
            if args.output.exists() or args.output.is_symlink():
                raise ValueError('Output must be a new directory: ' + str(args.output))
            strands = list(STRANDS) if args.strand == 'both' else [args.strand]
            scopes = [(c, s) for c in args.chromosomes for s in strands]
            settings = {k: v for k, v in vars(args).items() if k in ('genome_id', 'species', 'jbrowser_pipeline',
                'public_base_url', 'ftp_target_dir', 'hub_email', 'ftp_project', 'publication', 'reference_genome')}
            if settings.get('reference_genome'):
                settings['reference_genome'] = str(Path(settings['reference_genome']).resolve(strict=True))
            with tempfile.TemporaryDirectory(prefix='leap-review-') as temporary:
                parts = Path(temporary) / 'partitions'
                parts.mkdir()
                for c in args.chromosomes:
                    annotation = one_file(args.leap_output_dir / 'raw_annotation', 'current_annotation_' + c,
                        ('.gff3', '.gff', '.gff3.gz', '.gff.gz'))
                    raw_partition(c, annotation, parts / ('raw_' + c))
                    for s in strands:
                        extended = one_file(args.leap_output_dir / 'extended_transcripts', c + '_' + s + '_extended',
                            ('.gff', '.gff3', '.gff.gz', '.gff3.gz'))
                        stats = one_file(args.leap_output_dir / 'extended_transcripts_stats', c + '_' + s + '_stats', ('.tsv', '.tsv.gz'))
                        extension_partition(c, s, annotation, extended, stats, parts / ('extended_' + c + '_' + s))
                evidence = None
                if args.reads:
                    helper = read_processor()
                    chunks = helper.partition_source(args.reads, args.reads_index, args.chromosomes,
                        Path(temporary) / 'read-inputs', args.read_threads)
                    evidence = Path(temporary) / 'read_evidence'
                    helper.prepare_reads(chunks, scopes, evidence, args.read_threads)
                    provenance = json.loads((evidence / 'read-summary.json').read_text())
                    provenance['original_input'] = dict(bam=str(args.reads.resolve()), index=str(args.reads_index.resolve()))
                    (evidence / 'read-summary.json').write_text(json.dumps(provenance, indent=2) + '\n')
                build_bundle(sorted(parts.iterdir()), scopes, settings, args.output, args.reference_lengths, evidence)
        return 0
    except (OSError, ValueError, KeyError, TypeError, EOFError, csv.Error, subprocess.SubprocessError) as error:
        print('ERROR: ' + str(error), file=sys.stderr)
        return 2


if __name__ == '__main__':
    sys.exit(main())
