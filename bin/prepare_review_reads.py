#!/usr/bin/env python3
"""Prepare mapped input alignment evidence for LEAP's manual-review bundle."""

import argparse
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile

sys.dont_write_bytecode = True
VIEWS = ('positive', 'negative', 'unsplit')


def command(*arguments):
    result = subprocess.run(['samtools', *map(str, arguments)], capture_output=True, text=True)
    if result.returncode:
        raise ValueError('Samtools {} failed: {}'.format(arguments[0], result.stderr.strip() or result.returncode))
    return result.stdout


def scopes_by_chromosome(scopes):
    result = {}
    for chromosome, strand in scopes:
        chromosome = str(chromosome)
        if (not re.fullmatch(r'[A-Za-z0-9_.-]+', chromosome) or chromosome in ('.', '..')
                or strand not in ('positive', 'negative')):
            raise ValueError('Invalid read-evidence scope')
        if strand in result.setdefault(chromosome, set()):
            raise ValueError('Duplicate read-evidence scope: ' + chromosome + ' ' + strand)
        result[chromosome].add(strand)
    if not result:
        raise ValueError('Read-evidence scope must not be empty')
    return result


def require_pair(bam, index):
    bam, index = Path(bam), Path(index)
    if not bam.name.endswith('.bam'):
        raise ValueError('Read evidence requires BAM inputs: ' + str(bam))
    for path in (bam, index):
        if not path.is_file() or not os.access(path, os.R_OK):
            raise ValueError('Missing or unreadable read-evidence input: ' + str(path))
    command('quickcheck', '-v', bam)


def header_context(bam):
    references, groups, names = [], {}, set()
    for line in command('view', '-H', '--no-PG', bam).splitlines():
        if not line.startswith(('@SQ\t', '@RG\t')):
            continue
        fields = dict(field.split(':', 1) for field in line.split('\t')[1:])
        if line.startswith('@SQ'):
            name = fields.get('SN', '')
            if not name or name in names or not fields.get('LN', '').isdigit() or int(fields['LN']) <= 0:
                raise ValueError('Invalid or duplicate BAM reference dictionary entry: ' + str(bam))
            names.add(name)
            references.append(tuple(sorted(fields.items())))
        else:
            identifier = fields.get('ID')
            if not identifier or identifier in groups:
                raise ValueError('Invalid or duplicate BAM read-group ID: ' + str(bam))
            groups[identifier] = fields
    if not references:
        raise ValueError('BAM has no reference dictionary: ' + str(bam))
    return references, groups, names


def filter_sort(source, destination, flags, threads):
    # A file-backed stderr avoids deadlocks on a large stream of diagnostics.
    with tempfile.TemporaryFile(mode='w+') as errors:
        producer = subprocess.Popen(['samtools', 'view', '--no-PG', '-u', *flags, str(source)],
                                    stdout=subprocess.PIPE, stderr=errors)
        try:
            result = subprocess.run(['samtools', 'sort', '--no-PG', '-@', str(threads),
                                     '-o', str(destination), '-'], stdin=producer.stdout,
                                    stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
        finally:
            producer.stdout.close()
            code = producer.wait()
        errors.seek(0)
        details = errors.read().strip()
        if code or result.returncode:
            raise ValueError('Samtools read filtering/sorting failed: ' + (details or result.stderr.strip()))


def prepare_reads(partitions, scopes, output, threads=1):
    """Each input BAM must contain mapped records from its declared chromosome only."""
    if not shutil.which('samtools'):
        raise ValueError('Required command not found: samtools')
    if threads < 1:
        raise ValueError('Read-evidence threads must be positive')
    output = Path(output).absolute()
    if output.exists() or output.is_symlink():
        raise ValueError('Read-evidence output must be a new directory: ' + str(output))
    expected = scopes_by_chromosome(scopes)
    actual = [str(row['chromosome']) for row in partitions]
    if len(set(actual)) != len(actual) or set(actual) != set(expected):
        raise ValueError('Missing, duplicate or unexpected chromosome BAM/index pairs')
    context, sources = None, []
    for row in sorted(partitions, key=lambda row: str(row['chromosome'])):
        chromosome = str(row['chromosome'])
        bam, index = Path(row['bam']), Path(row['index'])
        require_pair(bam, index)
        references, groups, names = header_context(bam)
        if chromosome not in names:
            raise ValueError('Declared chromosome absent from BAM header: ' + chromosome)
        if context is not None and context != (references, groups):
            raise ValueError('Incompatible BAM reference dictionaries or read groups')
        context = references, groups
        total = int(command('view', '-c', '-F', '4', bam))
        indexed = int(command('view', '-c', '-F', '4', '-X', bam, index, chromosome))
        if total != indexed:
            raise ValueError('BAM contains mapped records outside chromosome {} or its index is inconsistent: {}'.format(chromosome, bam))
        sources.append(dict(chromosome=chromosome, bam=str(bam.resolve()), index=str(index.resolve()),
                            mapped_alignment_records=total))
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.' + output.name + '-', dir=output.parent) as directory:
        temporary = Path(directory)
        staging = temporary / 'reads'
        staging.mkdir()
        chunks = []
        for n, row in enumerate(sorted(partitions, key=lambda row: str(row['chromosome']))):
            chromosome = str(row['chromosome'])
            selected = expected[chromosome]
            flags = ['-F', '4'] if len(selected) == 2 else (
                ['-F', '20'] if 'positive' in selected else ['-F', '4', '-f', '16'])
            chunk = temporary / ('scope_' + str(n) + '.bam')
            filter_sort(row['bam'], chunk, flags, threads)
            sources[n]['selected_alignment_records'] = int(command('view', '-c', chunk))
            chunks.append(chunk)
        # All chunks derive from compatible headers. Keep RG and record tag IDs;
        # @PG processing-history lines may differ between chromosome tasks.
        command('merge', '--no-PG', '-c', '-p', '-@', threads, '-o', staging / 'unsplit.bam', *chunks)
        for strand, flags in (('positive', ['-F', '20']), ('negative', ['-F', '4', '-f', '16'])):
            filter_sort(staging / 'unsplit.bam', staging / (strand + '.bam'), flags, threads)
        counts = {}
        for view in VIEWS:
            bam = staging / (view + '.bam')
            index = Path(str(bam) + '.csi')
            command('quickcheck', '-v', bam)
            command('index', '-c', '-@', threads, bam, index)
            counts[view] = int(command('view', '-c', bam))
            indexed = int(command('view', '-c', '-X', bam, index, *sorted(expected)))
            if indexed != counts[view]:
                raise ValueError('Prepared read-evidence index does not reproduce alignment counts: ' + view)
        if counts['positive'] + counts['negative'] != counts['unsplit']:
            raise ValueError('Prepared strand alignment counts do not sum to unsplit')
        if counts['unsplit'] != sum(row['selected_alignment_records'] for row in sources):
            raise ValueError('Prepared alignment counts differ from the selected input records')
        summary = dict(scopes=sorted((c, s) for c, strands in expected.items() for s in strands),
            sources=sources, counts=counts, samtools_version=command('--version').splitlines()[0],
            description='All mapped input alignment records in scope; orientation uses BAM flag 0x10. Counts are records, not unique molecules or selected-extension support.')
        (staging / 'read-summary.json').write_text(json.dumps(summary, indent=2) + '\n')
        staging.rename(output)
    print('Prepared read evidence: ' + str(counts))


def partition_source(bam, index, chromosomes, output, threads=1):
    """Select chromosome copies from a whole-run BAM using its explicit index."""
    if not shutil.which('samtools'):
        raise ValueError('Required command not found: samtools')
    require_pair(bam, index)
    _, _, names = header_context(bam)
    if not set(chromosomes) <= names:
        raise ValueError('Selected chromosome absent from read-evidence BAM header')
    output = Path(output)
    output.mkdir()
    result = []
    for n, chromosome in enumerate(chromosomes):
        path = output / ('chromosome_' + str(n) + '.bam')
        command('view', '--no-PG', '-b', '-F', '4', '-o', path, '-X', bam, index, chromosome)
        sorted_path = output / ('sorted_' + str(n) + '.bam')
        command('sort', '--no-PG', '-@', threads, '-o', sorted_path, path)
        path.unlink()
        output_index = Path(str(sorted_path) + '.csi')
        command('index', '-c', sorted_path, output_index)
        result.append(dict(chromosome=chromosome, bam=str(sorted_path), index=str(output_index)))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('manifest', help='JSON object containing scopes and chromosome BAM/index pairs')
    parser.add_argument('output')
    parser.add_argument('--threads', type=int, default=1)
    args = parser.parse_args()
    try:
        manifest = json.loads(Path(args.manifest).read_text())
        prepare_reads(manifest['partitions'], manifest['scopes'], args.output, args.threads)
        return 0
    except (OSError, ValueError, KeyError, TypeError, subprocess.SubprocessError) as error:
        print('ERROR: ' + str(error), file=sys.stderr)
        return 2


if __name__ == '__main__':
    sys.exit(main())
