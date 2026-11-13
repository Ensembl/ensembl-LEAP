#!/usr/bin/env python3
"""Nextflow adapter: retain every check result, then summarize and gate the run."""

import argparse
import csv
import gzip
import html
import json
import os
import re
import subprocess
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

CHECKS = {
    'boundary': 'check_for_shorter_transcripts.py',
    'internal_exons': 'check_internal_exon_positions.py',
    'preservation': 'check_transcript_preservation.py',
    'stats_gff': 'check_stats_gff_consistency.py',
    'structure': 'check_gff_structural_integrity.py',
    'completeness': 'check_output_completeness.py',
    'partial_protection': 'check_partial_protection.py',
    'reference_coordinates': 'check_reference_coordinates.py',
    'cds_overlaps': 'check_cds_overlaps.py',
}
LEGACY_EMPTY = {'boundary', 'internal_exons', 'preservation', 'stats_gff', 'structure', 'cds_overlaps'}
SUMMARY_FIELDS = ['chromosome', 'strand', 'check', 'status', 'exit_code', 'message',
                  'report', 'stdout', 'stderr', 'elapsed_seconds']
VALID_STATUSES = {0: {'PASS', 'NOT_APPLICABLE'}, 1: {'FAIL'}, 2: {'INPUT_ERROR'}}


def safe_chromosome(value):
    if not re.fullmatch(r'[A-Za-z0-9_.-]+', value) or value in ('.', '..'):
        raise ValueError('Invalid chromosome name: ' + value)
    return value


def prefix(chromosome, strand, check):
    return '{}_{}_{}'.format(safe_chromosome(chromosome), strand, check)


def one_file(directory, suffixes, stem):
    names = {stem + suffix for suffix in suffixes}
    paths = sorted(p for p in directory.iterdir() if p.is_file() and p.name in names)
    if len(paths) != 1:
        raise ValueError('{}: expected exactly one corresponding input file; found {}'.format(directory, len(paths)))
    return paths[0]


def open_text(path):
    return (gzip.open if str(path).endswith('.gz') else open)(path, 'rt', encoding='utf-8', newline='')


def empty_gff(path):
    with open_text(path) as handle:
        for line in handle:
            if line.strip() == '##FASTA':
                break
            if line.strip() and not line.startswith('#'):
                return False
    return True


def invoke(command, stdout, stderr):
    """Checker failures are data, not task exit codes; final gate handles failure."""
    with stdout.open('w', encoding='utf-8') as out, stderr.open('w', encoding='utf-8') as err:
        result = subprocess.run(command, stdout=out, stderr=err,
                                env=dict(os.environ, PYTHONDONTWRITEBYTECODE='1'))
    return result.returncode


def completeness_command(context, checkers, chromosome, strand, report):
    return [sys.executable, str(checkers / CHECKS['completeness']),
            str(context / 'raw_annotation'), str(context / 'extended_transcripts'), str(report),
            '--extension-stats', str(context / 'extended_transcripts_stats'),
            '--chromosomes', chromosome, '--strand', strand]


def run_check(args):
    safe_chromosome(args.chromosome)
    name = prefix(args.chromosome, args.strand, args.check)
    output = Path(args.output)
    if output.exists():
        raise ValueError('Refusing to overwrite an existing check result: ' + str(output))
    context, checkers = Path(args.context), Path(args.checkers)
    output.mkdir(parents=True)
    report, stdout, stderr = (output / n for n in ('report.csv', 'stdout.txt', 'stderr.txt'))
    stdout.write_text('')
    stderr.write_text('')
    status = dict(schema_version=1, chromosome=args.chromosome, strand=args.strand, check=args.check,
                  status='INPUT_ERROR', exit_code=2, checker_exit_code=None, message='', command=[],
                  report=name + '/report.csv', stdout=name + '/stdout.txt', stderr=name + '/stderr.txt')
    started = time.monotonic()
    try:
        for folder in ('raw_annotation', 'extended_transcripts', 'extended_transcripts_stats',
                       'fixed_extensions_before_partial_protection', 'fixed_extensions',
                       'partial_protection', 'natural_overlaps'):
            (context / folder).mkdir(parents=True, exist_ok=True)
        checker = checkers / CHECKS[args.check]
        if not checker.is_file():
            raise ValueError('Missing staged datacheck script: ' + str(checker))
        if args.check == 'completeness':
            command = completeness_command(context, checkers, args.chromosome, args.strand, report)
        elif args.check == 'partial_protection':
            command = [sys.executable, str(checker), str(context), str(report),
                       '--chromosomes', args.chromosome, '--strand', args.strand]
        else:
            raw = one_file(context / 'raw_annotation', ('.gff', '.gff3', '.gff.gz', '.gff3.gz'), 'current_annotation_' + args.chromosome)
            extended = one_file(context / 'extended_transcripts', ('.gff', '.gff3', '.gff.gz', '.gff3.gz'), args.chromosome + '_' + args.strand + '_extended')
            if args.check in LEGACY_EMPTY and empty_gff(extended):
                preflight = output / 'empty_partition_completeness.csv'
                pre_command = completeness_command(context, checkers, args.chromosome, args.strand, preflight)
                result = invoke(pre_command, output / 'empty_partition_stdout.txt', output / 'empty_partition_stderr.txt')
                if result == 0:
                    with preflight.open(newline='') as handle:
                        rows = list(csv.DictReader(handle))
                    if not any(row['expected_retention'] == 'True' or row['gff_present'] == 'True'
                               or row['stats_present'] == 'True' for row in rows):
                        status.update(status='NOT_APPLICABLE', exit_code=0, command=pre_command,
                                      message='No output features and zero eligible transcripts; verified by completeness check')
                        stdout.write_text(status['message'] + '\n')
                        command = None
                    else:
                        raise ValueError('Empty-partition completeness preflight gave inconsistent membership')
                else:
                    command = []  # Run the original checker below; never hide malformed/missing models.
            else:
                command = []
            if status['status'] != 'NOT_APPLICABLE':
                if args.check == 'reference_coordinates':
                    reference = Path(args.reference_lengths)
                    if not reference.is_file():
                        raise ValueError('Missing reference lengths: ' + str(reference))
                    combined = output.with_name(output.name + '_reference_inputs')
                    combined.mkdir()
                    for path in (raw, extended):
                        # Links stay within declared Nextflow inputs; no published-output dependency.
                        (combined / path.name).symlink_to(os.path.relpath(path, combined))
                    command = [sys.executable, str(checker), str(combined), str(reference), str(report)]
                elif args.check == 'structure':
                    command = [sys.executable, str(checker), str(extended), str(report)]
                else:
                    command = [sys.executable, str(checker), str(raw), str(extended), str(report)]
                    if args.check in ('stats_gff', 'cds_overlaps'):
                        stats = one_file(context / 'extended_transcripts_stats', ('.tsv', '.tsv.gz'), args.chromosome + '_' + args.strand + '_stats')
                        command += ['--extension-stats', str(stats), '--leap-output-dir', str(context)]
        if command is not None:
            status['command'] = command
            result = invoke(command, stdout, stderr)
            code = result if result in (0, 1, 2) else 2
            status.update(checker_exit_code=result, exit_code=code,
                          status={0: 'PASS', 1: 'FAIL', 2: 'INPUT_ERROR'}[code])
            log = (stderr if code == 2 else stdout).read_text(errors='replace').strip()
            status['message'] = log.splitlines()[0][:2000] if log else 'Datacheck exited with code {}'.format(result)
    except Exception as error:
        # The adapter always records unexpected checker/adapter execution errors too.
        status.update(status='INPUT_ERROR', exit_code=2, message=str(error))
        with stderr.open('a', encoding='utf-8') as handle:
            handle.write('ERROR: ' + str(error) + '\n')
    if status['status'] in ('PASS', 'FAIL') and not report.is_file():
        status.update(status='INPUT_ERROR', exit_code=2, message='Checker did not produce its required CSV report; see logs')
    if status['exit_code'] and 'Traceback (most recent call last):' in stderr.read_text(errors='replace'):
        status.update(status='INPUT_ERROR', exit_code=2, message='Unexpected checker execution error; see stderr')
    if not report.exists():
        with report.open('w', newline='', encoding='utf-8') as handle:
            writer = csv.DictWriter(handle, fieldnames=['status', 'details'])
            writer.writeheader()
            writer.writerow(dict(status=status['status'], details=status['message']))
    status['elapsed_seconds'] = round(time.monotonic() - started, 3)
    (output / (name + '.status.json')).write_text(json.dumps(status, indent=2) + '\n')
    print('{} {}: {}'.format(name, status['status'], status['message']))
    return 0


def summarize(args):
    output = Path(args.output)
    if output.exists():
        raise ValueError('Refusing to overwrite an existing summary: ' + str(output))
    manifest = json.loads(Path(args.expected).read_text())
    scopes = manifest['scopes']
    records, findings = defaultdict(list), []
    expected = set()
    for chrom, strand in scopes:
        safe_chromosome(chrom)
        if strand not in ('positive', 'negative'):
            raise ValueError('Invalid strand in expected scope')
        for name in CHECKS:
            key = (chrom, strand, name)
            if key in expected:
                raise ValueError('Duplicate expected scope: ' + str(key))
            expected.add(key)
    if not expected:
        findings.append(dict(check='workflow_scope', status='INPUT_ERROR', exit_code=2,
                             message='No expected chromosome/strand partitions were selected'))
    for issue in manifest.get('issues', []):
        findings.append(dict(check='workflow_context', status='INPUT_ERROR', exit_code=2, message=str(issue)))
    directory = Path(args.statuses)
    for path in sorted(directory.glob('*.status.json')):
        try:
            record = json.loads(path.read_text())
            key = tuple(record[k] for k in ('chromosome', 'strand', 'check'))
            code = record['exit_code']
            if type(code) is not int or code not in VALID_STATUSES or record['status'] not in VALID_STATUSES[code]:
                raise ValueError('Invalid exit code/status combination')
            if key not in expected:
                raise ValueError('Unexpected result scope: ' + str(key))
            for field in ('report', 'stdout', 'stderr', 'message', 'elapsed_seconds'):
                if field not in record:
                    raise ValueError('Missing result field: ' + field)
            records[key].append(record)
        except (ValueError, KeyError, TypeError) as error:
            findings.append(dict(check='workflow_result', status='INPUT_ERROR', exit_code=2,
                                 message='{}: {}'.format(path.name, error)))
    rows = []
    for key in sorted(expected):
        items = records.get(key, [])
        if len(items) == 1:
            rows.append(items[0])
        else:
            rows.append(dict(zip(('chromosome', 'strand', 'check'), key), status='INPUT_ERROR', exit_code=2,
                             message='Expected exactly one result; received {}'.format(len(items))))
    rows += findings
    counts = Counter(row['status'] for row in rows)
    code = max((row['exit_code'] for row in rows), default=2)
    verdict = dict(schema_version=1, status={0: 'PASS', 1: 'FAIL', 2: 'INPUT_ERROR'}[code], exit_code=code,
                   expected_checks=len(expected), received_status_files=sum(len(v) for v in records.values()),
                   counts=dict(counts), results=rows)
    output.mkdir(parents=True)
    with (output / 'summary.csv').open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_FIELDS, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)
    (output / 'summary.json').write_text(json.dumps(verdict, indent=2) + '\n')
    (output / 'verdict.json').write_text(json.dumps({k: v for k, v in verdict.items() if k != 'results'}, indent=2) + '\n')
    page = ['<!doctype html><html lang="en"><meta charset="utf-8"><title>LEAP datachecks</title>',
            '<style>body{font:16px system-ui;margin:32px;line-height:1.5}td,th{padding:8px;border-bottom:1px solid #ddd;text-align:left}.PASS{color:#08633c}.FAIL,.INPUT_ERROR{color:#a72020}td:last-child{max-width:750px;overflow-wrap:anywhere}</style>',
            '<h1>LEAP datachecks: {}</h1><p>{} checks expected. {}</p>'.format(verdict['status'], len(expected), html.escape(str(dict(counts)))),
            '<p><a href="summary.csv">CSV summary</a> · <a href="summary.json">Full JSON summary</a></p>',
            '<table><tr><th>Chromosome</th><th>Strand</th><th>Check</th><th>Result</th><th>Files</th><th>Details</th></tr>']
    for row in rows:
        links = []
        for field in ('report', 'stdout', 'stderr'):
            if row.get(field):
                # Only adapter-created relative output paths may become links.
                target = row[field]
                if not re.fullmatch(r'[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+', target):
                    continue
                links.append('<a href="../{}">{}</a>'.format(html.escape(target, quote=True), field))
        cells = [html.escape(str(row.get(k, ''))) for k in ('chromosome', 'strand', 'check', 'status')]
        page.append('<tr class="{}"><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td></tr>'.format(
            row['status'], *cells, ' · '.join(links), html.escape(row.get('message', ''))))
    page.append('</table></html>')
    (output / 'summary.html').write_text(''.join(page))
    print('LEAP datachecks: {} — {}'.format(verdict['status'], dict(counts)))
    return 0  # Publish the summary before the separate final gate fails the workflow.


def gate(args):
    verdict = json.loads(Path(args.verdict).read_text())
    code = verdict['exit_code']
    if code not in (0, 1, 2):
        raise ValueError('Invalid verdict exit code')
    if code:
        print('LEAP datachecks failed: {}. All check reports and the combined summary were retained.'.format(verdict['counts']), file=sys.stderr)
        return code
    Path('datachecks.passed').write_text('All required datachecks passed or were explicitly not applicable.\n')
    print('All LEAP datachecks passed.')
    return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='action', required=True)
    check = commands.add_parser('check')
    check.add_argument('--check', choices=list(CHECKS), required=True)
    check.add_argument('--chromosome', required=True)
    check.add_argument('--strand', choices=('positive', 'negative'), required=True)
    check.add_argument('--context', required=True)
    check.add_argument('--checkers', required=True)
    check.add_argument('--reference-lengths', required=True)
    check.add_argument('--output', required=True)
    summary = commands.add_parser('summarize')
    summary.add_argument('--expected', required=True)
    summary.add_argument('--statuses', required=True)
    summary.add_argument('--output', required=True)
    final = commands.add_parser('gate')
    final.add_argument('--verdict', required=True)
    args = parser.parse_args()
    try:
        return {'check': run_check, 'summarize': summarize, 'gate': gate}[args.action](args)
    except (OSError, ValueError, KeyError, csv.Error, TypeError) as error:
        print('ERROR: ' + str(error), file=sys.stderr)
        return 2


if __name__ == '__main__':
    sys.exit(main())
