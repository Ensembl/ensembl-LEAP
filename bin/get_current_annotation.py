#!/usr/bin/env python3
"""Export a chromosome and apply LEAP's existing readthrough/OTT/blacklist rules."""

import argparse
from contextlib import redirect_stderr, redirect_stdout
import io
from pathlib import Path
import re
import shutil
import sys
import tempfile

sys.dont_write_bytecode = True

from annotation_db import load_connection, redact_connection_text
from gxf_dump import exporter
from gxf_dump.regions import fetch_regions


def exclusion_ids(secrets, chromosome):
    from sqlalchemy import bindparam, text
    with exporter.sql_connection(secrets) as connection:
        regions = fetch_regions(connection, chromosome=chromosome)
        region_ids = tuple(row['seq_region_id'] for row in regions)
        readthroughs = text('''
            SELECT DISTINCT g.stable_id FROM transcript t
            JOIN transcript_attrib ta ON t.transcript_id = ta.transcript_id
            JOIN gene g ON t.gene_id = g.gene_id
            WHERE t.is_current = 1 AND t.seq_region_id IN :region_ids
              AND LOWER(ta.value) = 'readthrough'
        ''').bindparams(bindparam('region_ids', expanding=True))
        blacklist = text('''
            SELECT DISTINCT t.stable_id FROM transcript t
            JOIN transcript_attrib ta ON t.transcript_id = ta.transcript_id
            JOIN attrib_type a ON a.attrib_type_id = ta.attrib_type_id
            WHERE t.is_current = 1 AND t.biotype = 'protein_coding'
              AND t.seq_region_id IN :region_ids AND a.code = 'remark'
              AND LOWER(ta.value) IN ('upstream atg', 'cage supported tss', 'rp supported tis')
        ''').bindparams(bindparam('region_ids', expanding=True))
        parameters = dict(region_ids=region_ids)
        return ({row[0] for row in connection.execute(readthroughs, parameters)},
                {row[0] for row in connection.execute(blacklist, parameters)})


def filter_annotation(source, destination, readthroughs=(), remove_ott=False):
    removed_genes, removed_transcripts = set(), set()
    def records():
        with Path(source).open(encoding='utf-8') as handle:
            for line in handle:
                if not line.strip() or line.startswith('#'):
                    continue
                fields = line.rstrip('\r\n').split('\t')
                if len(fields) != 9:
                    raise ValueError('Exporter produced a malformed GFF3 record')
                attrs = dict(item.split('=', 1) for item in fields[8].split(';') if item)
                yield line, fields, attrs
    if remove_ott:
        for line, fields, attrs in records():
            if fields[2] == 'gene' and attrs['gene_id'].startswith('OTT'):
                removed_genes.add(attrs['gene_id'])
            elif fields[2] == 'transcript' and attrs['transcript_id'].startswith('OTT'):
                removed_transcripts.add(attrs['transcript_id'])
    kept = 0
    with Path(destination).open('w', encoding='utf-8') as handle:
        for line, fields, attrs in records():
            gene_id = attrs['gene_id']
            if re.sub(r'\.\d+$', '', gene_id) in readthroughs:
                continue
            if (gene_id in removed_genes or attrs.get('transcript_id') in removed_transcripts
                    or attrs.get('Parent') in removed_genes | removed_transcripts):
                continue
            handle.write(line)
            kept += 1
    return kept


def run(db_secrets, chromosome, output_directory=Path('.')):
    if not re.fullmatch(r'[A-Za-z0-9_.-]+', chromosome) or chromosome in ('.', '..'):
        raise ValueError('Invalid chromosome name')
    secrets = load_connection(db_secrets)
    output_directory = Path(output_directory)
    names = ['db_out.gff3', 'db_out_RR.gff3', 'current_annotation_' + chromosome + '.gff3',
             'readthroughs.txt', 'do_not_extend_transcripts_positive.txt', 'do_not_extend_transcripts_negative.txt']
    for name in names:
        output = output_directory / name
        if output.resolve() == Path(db_secrets).resolve() or (output.exists() and output.samefile(db_secrets)):
            raise ValueError('Annotation output would overwrite the database secrets file')
    with tempfile.TemporaryDirectory(prefix='.annotation-', dir=output_directory) as temporary:
        staging = Path(temporary)
        # Exporter errors can contain database-driver text. Capture them rather
        # than leaking connection values or SQL parameters into task logs.
        errors = io.StringIO()
        try:
            with redirect_stdout(io.StringIO()), redirect_stderr(errors):
                exporter.export_main(['--db-secrets', str(db_secrets), '--chromosome', chromosome,
                                      '--vega', '--format', 'gff3', '--output', str(staging / names[0])])
                readthroughs, blacklist = exclusion_ids(secrets, chromosome)
        except (Exception, SystemExit) as error:
            details = errors.getvalue().strip() if isinstance(error, SystemExit) else str(error)
            raise ValueError('Database export/filter query failed: ' + redact_connection_text(details, secrets)) from None
        if errors.getvalue().strip():
            print(redact_connection_text(errors.getvalue().strip(), secrets), file=sys.stderr)
        filter_annotation(staging / names[0], staging / names[1], readthroughs)
        count = filter_annotation(staging / names[1], staging / names[2], remove_ott=True)
        if not count:
            raise ValueError('Current annotation is empty after readthrough/OTT filtering')
        (staging / names[3]).write_text('stable_id\n' + ''.join(str(v) + '\n' for v in sorted(readthroughs)))
        (staging / names[4]).write_text('stable_id\n' + ''.join(str(v) + '\n' for v in sorted(blacklist)))
        shutil.copyfile(staging / names[4], staging / names[5])
        for name in names:
            (staging / name).replace(output_directory / name)
    print('Prepared chromosome {}: {} annotation feature records'.format(chromosome, count))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('db_secrets')
    parser.add_argument('chromosome')
    args = parser.parse_args()
    try:
        run(args.db_secrets, args.chromosome)
        return 0
    except ValueError as error:
        print('ERROR: ' + str(error), file=sys.stderr)
        return 2
    except OSError:
        print('ERROR: Cannot read/write annotation task files; check input and output permissions', file=sys.stderr)
        return 2


if __name__ == '__main__':
    sys.exit(main())
