#!/usr/bin/env python3
"""Prepare CAGE/poly(A) interval views for manual JBrowser review."""

import argparse
from collections import Counter
import csv
import gzip
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
from urllib.parse import unquote

sys.dont_write_bytecode = True

STRANDS = {'positive': '+', 'negative': '-'}
VIEWS = ('positive', 'negative', 'unsplit')


def cage_preset_path(name):
    candidates = (Path(__file__).with_name(name+'.tsv'),
                  Path(__file__).resolve().parents[1]/'assets/jbrowser_region_presets'/(name+'.tsv'))
    path = next((p for p in candidates if p.is_file()), None)
    if path is None:
        raise ValueError('Missing CAGE browser format definition: '+name)
    return path


def scopes_by_chromosome(scopes):
    result = {}
    for chromosome, strand in scopes:
        if (not isinstance(chromosome, str) or not re.fullmatch(r'[A-Za-z0-9_.-]+', chromosome)
                or chromosome in ('.', '..') or strand not in STRANDS):
            raise ValueError('Invalid site-evidence scope')
        if strand in result.setdefault(chromosome, set()):
            raise ValueError('Duplicate site-evidence scope')
        result[chromosome].add(strand)
    if not result:
        raise ValueError('Site evidence requires a nonempty scope')
    return result


def source_records(path):
    name = path.name.lower().removesuffix('.gz')
    kind = 'bed' if name.endswith('.bed') else 'gff' if name.endswith(('.gff', '.gff3')) else None
    if kind is None:
        raise ValueError('Site evidence requires BED, GFF or GFF3 (optionally gzip): ' + str(path))
    with (gzip.open if path.name.endswith('.gz') else open)(path, 'rt', encoding='utf-8-sig') as handle:
        for number, line in enumerate(handle, 1):
            if line.strip() == '##FASTA':
                break
            if not line.strip() or line.startswith(('#', 'track ', 'browser ')):
                continue
            f = line.rstrip('\r\n').split('\t')
            context = '{}:{}'.format(path, number)
            try:
                if kind == 'bed':
                    if len(f) < 3:
                        raise ValueError()
                    start, end = int(f[1]), int(f[2])
                    strand = f[5] if len(f) >= 6 else '.'
                    label = f[3] if len(f) >= 4 and f[3] not in ('', '.') else 'site_' + str(number)
                    score = f[4] if len(f) >= 5 else '.'
                else:
                    if len(f) != 9:
                        raise ValueError()
                    start, end = int(f[3]) - 1, int(f[4])
                    strand, score = f[6], f[5]
                    attrs = {}
                    if f[8] != '.':
                        for item in f[8].split(';'):
                            if not item:
                                continue
                            key, separator, value = item.partition('=')
                            if not separator or not key or key in attrs:
                                raise ValueError()
                            attrs[key] = unquote(value)
                    label = attrs.get('Name') or attrs.get('ID') or f[2] + '_' + str(number)
                if not 0 <= start < end or strand not in ('+', '-', '.', '?') or f[0] in ('', '.'):
                    raise ValueError()
                if score != '.' and not math.isfinite(float(score)):
                    raise ValueError()
            except ValueError:
                raise ValueError(context + ': malformed site record, coordinates, score or strand') from None
            chromosome = f[0][3:] if f[0].startswith('chr') else f[0]
            if not chromosome:
                raise ValueError(context + ': missing chromosome')
            yield dict(chromosome=chromosome, start=start, end=end, strand='.' if strand=='?' else strand,
                name=re.sub(r'[\t\r\n]', ' ', label), score=score, fields=f, line=number, format=kind)


def selected(record, scopes):
    if record['chromosome'] not in scopes:
        return False
    return record['strand']=='.' or any(record['strand']==STRANDS[s] for s in scopes[record['chromosome']])


def display_fields(record, region_format):
    if region_format in ('cageBed6', 'cageBed9'):
        base = 'bed9' if region_format=='cageBed9' else 'bed6'
        return display_fields(record, base)+[record['score']]
    f = list(record['fields'])
    score = '0' if record['score']=='.' else str(max(0, min(1000, int(float(record['score'])))))
    if region_format in ('bed9', 'polyAtlas'):
        f[0], f[5] = record['chromosome'], record['strand']
        if region_format == 'bed9':
            try:
                if not record['start'] <= int(f[6]) <= int(f[7]) <= record['end']:
                    raise ValueError()
                if f[8]!='0' and (not re.fullmatch(r'[0-9]+,[0-9]+,[0-9]+',f[8]) or
                                  any(int(v)>255 for v in f[8].split(','))):
                    raise ValueError()
            except ValueError:
                raise ValueError('Invalid BED9 thick bounds or RGB at source line '+str(record['line'])) from None
            f[4] = score
        else:
            try:
                if not all(math.isfinite(float(f[i])) for i in (4,6,8)) or not re.fullmatch(r'[0-9]+',f[7]):
                    raise ValueError()
            except ValueError:
                raise ValueError('Invalid PolyA Atlas signal/support at source line '+str(record['line'])) from None
        return f
    return [record['chromosome'], str(record['start']), str(record['end']), record['name'], score, record['strand']]


def prepare_sites(sources, scopes, output, reference_lengths=None):
    expected = scopes_by_chromosome(scopes)
    if not sources or set(sources)-{'cage', 'polya'}:
        raise ValueError('Supply CAGE and/or poly(A) evidence sources')
    if not shutil.which('sort'):
        raise ValueError('Site evidence requires local sort')
    output = Path(output).absolute()
    if output.exists() or output.is_symlink():
        raise ValueError('Site-evidence output must be a new directory')
    lengths = None
    if reference_lengths:
        helper = next(p for p in (Path(__file__).with_name('check_reference_coordinates.py'),
            Path(__file__).resolve().parents[1]/'assets/check_reference_coordinates.py') if p.is_file())
        spec = importlib.util.spec_from_file_location('review_site_reference', helper)
        checker = importlib.util.module_from_spec(spec); spec.loader.exec_module(checker)
        lengths = checker.reference_lengths(Path(reference_lengths))
        if set(expected)-set(lengths):
            raise ValueError('Selected site-evidence chromosomes are absent from the reference lengths')
    output.parent.mkdir(parents=True, exist_ok=True)
    summary = dict(scopes=sorted(scopes), sources={}, description='Original CAGE 5-prime and poly(A) 3-prime evidence intervals in the processed chromosome/strand scope; not restricted to sites selected for a transcript extension.',
        coordinates='Output BED uses zero-based starts and exclusive ends. GFF starts are converted once; BED starts are retained. A leading chr prefix is removed, matching LEAP evidence preparation; no other aliases or coordinate changes are inferred.')
    with tempfile.TemporaryDirectory(prefix='.'+output.name+'-', dir=output.parent) as temp:
        staging = Path(temp)/'sites'; staging.mkdir()
        for kind, source in sorted(sources.items()):
            source = Path(source)
            widths, formats = set(), set()
            counts = Counter(scanned=0, selected=0, excluded=0, unstranded=0, score_display_changes=0,
                             positive=0, negative=0, unsplit=0)
            for record in source_records(source):
                counts['scanned'] += 1
                if not selected(record, expected):
                    counts['excluded'] += 1; continue
                if lengths and record['end'] > lengths[record['chromosome']]:
                    raise ValueError('{}:{}: site exceeds reference length'.format(source, record['line']))
                widths.add(len(record['fields'])); formats.add(record['format'])
            preset = ('polyAtlas' if kind=='polya' and widths=={11} and formats=={'bed'} else
                      'bed9' if widths=={9} and formats=={'bed'} else 'bed6')
            if kind=='cage':
                preset = 'cageBed9' if preset=='bed9' else 'cageBed6'
                presets = staging/'region_presets';presets.mkdir(exist_ok=True)
                shutil.copyfile(cage_preset_path(preset),presets/(preset+'.tsv'))
            folder=staging/kind; folder.mkdir()
            handles={v:(folder/(v+'.unsorted.bed')).open('w', encoding='utf-8') for v in VIEWS}
            try:
                with (folder/'source_records.tsv').open('w', newline='', encoding='utf-8') as audit:
                    writer=csv.writer(audit, delimiter='\t');writer.writerow(['source_line','chromosome','bed_start','bed_end','source_fields_json'])
                    for record in source_records(source):
                        if not selected(record, expected):
                            continue
                        fields=display_fields(record, preset)
                        counts['selected'] += 1; counts['unsplit'] += 1
                        if preset!='polyAtlas' and fields[4]!=record['score']:
                            counts['score_display_changes'] += 1
                        line='\t'.join(fields)+'\n';handles['unsplit'].write(line)
                        for view, sign in STRANDS.items():
                            if record['strand']==sign:
                                handles[view].write(line);counts[view] += 1
                        counts['unstranded'] += record['strand']=='.'
                        writer.writerow([record['line'],record['chromosome'],record['start'],record['end'],json.dumps(record['fields'],ensure_ascii=False)])
            finally:
                for handle in handles.values():
                    handle.close()
            for view in VIEWS:
                unsorted=folder/(view+'.unsorted.bed')
                with (folder/(view+'.bed')).open('w', encoding='utf-8') as target:
                    subprocess.run(['sort','-s','-t','\t','-k1,1','-k2,2n','-k3,3n',str(unsorted)],
                        stdout=target, stderr=subprocess.PIPE, text=True, env=dict(os.environ,LC_ALL='C'),check=True)
                unsorted.unlink()
            summary['sources'][kind]=dict(input=str(source.resolve()), region_format=preset, counts=dict(counts),
                description=('CAGE clusters supporting 5-prime ends.' if kind=='cage' else 'Poly(A) sites supporting 3-prime ends.')+
                    ' Interval views contain all supplied evidence in scope. Original fields are retained in source_records.tsv. '+
                    ('PolyA Atlas signal/support columns are preserved unchanged.' if preset=='polyAtlas' else
                     'BED display scores are limited to 0..1000; original scores remain in the source-record table.'+
                     (' The browser-visible originalCageScore field preserves each supplied score exactly.' if kind=='cage' else '')))
        (staging/'site-summary.json').write_text(json.dumps(summary,indent=2)+'\n')
        staging.rename(output)
    return summary


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('manifest', type=Path);parser.add_argument('output', type=Path)
    parser.add_argument('--reference-lengths',type=Path)
    args=parser.parse_args()
    try:
        manifest=json.loads(args.manifest.read_text())
        prepare_sites(manifest['sources'],manifest['scopes'],args.output,args.reference_lengths)
        return 0
    except (OSError, ValueError, KeyError, EOFError, csv.Error, subprocess.SubprocessError) as error:
        parser.exit(2,'ERROR: '+str(error)+'\n')


if __name__=='__main__':
    raise SystemExit(main())
