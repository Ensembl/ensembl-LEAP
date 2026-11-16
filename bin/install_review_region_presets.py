#!/usr/bin/env python3
"""Install bundled CAGE field definitions into a JBrowser checkout without overwriting presets."""

import argparse
import csv
import os
from pathlib import Path
import tempfile

CAGE_PRESETS = {'cageBed6', 'cageBed9'}


def install_presets(review, jbrowser):
    review, jbrowser = Path(review), Path(jbrowser)
    if not (jbrowser/'main.nf').is_file():
        raise ValueError('JBrowser checkout must contain main.nf')
    with (review/'samplesheet.csv').open(newline='') as handle:
        required = {r.get('region_format') for r in csv.DictReader(handle)} & CAGE_PRESETS
    # Check all definitions and existing destinations before installing anything.
    definitions = {}
    target = jbrowser/'assets/region_presets'
    for name in sorted(required):
        source = review/'region_presets'/(name+'.tsv')
        content = source.read_bytes()
        destination = target/source.name
        if destination.is_symlink() or (destination.exists() and destination.read_bytes()!=content):
            raise ValueError('Conflicting JBrowser preset; review it before replacing: '+str(destination))
        definitions[destination] = content
    if definitions:
        target.mkdir(parents=True, exist_ok=True)
    for destination, content in definitions.items():
        if destination.exists():
            continue
        with tempfile.NamedTemporaryFile(dir=target, prefix='.leap-cage-', delete=False) as handle:
            temporary = Path(handle.name)
            handle.write(content)
        try:
            temporary.chmod(0o644)
            # A complete file becomes visible atomically; an existing preset is never replaced.
            try:
                os.link(temporary, destination)
            except FileExistsError:
                if destination.is_symlink() or destination.read_bytes()!=content:
                    raise ValueError('Conflicting JBrowser preset: '+str(destination)) from None
        finally:
            temporary.unlink()
    return sorted(required)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('review_dir', type=Path)
    parser.add_argument('jbrowser_dir', type=Path)
    args = parser.parse_args()
    try:
        installed = install_presets(args.review_dir, args.jbrowser_dir)
    except (OSError, ValueError, csv.Error) as error:
        parser.exit(2, 'ERROR: '+str(error)+'\n')
    if installed:
        print('JBrowser CAGE field definitions ready: '+', '.join(installed))


if __name__=='__main__':
    main()
