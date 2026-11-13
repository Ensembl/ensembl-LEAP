"""Annotation-region selection, independent of DNA assembly mapping routes."""
from __future__ import annotations


def parse_chromosomes(value):
    """Normalize explicit region names without changing their case or spelling."""
    if value is None:
        return None
    names = value.split(',') if isinstance(value, str) else value
    names = tuple(dict.fromkeys(name.strip() for name in names))
    if not names or any(not name for name in names):
        raise ValueError('Chromosome lists must contain non-empty region names')
    return names


def validate_region_selection(chromosome=None, region_set=None):
    names = parse_chromosomes(chromosome)
    if names is not None and region_set is not None:
        raise ValueError('Use either --chromosome or --region-set, not both')
    if region_set is not None and region_set not in {'CHR', 'PRI', 'ALL'}:
        raise ValueError('Region set must be CHR, PRI or ALL')
    return names


def validate_requested_regions(names, available, context):
    missing = set(names or ()) - set(available)
    if missing:
        raise ValueError(f"No region(s) {', '.join(sorted(missing))} in {context}")


def validate_slice_selection(region_count, start, end):
    if start is None and end is None:
        return
    if region_count != 1:
        raise ValueError('--start/--end require selection of exactly one region')
    if (start is not None and start < 1) or (end is not None and end < (start or 1)):
        raise ValueError('Invalid slice bounds: require 1 <= start <= end')


def select_regions(systems, defaults, regions, chromosome=None, region_set=None):
    requested = validate_region_selection(chromosome, region_set)
    preset = region_set or 'CHR'
    species = {row['species_id'] for row in systems}
    if len(species) != 1:
        raise ValueError('Annotation selection requires one unambiguous species')
    species_id = next(iter(species))
    assemblies = {row['meta_value'] for row in defaults if row['species_id'] == species_id}
    if len(assemblies) != 1:
        raise ValueError('Missing or ambiguous assembly.default for annotation selection')
    assembly = next(iter(assemblies))
    selected_systems = {}
    for row in systems:
        attributes = set((row['attrib'] or '').split(','))
        if 'default_version' in attributes and row['version'] in (None, '', assembly):
            selected_systems[row['coord_system_id']] = row
    if not any(row['version'] == assembly for row in selected_systems.values()):
        raise ValueError(f'No default coordinate system for assembly {assembly}')
    selected = []
    for row in regions:
        system = selected_systems.get(row['coord_system_id'])
        if system is None or not row['toplevel']:
            continue
        if requested is not None:
            if row['name'] not in requested:
                continue
        else:
            if preset != 'ALL' and row['non_ref']:
                continue
            if preset == 'CHR' and system['name'] != 'chromosome':
                continue
        selected.append(dict(row, coordinate_system=system['name']))
    names = [row['name'] for row in selected]
    if len(names) != len(set(names)):
        raise ValueError(f'Ambiguous sequence-region names in default assembly {assembly}')
    validate_requested_regions(requested, names, f'default assembly {assembly}')
    return sorted(selected, key=lambda row: row['name'])


def fetch_regions(connection, chromosome=None, region_set=None):
    from sqlalchemy import text

    def rows(sql):
        return [dict(row._mapping) for row in connection.execute(text(sql))]

    systems = rows('SELECT coord_system_id, species_id, name, version, attrib FROM coord_system')
    defaults = rows("SELECT species_id, meta_value FROM meta WHERE meta_key = 'assembly.default'")
    regions = rows('''
        SELECT sr.seq_region_id, sr.name, sr.length, sr.coord_system_id,
            EXISTS (SELECT 1 FROM seq_region_attrib sa JOIN attrib_type a USING (attrib_type_id)
                    WHERE sa.seq_region_id = sr.seq_region_id AND a.code = 'toplevel') AS toplevel,
            EXISTS (SELECT 1 FROM seq_region_attrib sa JOIN attrib_type a USING (attrib_type_id)
                    WHERE sa.seq_region_id = sr.seq_region_id AND a.code = 'non_ref') AS non_ref
        FROM seq_region sr
        WHERE EXISTS (SELECT 1 FROM seq_region_attrib sa JOIN attrib_type a USING (attrib_type_id)
                      WHERE sa.seq_region_id = sr.seq_region_id AND a.code = 'toplevel')
    ''')
    return select_regions(systems, defaults, regions, chromosome, region_set)
