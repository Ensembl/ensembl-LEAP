"""Interval DNA retrieval with Ensembl assembly-exception normalization."""

from dataclasses import dataclass
from functools import lru_cache

from .mapping import MappingError, SequenceMappingPaths


class SequenceError(ValueError):
    """A requested interval cannot be reconstructed reliably."""


class _SequenceConflictError(SequenceError):
    """Components on the selected path disagree at a mapped position."""


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGTRYMKBDHVN", "TGCAYRKMVHDBN"))[::-1]


@dataclass(frozen=True)
class _ResolvedSequence:
    sequence: str
    # A DNA-backed N is mapped; an assembly-gap N is not. Preserve this through
    # recursive projections so padding never creates a component conflict.
    mapped: bytes

    @classmethod
    def gap(cls, length):
        return cls('N' * length, bytes(length))

    @classmethod
    def join(cls, chunks):
        return cls(''.join(chunk.sequence for chunk in chunks), b''.join(chunk.mapped for chunk in chunks))

    def reversed(self):
        return type(self)(reverse_complement(self.sequence), self.mapped[::-1])


class SQLSequenceStore:
    def __init__(self, connection):
        self.connection = connection
        # Caches belong to this connection/export, not a global connection pool.
        self.region = lru_cache(maxsize=4096)(self.region)
        self.exceptions = lru_cache(maxsize=4096)(self.exceptions)
        self.assembly = lru_cache(maxsize=1024)(self.assembly)
        self.mapping_paths = lru_cache(maxsize=16)(self.mapping_paths)

    def rows(self, sql, **params):
        from sqlalchemy import text
        return [dict(row) for row in self.connection.execute(text(sql), params).mappings()]

    def region(self, region_id):
        rows = self.rows("""SELECT sr.name, sr.length, sr.coord_system_id,
            cs.attrib, cs.species_id FROM seq_region sr JOIN coord_system cs
            ON cs.coord_system_id = sr.coord_system_id WHERE sr.seq_region_id = :id""", id=region_id)
        if not rows:
            raise SequenceError(f"Unknown seq_region_id {region_id}")
        return rows[0]

    def exceptions(self, region_id):
        return self.rows("""SELECT seq_region_start AS start, seq_region_end AS end,
            exc_type AS type, exc_seq_region_id AS target,
            exc_seq_region_start AS target_start, exc_seq_region_end AS target_end, ori
            FROM assembly_exception WHERE seq_region_id = :id
            ORDER BY seq_region_start, seq_region_end""", id=region_id)

    def mapping_paths(self, species_id):
        systems = self.rows("""SELECT coord_system_id, species_id, name, version, attrib
            FROM coord_system WHERE species_id = :species""", species=species_id)
        declarations = self.rows("""SELECT meta_value FROM meta
            WHERE species_id = :species AND meta_key = 'assembly.mapping'
            ORDER BY meta_id""", species=species_id)
        return SequenceMappingPaths(systems, [row['meta_value'] for row in declarations])

    def assembly(self, region_id, target_cs):
        forward = self.rows("""SELECT a.asm_start AS start, a.asm_end AS end,
            a.cmp_seq_region_id AS target, a.cmp_start AS target_start,
            a.cmp_end AS target_end, a.ori, sr.coord_system_id, cs.attrib
            FROM assembly a JOIN seq_region sr ON sr.seq_region_id = a.cmp_seq_region_id
            JOIN coord_system cs ON cs.coord_system_id = sr.coord_system_id
            WHERE a.asm_seq_region_id = :id AND sr.coord_system_id = :cs
            ORDER BY a.asm_start, a.asm_end""", id=region_id, cs=target_cs)
        reverse = self.rows("""SELECT a.cmp_start AS start, a.cmp_end AS end,
            a.asm_seq_region_id AS target, a.asm_start AS target_start,
            a.asm_end AS target_end, a.ori, sr.coord_system_id, cs.attrib
            FROM assembly a JOIN seq_region sr ON sr.seq_region_id = a.asm_seq_region_id
            JOIN coord_system cs ON cs.coord_system_id = sr.coord_system_id
            WHERE a.cmp_seq_region_id = :id AND sr.coord_system_id = :cs
            ORDER BY a.cmp_start, a.cmp_end""", id=region_id, cs=target_cs)
        return sorted(forward + reverse, key=lambda row: (row['start'], row['end'], row['target']))

    def dna(self, region_id, start, end):
        rows = self.rows("""SELECT SUBSTR(sequence, :start, :length) AS sequence
            FROM dna WHERE seq_region_id = :id""", id=region_id, start=start, length=end-start+1)
        return None if not rows or rows[0]["sequence"] is None else rows[0]["sequence"].upper()


class SequenceResolver:
    def __init__(self, store):
        self.store = store

    def fetch(self, region_id, start, end, path=()):
        return self._fetch(region_id, start, end, path).sequence

    def _fetch(self, region_id, start, end, path, steps=None):
        region_id = int(region_id)
        region = self.store.region(region_id)
        context = f"{region['name']} ({region_id}):{start}-{end}"
        if region_id in path or len(path) >= 64:
            raise SequenceError(f"Projection cycle/depth limit at {context}; route {path}")
        if start > end:
            raise SequenceError(f"Invalid interval {context}")
        try:
            # Ensembl pads slice overhangs, including unequal terminal patch flanks.
            lo, hi = max(1, start), min(region['length'], end)
            if lo > hi:
                return _ResolvedSequence.gap(end-start+1)
            sequence = _ResolvedSequence.join([
                _ResolvedSequence.gap(lo-start),
                self._normalized(region_id, lo, hi, path + (region_id,), steps),
                _ResolvedSequence.gap(end-hi),
            ])
            if len(sequence.sequence) != end-start+1:
                raise SequenceError("Reconstructed sequence has incorrect length")
            return sequence
        except SequenceError as error:
            raise type(error)(f"{context}: {error}") from error

    def _normalized(self, region_id, start, end, path, steps=None):
        exceptions = self.store.exceptions(region_id)
        projections = [row for row in exceptions if row['type'] == 'PAR']
        patches = [row for row in exceptions if row['type'] != 'PAR']
        if patches:
            targets = {row['target'] for row in patches}
            if len(targets) != 1 or any(row['ori'] != 1 for row in patches):
                raise SequenceError("Ambiguous/non-forward patch exception mapping")
            source_pos = target_pos = 1
            target = patches[0]['target']
            for row in patches:
                if row['start'] < source_pos or row['target_start'] < target_pos:
                    raise SequenceError("Overlapping or out-of-order patch exceptions")
                if row['start']-source_pos != row['target_start']-target_pos:
                    raise SequenceError("Patch flanks have inconsistent lengths")
                if source_pos < row['start']:
                    projections.append(dict(start=source_pos, end=row['start']-1,
                        target=target, target_start=target_pos, target_end=row['target_start']-1, ori=1))
                source_pos, target_pos = row['end']+1, row['target_end']+1
            length = self.store.region(region_id)['length']
            if source_pos <= length:
                projections.append(dict(start=source_pos, end=length, target=target,
                    target_start=target_pos, target_end=target_pos+length-source_pos, ori=1))
        chunks, cursor = [], start
        last_end = 0
        for row in sorted(projections, key=lambda item: item['start']):
            self._validate_mapping(row)
            if row['start'] <= last_end:
                raise SequenceError("Overlapping exception projections")
            last_end = row['end']
            lo, hi = max(start, row['start']), min(end, row['end'])
            if lo > hi:
                continue
            if lo > cursor:
                chunks.append(self._assembled(region_id, cursor, lo-1, path, steps))
            chunks.append(self._project(row, lo, hi, path))
            cursor = hi+1
        if cursor <= end:
            chunks.append(self._assembled(region_id, cursor, end, path, steps))
        return _ResolvedSequence.join(chunks)

    @staticmethod
    def _validate_mapping(row):
        if (row['ori'] not in (-1, 1) or row['start'] < 1 or row['end'] < row['start']
                or row['target_start'] < 1 or row['end']-row['start'] != row['target_end']-row['target_start']):
            raise SequenceError(f"Invalid projection {row}")

    def _project(self, row, start, end, path, steps=None):
        self._validate_mapping(row)
        if row['ori'] == 1:
            lo = row['target_start'] + start-row['start']
            hi = row['target_start'] + end-row['start']
        else:
            lo = row['target_end'] - (end-row['start'])
            hi = row['target_end'] - (start-row['start'])
        sequence = self._fetch(row['target'], lo, hi, path, steps)
        return sequence if row['ori'] == 1 else sequence.reversed()

    @staticmethod
    def _route_description(coord_system_id, rows, start, end):
        mappings = []
        for row in rows:
            lo, hi = max(start, row['start']), min(end, row['end'])
            if row['ori'] == 1:
                target_lo = row['target_start'] + lo-row['start']
                target_hi = row['target_start'] + hi-row['start']
            else:
                target_lo = row['target_end'] - (hi-row['start'])
                target_hi = row['target_end'] - (lo-row['start'])
            mappings.append(f"{lo}-{hi} -> component {row['target']}:{target_lo}-{target_hi} (ori={row['ori']})")
        return f"coord_system_id={coord_system_id} [" + '; '.join(mappings) + ']'

    def _assembled(self, region_id, start, end, path, steps=None):
        region = self.store.region(region_id)
        try:
            routing = self.store.mapping_paths(region['species_id'])
            if steps is None:
                steps = routing.to_sequence(region['coord_system_id'])
        except MappingError as error:
            raise SequenceError(str(error)) from error
        description = routing.describe(steps)
        try:
            if not steps:
                if region['coord_system_id'] != routing.sequence_level:
                    raise SequenceError('Selected path did not reach sequence_level')
                dna = self.store.dna(region_id, start, end)
                if dna is None:
                    raise SequenceError('No DNA available at sequence_level')
                if len(dna) != end-start+1:
                    raise SequenceError('Truncated DNA sequence')
                return _ResolvedSequence(dna, b'\x01' * len(dna))
            step = steps[0]
            if step.source != region['coord_system_id']:
                raise SequenceError('Selected path does not match current coordinate system')
            components = self.store.assembly(region_id, step.target)
            if not components:
                raise SequenceError('No assembly components available on selected path')
            rows = [row for row in components
                    if row['start'] <= end and row['end'] >= start]
            if not rows:
                # With an authoritative route, an interval between its components
                # is a genuine gap. Missing DNA in a mapped component still fails.
                return _ResolvedSequence.gap(end-start+1)
            sequence = bytearray(b'N' * (end-start+1))
            mapped = bytearray(end-start+1)
            last_end = 0
            previous = []
            for row in rows:
                self._validate_mapping(row)
                if row['start'] <= last_end and not step.multiple:
                    raise SequenceError('Overlapping components on a single-component | mapping')
                last_end = max(last_end, row['end'])
                if row['target_end'] > self.store.region(row['target'])['length']:
                    raise SequenceError('Assembly component mapping exceeds target region length')
                lo, hi = max(start, row['start']), min(end, row['end'])
                piece = self._project(row, lo, hi, path, steps[1:])
                data = piece.sequence.encode('ascii')
                offset = lo-start
                for i, covered in enumerate(piece.mapped):
                    if not covered:
                        continue
                    pos = offset+i
                    if mapped[pos] and sequence[pos] != data[i]:
                        mappings = self._route_description(step.target, previous + [row], start, end)
                        raise _SequenceConflictError(
                            f'Conflicting components on selected path; first conflict at '
                            f'seq_region_id={region_id}:{start+pos} '
                            f'({chr(sequence[pos])} vs {chr(data[i])}); {mappings}'
                        )
                    sequence[pos] = data[i]
                    mapped[pos] = 1
                previous.append(row)
            return _ResolvedSequence(sequence.decode('ascii'), bytes(mapped))
        except SequenceError as error:
            raise type(error)(f'Selected path {description}: {error}') from error
