"""Select declared Ensembl coordinate-system paths, not arbitrary assembly edges."""

import re
from dataclasses import dataclass


class MappingError(ValueError):
    pass


@dataclass(frozen=True)
class MappingStep:
    source: int
    target: int
    multiple: bool = False

    def reversed(self):
        return MappingStep(self.target, self.source, self.multiple)


class SequenceMappingPaths:
    def __init__(self, coordinate_systems, declarations):
        self.systems = {row['coord_system_id']: row for row in coordinate_systems}
        sequence_levels = [id for id, row in self.systems.items()
                           if 'sequence_level' in (row['attrib'] or '').split(',')]
        if len(sequence_levels) != 1:
            raise MappingError(f'Expected one sequence_level coordinate system; found {sequence_levels}')
        self.sequence_level = sequence_levels[0]
        self.paths = {}
        self.edges = {}
        for declaration in declarations:
            parts = re.split(r'([|#])', declaration)
            if len(parts) < 3 or any(not part.strip() for part in parts):
                raise MappingError(f'Invalid assembly.mapping: {declaration!r}')
            ids = [self._resolve(name.strip()) for name in parts[::2]]
            if len(ids) != len(set(ids)):
                raise MappingError(f'Cyclic assembly.mapping: {declaration!r}')
            steps = tuple(MappingStep(a, b, separator == '#')
                          for a, b, separator in zip(ids, ids[1:], parts[1::2]))
            reverse = tuple(step.reversed() for step in reversed(steps))
            self.paths.setdefault((ids[0], ids[-1]), set()).add(steps)
            self.paths.setdefault((ids[-1], ids[0]), set()).add(reverse)
            if len(steps) == 1:
                self.edges.setdefault(ids[0], {}).setdefault(ids[-1], set()).add(steps[0])
                self.edges.setdefault(ids[-1], {}).setdefault(ids[0], set()).add(reverse[0])

    def label(self, id):
        row = self.systems[id]
        version = f":{row['version']}" if row['version'] else ''
        return f"{row['name']}{version} (coord_system_id={id})"

    def describe(self, steps):
        if not steps:
            return self.label(self.sequence_level)
        return self.label(steps[0].source) + ''.join(
            (' # ' if step.multiple else ' | ') + self.label(step.target) for step in steps)

    def _resolve(self, name):
        fields = name.split(':', 1)
        candidates = [id for id, row in self.systems.items() if row['name'].lower() == fields[0].lower()]
        if len(fields) == 2 and fields[1]:
            candidates = [id for id in candidates
                          if (self.systems[id]['version'] or '').lower() == fields[1].lower()]
        else:
            defaults = [id for id in candidates
                        if 'default_version' in (self.systems[id]['attrib'] or '').split(',')]
            if defaults:
                candidates = defaults
        if len(candidates) != 1:
            raise MappingError(f'Unknown or ambiguous coordinate system {name!r}: {candidates}')
        return candidates[0]

    def _choose(self, paths, source):
        # Ensembl represents a two-system '#' declaration as a longer mapper
        # path, unlike an ordinary direct '|' declaration.
        def cost(path):
            return len(path) + int(len(path) == 1 and path[0].multiple)

        shortest = min(map(cost, paths))
        choices = {path for path in paths if cost(path) == shortest}
        if len(choices) != 1:
            descriptions = sorted(self.describe(path) for path in choices)
            raise MappingError(f'Ambiguous sequence path from {self.label(source)}: ' + '; '.join(descriptions))
        return next(iter(choices))

    def to_sequence(self, source):
        if source not in self.systems:
            raise MappingError(f'Coordinate system {source} is outside the selected species')
        if source == self.sequence_level:
            return ()
        explicit = self.paths.get((source, self.sequence_level))
        if explicit:
            return self._choose(explicit, source)
        inferred = set()
        for middle, first_steps in self.edges.get(source, {}).items():
            if middle in (source, self.sequence_level):
                continue
            for first in first_steps:
                for second in self.edges.get(middle, {}).get(self.sequence_level, ()):
                    inferred.add((first, second))
        if not inferred:
            raise MappingError(f'No declared sequence path from {self.label(source)} to {self.label(self.sequence_level)}')
        return self._choose(inferred, source)
