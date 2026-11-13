from __future__ import annotations

import gzip
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable
from urllib.parse import quote, unquote


@dataclass
class Feature:
    chrom: str
    source: str
    feature_type: str
    start: int
    end: int
    score: str
    strand: str
    phase: str
    attrs: dict[str, str | list[str]] = field(default_factory=dict)


GFF3_ATTR_ORDER = [
    "ID",
    "Parent",
    "gene_id",
    "transcript_id",
    "gene_type",
    "gene_name",
    "transcript_type",
    "transcript_name",
    "exon_number",
    "exon_id",
    "tag",
    "protein_id",
    "Dbxref",
]

GTF_ATTR_ORDER = [
    "gene_id",
    "transcript_id",
    "gene_type",
    "gene_name",
    "transcript_type",
    "transcript_name",
    "exon_number",
    "exon_id",
    "tag",
    "protein_id",
    "db_xref",
]


def open_text(path: str | Path):
    path = Path(path)
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def parse_gtf_attrs(raw: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    pattern = re.compile(r'\s*([^;\s]+)(?:\s+(?:"((?:\\.|[^"])*)"|([^;]*)))?\s*;?')
    pos = 0
    text = raw.strip()
    while pos < len(text):
        match = pattern.match(text, pos)
        if not match or not match.group(1):
            break
        key = match.group(1)
        quoted_value = match.group(2)
        bare_value = match.group(3)
        value = quoted_value if quoted_value is not None else (bare_value.strip() if bare_value is not None else "")
        value = value.replace(r'\"', '"')
        attrs[key] = f"{attrs[key]},{value}" if key in attrs else value
        pos = match.end()
    return attrs


def parse_gff3_attrs(raw: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for item in raw.strip().split(";"):
        if not item:
            continue
        if "=" not in item:
            attrs[unquote(item)] = ""
            continue
        key, value = item.split("=", 1)
        attrs[unquote(key)] = unquote(value)
    return attrs


def ordered_attr_items(attrs: dict[str, object], preferred: list[str]) -> Iterable[tuple[str, object]]:
    seen = set()
    for key in preferred:
        if key in attrs and attrs[key] not in (None, ""):
            seen.add(key)
            yield key, attrs[key]
    for key, value in attrs.items():
        if key not in seen and value not in (None, ""):
            yield key, value


def format_gtf_attrs(attrs: dict[str, object]) -> str:
    fields = []
    for key, value in ordered_attr_items(attrs, GTF_ATTR_ORDER):
        if key == "exon_number":
            fields.append(f"{key} {value};")
        elif key == "tag":
            fields.extend(f'{key} "{part}";' for part in str(value).split(",") if part)
        elif key == "hidden_remark" and isinstance(value, list):
            for remark in value:
                escaped = (remark.replace('\\', '\\\\').replace('"', '\\"')
                           .replace('\t', '\\t').replace('\r', '\\r').replace('\n', '\\n'))
                fields.append(f'{key} "{escaped}";')
        else:
            fields.append(f'{key} "{value}";')
    return " ".join(fields)


def escape_gff3_remark(value: str) -> str:
    return re.sub(r'[%;&=,\x00-\x1f\x7f]', lambda match: f'%{ord(match[0]):02X}', value)


def format_gff3_attrs(attrs: dict[str, object]) -> str:
    fields = []
    for key, value in ordered_attr_items(attrs, GFF3_ATTR_ORDER):
        # Live remarks are lists; encode each paragraph before adding list separators.
        if key == "hidden_remark" and isinstance(value, list):
            encoded = ','.join(escape_gff3_remark(remark) for remark in value)
        else:
            encoded = quote(str(value), safe=':,.|_-')
        fields.append(f"{quote(str(key), safe='')}={encoded}")
    return ';'.join(fields)


def format_gff3_attrs_perl(attrs: dict[str, object], feature_type: str) -> str:
    text = format_gff3_attrs(attrs)
    ordered = dict(ordered_attr_items(attrs, GFF3_ATTR_ORDER))
    if feature_type != "gene" and text and next(reversed(ordered)) not in {"tag", "protein_id"}:
        return f"{text};"
    return text


def parse_feature(line: str, fmt: str | None = None) -> Feature:
    fields = line.rstrip("\n").split("\t")
    if len(fields) != 9:
        raise ValueError(f"Expected 9 tab-delimited fields, found {len(fields)}: {line[:120]}")
    # GTF separates a key and value with whitespace; GFF3 uses '=' even with quoted text.
    inferred = fmt or ("gtf" if re.match(r'\s*[^=;\s]+\s+', fields[8]) else "gff3")
    attrs = parse_gtf_attrs(fields[8]) if inferred == "gtf" else parse_gff3_attrs(fields[8])
    return Feature(
        chrom=fields[0],
        source=fields[1],
        feature_type=fields[2],
        start=int(fields[3]),
        end=int(fields[4]),
        score=fields[5],
        strand=fields[6],
        phase=fields[7],
        attrs=attrs,
    )


def read_features(path: str | Path, fmt: str | None = None) -> Iterable[Feature]:
    with open_text(path) as handle:
        for line in handle:
            if line.strip() and not line.startswith("#"):
                yield parse_feature(line, fmt)


def fail(message: str, code: int = 2) -> None:
    print(f"ERROR: {message}", file=sys.stderr)
    raise SystemExit(code)
