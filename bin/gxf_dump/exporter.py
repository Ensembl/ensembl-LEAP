from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

from .annotation import Feature, fail, format_gff3_attrs_perl, format_gtf_attrs, open_text, parse_feature
from .regions import parse_chromosomes, validate_region_selection, validate_requested_regions, validate_slice_selection


# LEAP additions to the pinned upstream exporter: expose curated MANE remarks
# and normalize the known selection tags without exporting free-text remarks.
SELECTION_REMARK_TAGS = {
    'mane_select': 'MANE_Select',
    'mane_plus_clinical': 'MANE_Plus_Clinical',
    'ensembl_canonical_extended': 'Ensembl_canonical_extended',
    'upstream atg': 'upstream_ATG',
    'cage supported tss': 'CAGE_supported_TSS',
    'rp supported tis': 'RP_supported_TIS',
}


CODING_BIOTYPES = {
    "protein_coding",
    "nonsense_mediated_decay",
    "non_stop_decay",
    "polymorphic_pseudogene",
    "protein_coding_LoF",
}

REMARK_TAGS = {
    "alternative 3' UTR", "alternative 5' UTR", "bicistronic", "CAGE supported TSS",
    "dotter confirmed", "downstream ATG", "Ensembl_canonical_extended",
    "inferred exon combination", "inferred transcript model", "low sequence quality",
    "NAGNAG splice site", "NMD exception", "NMD likely if extended", "non-ATG start",
    "non canonical conserved", "non canonical genome sequence error", "non canonical other",
    "non canonical polymorphism", "non canonical TEC", "non canonical U12",
    "non-submitted evidence", "not best-in-genome evidence", "not for VEGA",
    "not organism-supported", "overlapping uORF", "retained intron CDS",
    "retained intron final", "retained intron first", "RNA-Seq supported only",
    "RNA-Seq supported partial", "RP supported TIS", "sequence error", "upstream ATG",
    "upstream uORF", "3' nested supported extension", "3' standard supported extension",
    "5' nested supported extension", "5' standard supported extension",
    "454 RNA-Seq supported", "nested 454 RNA-Seq supported",
}


def perl_tag(value: str) -> str:
    if re.fullmatch(r'[Nn]o.+org.+sup.+ted', value):
        return 'not_organism_supported'
    return value.replace("'", "").replace("-", "_").replace(" ", "_")


def annotation_name(row: dict, kind: str, stable_id: str,
                    attributes: Iterable[tuple[str, str]] = ()) -> str:
    name = row.get(f'{kind}_display_name') or row.get(f'{kind}_external_name')
    if name:
        return name
    # Perl uses get_name for genes, but _get_name (without this fallback) for transcripts.
    if kind == 'gene':
        for code, value in attributes:
            if code == 'name':
                return value
    return stable_id


def unique_preserve_order(values: list[str]) -> list[str]:
    seen = set()
    out = []
    for value in values:
        if value in seen:
            continue
        seen.add(value)
        out.append(value)
    return out


@dataclass
class ExonRecord:
    exon_db_id: str
    exon_id: str
    start: int
    end: int
    rank: int
    seq_region_id: str | None = None
    version: Any = None
    phase: str = "."
    end_phase: str = "."
    sequence: str | None = None


@dataclass
class TranscriptRecord:
    chrom: str
    source: str
    gene_db_id: int
    gene_id: str
    gene_stable_id: str
    gene_version: Any
    gene_name: str
    gene_biotype: str
    gene_start: int
    gene_end: int
    gene_strand: str
    transcript_id: str
    transcript_db_id: int
    transcript_stable_id: str
    transcript_version: Any
    transcript_name: str
    transcript_biotype: str
    transcript_start: int
    transcript_end: int
    transcript_strand: str
    transcript_source: str
    protein_id: str | None = None
    protein_stable_id: str | None = None
    protein_version: Any = None
    translation_start_exon_id: str | None = None
    translation_end_exon_id: str | None = None
    seq_start: int | None = None
    seq_end: int | None = None
    peptide_sequence: str | None = None
    tags: list[str] = field(default_factory=list)
    gene_tags: list[str] = field(default_factory=list)
    transcript_remarks: list[str] = field(default_factory=list)
    transcript_attributes: list[tuple[str, str]] = field(default_factory=list)
    cds_start_nf: bool = False
    cds_end_nf: bool = False
    mrna_start_nf: bool = False
    mrna_end_nf: bool = False
    split_transcript: bool = False
    exons: list[ExonRecord] = field(default_factory=list)


def normalise_secrets(raw: dict[str, Any]) -> dict[str, Any]:
    database = raw.get("database") or raw.get("dbname")
    password = raw.get("password") if "password" in raw else raw.get("pass", "")
    required = ["host", "port", "user"]
    missing = [key for key in required if raw.get(key) in (None, "")]
    if not database:
        missing.append("database/dbname")
    if missing:
        fail(f"DB secrets JSON is missing required field(s): {', '.join(missing)}")
    return {
        "host": raw["host"],
        "port": int(raw.get("port", 3306)),
        "user": raw["user"],
        "password": password or "",
        "database": database,
        "dialect": raw.get("dialect", "mysql"),
    }


def load_secrets(path: str | Path) -> dict[str, Any]:
    with Path(path).open("r", encoding="utf-8") as handle:
        return normalise_secrets(json.load(handle))


def sql_connection(secrets: dict[str, Any]):
    try:
        from sqlalchemy import create_engine
        from sqlalchemy.engine import URL
    except ImportError as exc:
        fail("SQLAlchemy is required for live DB export. Install sqlalchemy and a MySQL driver.")
        raise exc

    driver = "mysql+pymysql" if secrets["dialect"] == "mysql" else secrets["dialect"]
    url = URL.create(
        drivername=driver,
        username=secrets["user"],
        password=secrets["password"],
        host=secrets["host"],
        port=secrets["port"],
        database=secrets["database"],
    )
    return create_engine(url).connect()


def stable_with_version(stable_id: Any, version: Any, no_version: bool) -> str:
    if stable_id in (None, ""):
        return ""
    if no_version or version in (None, ""):
        return str(stable_id)
    return f"{stable_id}.{version}"


def strand_symbol(value: Any) -> str:
    if str(value) in {"1", "+", "1.0"}:
        return "+"
    if str(value) in {"-1", "-", "-1.0"}:
        return "-"
    return "."


def reverse_complement(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def orient_exon_sequence(seq: Any, strand: str) -> str | None:
    if seq in (None, ""):
        return None
    text = str(seq).upper()
    return reverse_complement(text).upper() if strand == "-" else text


def read_id_list(path: str | None) -> set[str] | None:
    if not path:
        return None
    with Path(path).open("r", encoding="utf-8") as handle:
        return {line.strip().split()[0] for line in handle if line.strip()}


def load_mane_map(path: str | None) -> dict[str, str]:
    if not path:
        return {}
    mapping: dict[str, str] = {}
    with Path(path).open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 9:
                rs_tr, rs_prot, ens_tr, ens_prot = fields[5], fields[6], fields[7], fields[8]
                if ens_tr:
                    mapping[ens_tr] = rs_tr
                if ens_prot:
                    mapping[ens_prot] = rs_prot
    return mapping


def fetch_chromosomes_from_db(secrets: dict[str, Any], include_chromosomes: bool, include_scaffolds: bool, include_patches: bool) -> list[str]:
    from .regions import fetch_regions

    if not any([include_chromosomes, include_scaffolds, include_patches]):
        fail("At least one of include_chromosomes/include_scaffolds/include_patches must be true")
    with sql_connection(secrets) as conn:
        regions = fetch_regions(conn, region_set='ALL' if include_patches else 'PRI')
    return [row['name'] for row in regions if
            (bool(row['non_ref']) and include_patches) or
            (not row['non_ref'] and
             ((row['coordinate_system'] == 'chromosome' and include_chromosomes) or
              (row['coordinate_system'] != 'chromosome' and include_scaffolds)))]


def fetch_annotation_rows(
    secrets: dict[str, Any],
    chromosome: str | tuple[str, ...] | None = None,
    start: int | None = None,
    end: int | None = None,
    transcript_ids: set[str] | None = None,
    gene_ids: set[str] | None = None,
    biotypes: set[str] | None = None,
    source: str | None = None,
    analysis: str | None = None,
    region_set: str | None = None,
    allow_artifacts: bool = False,
    vega_only: bool = False,
    regions: list[dict] | None = None,
) -> list[dict[str, Any]]:
    """Fetch an Ensembl-core structural annotation subset for db2gxf parity.

    The query intentionally keeps one row per transcript exon and includes
    translation metadata. Additional attributes are fetched separately.
    """
    from sqlalchemy import bindparam, text
    from .regions import fetch_regions

    if regions is None:
        with sql_connection(secrets) as conn:
            regions = fetch_regions(conn, chromosome, region_set)

    clauses = ["g.is_current = 1", "g.seq_region_id IN :region_ids"]
    transcript_clauses = ['t.gene_id = g.gene_id']
    params: dict[str, Any] = {'region_ids': tuple(row['seq_region_id'] for row in regions)}
    if start is not None:
        clauses.append("g.seq_region_end >= :start")
        params["start"] = start
    if end is not None:
        clauses.append("g.seq_region_start <= :end")
        params["end"] = end
    if transcript_ids:
        clauses.append("EXISTS (SELECT 1 FROM transcript selected WHERE selected.gene_id = g.gene_id AND selected.stable_id IN :transcript_ids)")
        transcript_clauses.append("t.stable_id IN :transcript_ids")
        params["transcript_ids"] = tuple(transcript_ids)
    if gene_ids:
        clauses.append("g.stable_id IN :gene_ids")
        params["gene_ids"] = tuple(gene_ids)
    if biotypes:
        clauses.append("(g.biotype IN :biotypes OR t.biotype IN :biotypes)")
        params["biotypes"] = tuple(biotypes)
    if source:
        clauses.append("(g.source = :source OR t.source = :source)")
        params["source"] = source
    if analysis:
        clauses.append("ga.logic_name = :analysis")
        params["analysis"] = analysis
    if not allow_artifacts:
        duplication = """EXISTS (SELECT 1 FROM gene_attrib dup JOIN attrib_type dat USING (attrib_type_id)
                                  WHERE dup.gene_id = g.gene_id AND dat.code = 'artef_dupl')"""
        clauses.append(f"(g.biotype != 'artifact' OR {duplication})")
        transcript_clauses.append(f"(t.biotype != 'artifact' OR {duplication})")
    if vega_only:
        transcript_clauses.append(
            """
            NOT EXISTS (
                SELECT 1
                FROM transcript_attrib ta
                JOIN attrib_type tat ON tat.attrib_type_id = ta.attrib_type_id
                WHERE ta.transcript_id = t.transcript_id
                  AND tat.code = 'remark'
                  AND ta.value = 'not for VEGA'
            )
            """
        )
        clauses.append(
            """
            NOT EXISTS (
                SELECT 1
                FROM gene_attrib gat
                JOIN attrib_type gaat ON gaat.attrib_type_id = gat.attrib_type_id
                WHERE gat.gene_id = g.gene_id
                  AND gaat.code = 'remark'
                  AND gat.value = 'not for VEGA'
            )
            """
        )
    where = " AND ".join(clauses)
    transcript_join = " AND ".join(transcript_clauses)

    sql = text(
        f"""
        SELECT
            sr.name AS chrom,
            g.gene_id AS gene_db_id,
            g.source AS gene_source,
            g.stable_id AS gene_id,
            g.version AS gene_version,
            g.seq_region_start AS gene_start,
            g.seq_region_end AS gene_end,
            g.seq_region_strand AS gene_strand,
            g.biotype AS gene_biotype,
            gx.display_label AS gene_display_name,
            t.transcript_id AS transcript_db_id,
            t.stable_id AS transcript_id,
            t.version AS transcript_version,
            t.seq_region_start AS transcript_start,
            t.seq_region_end AS transcript_end,
            t.seq_region_strand AS transcript_strand,
            t.biotype AS transcript_biotype,
            t.source AS transcript_source,
            tx.display_label AS transcript_display_name,
            e.exon_id AS exon_db_id,
            e.stable_id AS exon_id,
            e.version AS exon_version,
            e.seq_region_id AS exon_seq_region_id,
            e.seq_region_start AS exon_start,
            e.seq_region_end AS exon_end,
            e.phase AS exon_phase,
            e.end_phase AS exon_end_phase,
            SUBSTRING(d.sequence, e.seq_region_start, e.seq_region_end - e.seq_region_start + 1) AS exon_sequence,
            et.rank AS exon_rank,
            tr.translation_id AS translation_db_id,
            tr.stable_id AS protein_id,
            tr.version AS protein_version,
            tr.start_exon_id AS translation_start_exon_id,
            tr.end_exon_id AS translation_end_exon_id,
            tr.seq_start AS translation_seq_start,
            tr.seq_end AS translation_seq_end
        FROM gene g
        JOIN seq_region sr ON sr.seq_region_id = g.seq_region_id
        LEFT JOIN xref gx ON gx.xref_id = g.display_xref_id
        LEFT JOIN analysis ga ON ga.analysis_id = g.analysis_id
        LEFT JOIN transcript t ON {transcript_join}
        LEFT JOIN xref tx ON tx.xref_id = t.display_xref_id
        LEFT JOIN exon_transcript et ON et.transcript_id = t.transcript_id
        LEFT JOIN exon e ON e.exon_id = et.exon_id
        LEFT JOIN dna d ON d.seq_region_id = e.seq_region_id
        LEFT JOIN translation tr ON tr.transcript_id = t.transcript_id
        WHERE {where}
        ORDER BY sr.name, g.seq_region_start, g.stable_id, t.seq_region_start, t.stable_id, et.rank
        """
    )
    expanding = [bindparam('region_ids', expanding=True)]
    if transcript_ids:
        expanding.append(bindparam("transcript_ids", expanding=True))
    if gene_ids:
        expanding.append(bindparam("gene_ids", expanding=True))
    if biotypes:
        expanding.append(bindparam("biotypes", expanding=True))
    if expanding:
        sql = sql.bindparams(*expanding)
    with sql_connection(secrets) as conn:
        result = conn.execute(sql, params)
        return [dict(row._mapping) for row in result]


def fetch_attribute_rows(secrets: dict[str, Any], gene_db_ids: set[int], transcript_db_ids: set[int]) -> tuple[dict[int, list[tuple[str, str]]], dict[int, list[tuple[str, str]]]]:
    from sqlalchemy import bindparam, text

    gene_attrs: dict[int, list[tuple[str, str]]] = defaultdict(list)
    transcript_attrs: dict[int, list[tuple[str, str]]] = defaultdict(list)
    with sql_connection(secrets) as conn:
        if gene_db_ids:
            sql = text(
                """
                SELECT ga.gene_id, at.code, ga.value
                FROM gene_attrib ga
                JOIN attrib_type at ON at.attrib_type_id = ga.attrib_type_id
                WHERE ga.gene_id IN :ids
                """
            ).bindparams(bindparam("ids", expanding=True))
            for row in conn.execute(sql, {"ids": list(gene_db_ids)}):
                gene_attrs[int(row[0])].append((row[1], row[2]))
        if transcript_db_ids:
            sql = text(
                """
                SELECT ta.transcript_id, at.code, ta.value
                FROM transcript_attrib ta
                JOIN attrib_type at ON at.attrib_type_id = ta.attrib_type_id
                WHERE ta.transcript_id IN :ids
                """
            ).bindparams(bindparam("ids", expanding=True))
            for row in conn.execute(sql, {"ids": list(transcript_db_ids)}):
                transcript_attrs[int(row[0])].append((row[1], row[2]))
    return gene_attrs, transcript_attrs


def fetch_exon_sequences(secrets: dict[str, Any], records: list[TranscriptRecord]) -> dict[str, str]:
    from .sequence import SQLSequenceStore, SequenceError, SequenceResolver

    required = required_sequence_exon_ids(records)
    sequences: dict[str, str] = {}
    if not required:
        return sequences
    with sql_connection(secrets) as conn:
        resolver = SequenceResolver(SQLSequenceStore(conn))
        for record in records:
            for exon in record.exons:
                if exon.exon_db_id not in required or exon.exon_db_id in sequences:
                    continue
                try:
                    if exon.seq_region_id is None:
                        raise SequenceError("Missing exon seq_region_id")
                    sequences[exon.exon_db_id] = resolver.fetch(exon.seq_region_id, exon.start, exon.end)
                except SequenceError as error:
                    raise ValueError(
                        f"DNA sequence retrieval failed for transcript {record.transcript_id}, "
                        f"exon {exon.exon_id}, region {exon.seq_region_id}:{exon.start}-{exon.end}: {error}"
                    ) from error
    return sequences


def rows_to_records(
    rows: list[dict[str, Any]],
    no_version: bool = False,
    add_chr: bool = False,
    mane_map: dict[str, str] | None = None,
    gene_attrs: dict[int, list[tuple[str, str]]] | None = None,
    transcript_attrs: dict[int, list[tuple[str, str]]] | None = None,
    exon_sequences: dict[str, str] | None = None,
) -> list[TranscriptRecord]:
    mane_map = mane_map or {}
    gene_attrs = gene_attrs or {}
    transcript_attrs = transcript_attrs or {}
    exon_sequences = exon_sequences or {}
    records: dict[str, TranscriptRecord] = {}
    seen_exons: set[tuple[str, str, int]] = set()
    for row in rows:
        if row.get('transcript_db_id') is None:
            continue
        gene_stable = str(row["gene_id"])
        transcript_stable = str(row["transcript_id"] or row["gene_id"])
        gene_id = stable_with_version(gene_stable, row.get("gene_version"), no_version)
        transcript_id = stable_with_version(transcript_stable, row.get("transcript_version"), no_version)
        protein_stable = row.get("protein_id")
        protein_id = stable_with_version(protein_stable, row.get("protein_version"), no_version) if protein_stable else None
        chrom = str(row["chrom"])
        if add_chr and not chrom.startswith("chr") and not chrom.startswith(("HG", "HS", "KI", "GL")):
            chrom = f"chr{chrom}"
        gidb = int(row["gene_db_id"])
        tidb = int(row["transcript_db_id"])
        gene_name = annotation_name(row, 'gene', gene_stable, gene_attrs.get(gidb, []))
        transcript_name = annotation_name(row, 'transcript', transcript_stable)
        if transcript_id not in records:
            tx_tags = []
            gene_tags = []
            remarks = []
            raw_transcript_attrs = []
            cds_start_nf = False
            cds_end_nf = False
            mrna_start_nf = False
            mrna_end_nf = False
            split_transcript = False
            for code, value in transcript_attrs.get(tidb, []):
                raw_transcript_attrs.append((code, value))
                if code in {'cds_start_NF', 'cds_end_NF', 'mRNA_start_NF', 'mRNA_end_NF'} and not float(value or 0):
                    continue
                if code.lower() in {'mane_select', 'mane_plus_clinical'}:
                    tx_tags.append(SELECTION_REMARK_TAGS[code.lower()])
                elif code == "cds_start_NF":
                    cds_start_nf = True
                    tx_tags.append(code)
                elif code == "cds_end_NF":
                    cds_end_nf = True
                    tx_tags.append(code)
                elif code == "mRNA_start_NF":
                    mrna_start_nf = True
                    tx_tags.append(code)
                elif code == "mRNA_end_NF":
                    mrna_end_nf = True
                    tx_tags.append(code)
                elif code == "split_tscript":
                    split_transcript = True
                elif code in {'remark', 'hidden_remark', 'havana_cv'}:
                    if code == 'remark':
                        remarks.append(value)
                    if value is not None and value.lower() in SELECTION_REMARK_TAGS:
                        # MANE status comes from explicit attributes or remarks,
                        # not from incidental mentions in hidden curation notes.
                        if code == 'remark' or not value.lower().startswith('mane_'):
                            tx_tags.append(SELECTION_REMARK_TAGS[value.lower()])
                    if value is not None and (value in REMARK_TAGS or re.fullmatch(r'[Nn]o.+org.+sup.+ted', value)):
                        tx_tags.append(perl_tag(value))
            for code, value in gene_attrs.get(gidb, []):
                if code == "remark":
                    if value == "not for VEGA":
                        gene_tags.append("not_for_VEGA")
            if transcript_stable in mane_map:
                tx_tags.append("MANE_RefSeq")
            records[transcript_id] = TranscriptRecord(
                chrom=chrom,
                source=row.get("gene_source") or row.get("transcript_source") or "ensembl",
                gene_db_id=gidb,
                gene_id=gene_id,
                gene_stable_id=gene_stable,
                gene_version=row.get("gene_version"),
                gene_name=gene_name,
                gene_biotype=row.get("gene_biotype") or "",
                gene_start=int(row["gene_start"]),
                gene_end=int(row["gene_end"]),
                gene_strand=strand_symbol(row["gene_strand"]),
                transcript_id=transcript_id,
                transcript_db_id=tidb,
                transcript_stable_id=transcript_stable,
                transcript_version=row.get("transcript_version"),
                transcript_name=transcript_name,
                transcript_biotype=row.get("transcript_biotype") or "",
                transcript_start=int(row["transcript_start"]),
                transcript_end=int(row["transcript_end"]),
                transcript_strand=strand_symbol(row["transcript_strand"]),
                transcript_source=row.get("transcript_source") or row.get("gene_source") or "ensembl",
                protein_id=protein_id,
                protein_stable_id=protein_stable,
                protein_version=row.get("protein_version"),
                translation_start_exon_id=str(row["translation_start_exon_id"]) if row.get("translation_start_exon_id") else None,
                translation_end_exon_id=str(row["translation_end_exon_id"]) if row.get("translation_end_exon_id") else None,
                seq_start=int(row["translation_seq_start"]) if row.get("translation_seq_start") is not None else None,
                seq_end=int(row["translation_seq_end"]) if row.get("translation_seq_end") is not None else None,
                peptide_sequence=str(row.get("peptide_sequence") or row.get("translation_sequence") or "") or None,
                tags=unique_preserve_order(tx_tags),
                gene_tags=unique_preserve_order(gene_tags),
                transcript_remarks=remarks,
                transcript_attributes=raw_transcript_attrs,
                cds_start_nf=cds_start_nf,
                cds_end_nf=cds_end_nf,
                mrna_start_nf=mrna_start_nf,
                mrna_end_nf=mrna_end_nf,
                split_transcript=split_transcript,
            )
        if row.get('exon_db_id') is None:
            continue
        exon_id = stable_with_version(row["exon_id"], row.get("exon_version"), no_version)
        exon_key = (transcript_id, exon_id, int(row["exon_rank"]))
        if exon_key not in seen_exons:
            seen_exons.add(exon_key)
            transcript_strand = records[transcript_id].transcript_strand
            records[transcript_id].exons.append(
                ExonRecord(
                    exon_db_id=str(row["exon_db_id"]),
                    exon_id=exon_id,
                    start=int(row["exon_start"]),
                    end=int(row["exon_end"]),
                    rank=int(row["exon_rank"]),
                    seq_region_id=str(row["exon_seq_region_id"]) if row.get("exon_seq_region_id") is not None else None,
                    version=row.get("exon_version"),
                    phase=str(row.get("exon_phase") if row.get("exon_phase") is not None else "."),
                    end_phase=str(row.get("exon_end_phase") if row.get("exon_end_phase") is not None else "."),
                    sequence=orient_exon_sequence(row.get("exon_sequence") or exon_sequences.get(str(row["exon_db_id"])), transcript_strand),
                )
            )
    for record in records.values():
        record.exons.sort(key=lambda exon: exon.rank)
    return list(records.values())


def required_sequence_exon_ids(records: list[TranscriptRecord]) -> set[str]:
    ids: set[str] = set()
    for record in records:
        if record.seq_start is None or record.seq_end is None or not record.translation_start_exon_id or not record.translation_end_exon_id:
            continue
        start_idx = exon_index_by_db_id(record, record.translation_start_exon_id)
        end_idx = exon_index_by_db_id(record, record.translation_end_exon_id)
        if start_idx is None or end_idx is None:
            continue
        lo, hi = min(start_idx, end_idx), max(start_idx, end_idx)
        ids.update(record.exons[idx].exon_db_id for idx in range(lo, hi + 1))
    return ids


def attach_exon_sequences(records: list[TranscriptRecord], exon_sequences: dict[str, str]) -> None:
    for record in records:
        for exon in record.exons:
            seq = exon_sequences.get(exon.exon_db_id)
            if seq:
                exon.sequence = orient_exon_sequence(seq, record.transcript_strand)


def attrs_for_feature(record: TranscriptRecord, feature_type: str, fmt: str, extra: dict[str, Any] | None = None) -> dict[str, Any]:
    extra = extra or {}
    if fmt == "gtf":
        attrs: dict[str, Any] = {
            "gene_id": record.gene_id,
            "transcript_id": record.transcript_id if feature_type != "gene" else None,
            "gene_type": record.gene_biotype,
            "gene_name": record.gene_name,
            "transcript_type": record.transcript_biotype if feature_type != "gene" else None,
            "transcript_name": record.transcript_name if feature_type != "gene" else None,
        }
        if record.protein_id and feature_type != "gene":
            attrs["protein_id"] = record.protein_id
        if record.tags and feature_type != "gene":
            attrs["tag"] = ",".join(record.tags)
        attrs.update(extra)
        return attrs
    attrs = {
        "gene_id": record.gene_id,
        "gene_type": record.gene_biotype,
        "gene_name": record.gene_name,
        "transcript_id": record.transcript_id if feature_type != "gene" else None,
        "transcript_type": record.transcript_biotype if feature_type != "gene" else None,
        "transcript_name": record.transcript_name if feature_type != "gene" else None,
    }
    if record.protein_id and feature_type != "gene":
        attrs["protein_id"] = record.protein_id
    if record.tags and feature_type != "gene":
        attrs["tag"] = ",".join(record.tags)
    attrs.update(extra)
    return attrs


def coding_segments(record: TranscriptRecord) -> list[tuple[int, int, ExonRecord]]:
    if record.seq_start is None or record.seq_end is None or not record.translation_start_exon_id or not record.translation_end_exon_id:
        return []
    start_exon_idx = None
    end_exon_idx = None
    for idx, exon in enumerate(record.exons):
        if exon.exon_db_id == record.translation_start_exon_id:
            start_exon_idx = idx
        if exon.exon_db_id == record.translation_end_exon_id:
            end_exon_idx = idx
    if start_exon_idx is None or end_exon_idx is None:
        return []
    lo = min(start_exon_idx, end_exon_idx)
    hi = max(start_exon_idx, end_exon_idx)
    segments = []
    for idx in range(lo, hi + 1):
        exon = record.exons[idx]
        seg_start, seg_end = exon.start, exon.end
        if idx == start_exon_idx:
            if record.transcript_strand == "-":
                seg_end = exon.end - int(record.seq_start) + 1
            else:
                seg_start = exon.start + int(record.seq_start) - 1
        if idx == end_exon_idx:
            if record.transcript_strand == "-":
                seg_start = exon.end - int(record.seq_end) + 1
            else:
                seg_end = exon.start + int(record.seq_end) - 1
        if seg_start <= seg_end:
            segments.append((seg_start, seg_end, exon))
    return sorted(segments)


def phased_coding_segments(record: TranscriptRecord) -> list[tuple[int, int, ExonRecord, str]]:
    segments = coding_segments(record)
    if not segments:
        return []
    transcript_order = segments if record.transcript_strand == "+" else list(reversed(segments))
    phase_by_segment: dict[tuple[int, int, str], str] = {}
    phase = 0
    for start, end, exon in transcript_order:
        if exon.phase in {"0", "1", "2"}:
            phase_by_segment[(start, end, exon.exon_db_id)] = str((3 - int(exon.phase)) % 3)
        else:
            phase_by_segment[(start, end, exon.exon_db_id)] = str(phase)
        phase = (3 - ((end - start + 1 - phase) % 3)) % 3
    return [(start, end, exon, phase_by_segment[(start, end, exon.exon_db_id)]) for start, end, exon in segments]


def subtract_intervals(start: int, end: int, intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    pieces = [(start, end)]
    for cut_start, cut_end in intervals:
        next_pieces = []
        for piece_start, piece_end in pieces:
            if cut_end < piece_start or cut_start > piece_end:
                next_pieces.append((piece_start, piece_end))
                continue
            if piece_start < cut_start:
                next_pieces.append((piece_start, cut_start - 1))
            if cut_end < piece_end:
                next_pieces.append((cut_end + 1, piece_end))
        pieces = next_pieces
    return pieces


def trim_stop_from_cds_segments(segments: list[tuple[int, int, ExonRecord]], record: TranscriptRecord) -> list[tuple[int, int, ExonRecord]]:
    if not segments or not should_emit_codon(record, "stop"):
        return list(segments)
    stop_intervals = [(start, end) for start, end, _, _ in codon_segments(record, "stop")]
    trimmed = []
    for start, end, exon in segments:
        for piece_start, piece_end in subtract_intervals(start, end, stop_intervals):
            trimmed.append((piece_start, piece_end, exon))
    return sorted(trimmed)


def trim_stop_from_phased_segments(segments: list[tuple[int, int, ExonRecord, str]], record: TranscriptRecord) -> list[tuple[int, int, ExonRecord, str]]:
    if not segments or not should_emit_codon(record, "stop"):
        return list(segments)
    stop_intervals = [(start, end) for start, end, _, _ in codon_segments(record, "stop")]
    trimmed = []
    for start, end, exon, phase in segments:
        for piece_start, piece_end in subtract_intervals(start, end, stop_intervals):
            trimmed.append((piece_start, piece_end, exon, phase))
    return sorted(trimmed)


def utr_segments_from_cds(record: TranscriptRecord, cds: list[tuple[int, int, ExonRecord]]) -> tuple[list[tuple[int, int, ExonRecord]], list[tuple[int, int, ExonRecord]]]:
    if not cds:
        return [], []
    cds_start = min(start for start, _, _ in cds)
    cds_end = max(end for _, end, _ in cds)
    five_prime = []
    three_prime = []
    for exon in record.exons:
        pieces = []
        if exon.start < cds_start:
            pieces.append((exon.start, min(exon.end, cds_start - 1), exon))
        if exon.end > cds_end:
            pieces.append((max(exon.start, cds_end + 1), exon.end, exon))
        for piece in pieces:
            if piece[0] > piece[1]:
                continue
            if record.transcript_strand == "+":
                (five_prime if piece[1] < cds_start else three_prime).append(piece)
            else:
                (three_prime if piece[1] < cds_start else five_prime).append(piece)
    return sorted(five_prime), sorted(three_prime)


def utr_segments(record: TranscriptRecord) -> tuple[list[tuple[int, int, ExonRecord]], list[tuple[int, int, ExonRecord]]]:
    return utr_segments_from_cds(record, coding_segments(record))


def codon_segments(record: TranscriptRecord, codon: str) -> list[tuple[int, int, ExonRecord, str]]:
    cds = coding_segments(record)
    if not cds:
        return []
    bases_needed = 3
    segments = cds if (record.transcript_strand == "+" and codon == "start") or (record.transcript_strand == "-" and codon == "stop") else list(reversed(cds))
    out: list[tuple[int, int, ExonRecord]] = []
    for start, end, exon in segments:
        length = end - start + 1
        take = min(length, bases_needed)
        if segments is cds:
            out.append((start, start + take - 1, exon))
        else:
            out.append((end - take + 1, end, exon))
        bases_needed -= take
        if bases_needed == 0:
            break
    genomic_order = sorted(out)
    transcript_order = genomic_order if record.transcript_strand == "+" else list(reversed(genomic_order))
    phase_by_segment: dict[tuple[int, int, str], str] = {}
    bases_consumed = 0
    for start, end, exon in transcript_order:
        phase_by_segment[(start, end, exon.exon_db_id)] = "0" if bases_consumed == 0 else str((3 - (bases_consumed % 3)) % 3)
        bases_consumed += end - start + 1
    return [(start, end, exon, phase_by_segment[(start, end, exon.exon_db_id)]) for start, end, exon in genomic_order]


def codon_boundary_exon(record: TranscriptRecord, codon: str) -> ExonRecord | None:
    boundary_id = record.translation_start_exon_id if codon == "start" else record.translation_end_exon_id
    if not boundary_id:
        return None
    for exon in record.exons:
        if exon.exon_db_id == boundary_id:
            return exon
    return None


def has_transcript_attribute(record: TranscriptRecord, code: str, value: str | None = None) -> bool:
    for attr_code, attr_value in record.transcript_attributes:
        if attr_code == code and (value is None or attr_value == value):
            return True
    return False


GENETIC_CODE = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}


def exon_index_by_db_id(record: TranscriptRecord, exon_db_id: str | None) -> int | None:
    if not exon_db_id:
        return None
    for idx, exon in enumerate(record.exons):
        if exon.exon_db_id == exon_db_id:
            return idx
    return None


def sequence_validation_available(record: TranscriptRecord) -> bool:
    if record.seq_start is None or record.seq_end is None or not record.translation_start_exon_id or not record.translation_end_exon_id:
        return False
    start_idx = exon_index_by_db_id(record, record.translation_start_exon_id)
    end_idx = exon_index_by_db_id(record, record.translation_end_exon_id)
    if start_idx is None or end_idx is None:
        return False
    lo, hi = min(start_idx, end_idx), max(start_idx, end_idx)
    return all(
        record.exons[idx].sequence is not None
        and len(record.exons[idx].sequence) == record.exons[idx].end - record.exons[idx].start + 1
        for idx in range(lo, hi + 1)
    )


def translateable_sequence(record: TranscriptRecord) -> str | None:
    if not sequence_validation_available(record):
        return None
    start_idx = exon_index_by_db_id(record, record.translation_start_exon_id)
    end_idx = exon_index_by_db_id(record, record.translation_end_exon_id)
    if start_idx is None or end_idx is None:
        return None
    lo, hi = min(start_idx, end_idx), max(start_idx, end_idx)
    chunks = []
    for idx in range(lo, hi + 1):
        exon_seq = record.exons[idx].sequence
        if exon_seq is None:
            return None
        start = 0
        end = len(exon_seq)
        if idx == start_idx and record.seq_start is not None:
            start = max(int(record.seq_start) - 1, 0)
        if idx == end_idx and record.seq_end is not None:
            end = min(int(record.seq_end), len(exon_seq))
        if start < end:
            chunks.append(exon_seq[start:end])
    return "".join(chunks).upper()


def translated_peptide(seq: str) -> str:
    amino_acids = []
    for idx in range(0, len(seq) - 2, 3):
        amino_acids.append(GENETIC_CODE.get(seq[idx : idx + 3], "X"))
    if amino_acids and amino_acids[-1] == "*":
        amino_acids.pop()
    return "".join(amino_acids)


def initial_coding_phase(record: TranscriptRecord) -> int:
    segments = phased_coding_segments(record)
    if not segments:
        return 0
    transcript_order = segments if record.transcript_strand == "+" else list(reversed(segments))
    phase = transcript_order[0][3]
    return int(phase) if phase in {"0", "1", "2"} else 0


def inferred_peptide_length(record: TranscriptRecord, seq: str, endseq: str) -> int:
    phase = initial_coding_phase(record)
    coding_bases = max(len(seq) - phase, 0)
    if endseq in {"TAG", "TGA", "TAA"} and coding_bases >= 3 and coding_bases % 3 == 0:
        coding_bases -= 3
    return coding_bases // 3


def has_non_atg_start_attribute(record: TranscriptRecord) -> bool:
    return any(
        code in {"remark", "hidden_remark", "havana_cv"} and value == "non-ATG start"
        for code, value in record.transcript_attributes
    )


def perl_start_and_stop_presence(record: TranscriptRecord) -> tuple[bool, bool] | None:
    seq = translateable_sequence(record)
    if seq is None or len(seq) < 3:
        return None
    startseq = seq[:3]
    endseq = seq[-3:]
    peptide = record.peptide_sequence if record.peptide_sequence is not None else translated_peptide(seq)
    peptide_length = len(record.peptide_sequence) if record.peptide_sequence is not None else inferred_peptide_length(record, seq, endseq)

    has_start = True
    if startseq != "ATG":
        has_start = has_non_atg_start_attribute(record)
    elif not peptide.startswith("M"):
        has_start = False

    has_end = True
    if endseq not in {"TAG", "TGA", "TAA"}:
        has_end = False
    elif len(seq) - 3 * peptide_length < 3:
        has_end = False
    return has_start, has_end


def should_emit_codon(record: TranscriptRecord, codon_name: str) -> bool:
    presence = perl_start_and_stop_presence(record)
    if codon_name == "start":
        if record.cds_start_nf:
            return False
        return True if presence is None else presence[0]
    if record.cds_end_nf:
        return False
    return True if presence is None else presence[1]


def require_sequence_validation(records: list[TranscriptRecord]) -> None:
    missing = [
        record.transcript_id
        for record in records
        if record.translation_start_exon_id and record.translation_end_exon_id and not sequence_validation_available(record)
    ]
    if missing:
        shown = ", ".join(missing[:10])
        suffix = f" and {len(missing) - 10} more" if len(missing) > 10 else ""
        fail(
            "Live DB export requires DNA sequence for Perl-style codon validation; "
            f"missing sequence context for {shown}{suffix}."
        )


def hidden_remark_values(attributes):
    return sorted(value for code, value in attributes
                  if code == 'hidden_remark' and value not in (None, ''))


def gene_features_from_rows(rows, gene_attrs, fmt, no_version=False, add_chr=False,
                            include_hidden_remarks=False):
    genes = {}
    for row in rows:
        db_id = int(row['gene_db_id'])
        if db_id in genes:
            continue
        stable_id = str(row['gene_id'])
        gene_id = stable_with_version(stable_id, row.get('gene_version'), no_version)
        attrs = {'ID': gene_id} if fmt == 'gff3' else {}
        attrs.update(gene_id=gene_id, gene_type=row['gene_biotype'],
                     gene_name=annotation_name(row, 'gene', stable_id, gene_attrs.get(db_id, [])))
        if ('remark', 'not for VEGA') in gene_attrs.get(db_id, []):
            attrs['tag'] = 'not_for_VEGA'
        if include_hidden_remarks:
            remarks = hidden_remark_values(gene_attrs.get(db_id, []))
            if remarks:
                attrs['hidden_remark'] = remarks
        chrom = str(row['chrom'])
        if add_chr and not chrom.startswith(('chr', 'HG', 'HS', 'KI', 'GL')):
            chrom = 'chr' + chrom
        genes[db_id] = Feature(chrom, row['gene_source'] or 'ensembl', 'gene',
                               int(row['gene_start']), int(row['gene_end']), '.',
                               strand_symbol(row['gene_strand']), '.', attrs)
    return genes


def records_to_features(records: list[TranscriptRecord], fmt: str, mane_map: dict[str, str] | None = None,
                        genes: dict[int, Feature] | None = None,
                        include_hidden_remarks: bool = False) -> list[Feature]:
    mane_map = mane_map or {}
    features: list[Feature] = []
    seen_genes: set[str] = set()
    for record in sorted(records, key=lambda r: (r.chrom, r.gene_db_id, r.transcript_db_id)):
        if record.gene_id not in seen_genes:
            seen_genes.add(record.gene_id)
            extra = {"ID": record.gene_id} if fmt == "gff3" else {}
            if fmt == "gff3" and record.gene_tags:
                extra["tag"] = ",".join(record.gene_tags)
            features.append(
                Feature(
                    record.chrom,
                    record.source,
                    "gene",
                    record.gene_start,
                    record.gene_end,
                    ".",
                    record.gene_strand,
                    ".",
                    attrs_for_feature(record, "gene", fmt, extra),
                )
            )
        tx_extra: dict[str, Any] = {}
        if fmt == "gff3":
            tx_extra.update({"ID": record.transcript_id, "Parent": record.gene_id})
        if record.tags:
            tx_extra["tag"] = ",".join(record.tags)
        if record.protein_id:
            tx_extra["protein_id"] = record.protein_id
        if record.transcript_stable_id in mane_map:
            tx_extra["Dbxref" if fmt == "gff3" else "db_xref"] = f"RefSeq:{mane_map[record.transcript_stable_id]}"
        if include_hidden_remarks:
            remarks = hidden_remark_values(record.transcript_attributes)
            if remarks:
                tx_extra['hidden_remark'] = remarks
        features.append(
            Feature(
                record.chrom,
                record.transcript_source,
                "transcript",
                record.transcript_start,
                record.transcript_end,
                ".",
                record.transcript_strand,
                ".",
                attrs_for_feature(record, "transcript", fmt, tx_extra),
            )
        )

        coding_children: dict[str, list[Feature]] = defaultdict(list)
        for exon in record.exons:
            exon_extra = {
                "exon_id": exon.exon_id,
                "exon_number": exon.rank,
            }
            if fmt == "gff3":
                exon_extra.update({"ID": f"exon:{record.transcript_id}:{exon.rank}", "Parent": record.transcript_id})
            coding_children[exon.exon_db_id].append(
                Feature(
                    record.chrom,
                    record.transcript_source,
                    "exon",
                    exon.start,
                    exon.end,
                    ".",
                    record.transcript_strand,
                    ".",
                    attrs_for_feature(record, "exon", fmt, exon_extra),
                )
            )
        cds_segments = coding_segments(record)
        feature_cds_segments = trim_stop_from_cds_segments(cds_segments, record) if fmt == "gtf" else cds_segments
        utr_source_segments = feature_cds_segments if fmt == "gtf" else cds_segments
        five_utr, three_utr = utr_segments_from_cds(record, utr_source_segments)
        if not (record.transcript_biotype in CODING_BIOTYPES or
                re.fullmatch(r'(IG|TR)_.+_gene', record.transcript_biotype)):
            five_utr, three_utr = [], []
        utr_features: list[Feature] = []
        for label, segments in [
            ("five_prime_UTR", five_utr),
            ("three_prime_UTR", three_utr),
        ]:
            for idx, (start, end, exon) in enumerate(sorted(segments, key=lambda item: item[2].rank), 1):
                feature_label = "UTR" if fmt == "gtf" else label
                extra = {}
                if fmt == "gff3":
                    utr_prefix = "UTR5" if label == "five_prime_UTR" else "UTR3"
                    extra = {"ID": f"{utr_prefix}:{record.transcript_id}", "Parent": record.transcript_id}
                extra["exon_number"] = exon.rank
                extra["exon_id"] = exon.exon_id
                utr_features.append(
                    Feature(record.chrom, record.transcript_source, feature_label, start, end, ".", record.transcript_strand, ".", attrs_for_feature(record, feature_label, fmt, extra))
                )
        phased_segments = phased_coding_segments(record)
        if fmt == "gtf":
            phased_segments = trim_stop_from_phased_segments(phased_segments, record)
        for start, end, exon, phase in phased_segments:
            extra = {}
            if fmt == "gff3":
                extra = {"ID": f"CDS:{record.transcript_id}", "Parent": record.transcript_id}
            extra["exon_number"] = exon.rank
            extra["exon_id"] = exon.exon_id
            if record.protein_id:
                extra["protein_id"] = record.protein_id
            coding_children[exon.exon_db_id].append(
                Feature(record.chrom, record.transcript_source, "CDS", start, end, ".", record.transcript_strand, phase, attrs_for_feature(record, "CDS", fmt, extra))
            )
        for codon_name in ["start", "stop"]:
            if not should_emit_codon(record, codon_name):
                continue
            boundary_exon = codon_boundary_exon(record, codon_name)
            for start, end, exon, phase in codon_segments(record, codon_name):
                feature_type = f"{codon_name}_codon"
                extra = {}
                if fmt == "gff3":
                    extra = {"ID": f"{feature_type}:{record.transcript_id}", "Parent": record.transcript_id}
                attr_exon = boundary_exon or exon
                extra["exon_number"] = exon.rank
                extra["exon_id"] = attr_exon.exon_id
                order_exon = boundary_exon or exon
                coding_children[order_exon.exon_db_id].append(
                    Feature(record.chrom, record.transcript_source, feature_type, start, end, ".", record.transcript_strand, phase, attrs_for_feature(record, feature_type, fmt, extra))
                )
        feature_order = {"exon": 0, "CDS": 1, "start_codon": 2, "stop_codon": 3}

        def child_order(feature: Feature) -> tuple[int, int, int]:
            start_order = -feature.start if feature.feature_type in {"start_codon", "stop_codon"} and record.transcript_strand == "-" else feature.start
            return feature_order.get(feature.feature_type, 99), start_order, feature.end

        for exon in record.exons:
            features.extend(sorted(coding_children[exon.exon_db_id], key=child_order))
        features.extend(sorted(utr_features, key=lambda feature: (int(feature.attrs.get("exon_number", 0)), feature.start, feature.end, feature.feature_type)))
    if genes is not None:
        # Preserve gene/transcript traversal while retaining genes without surviving transcripts.
        children = defaultdict(list)
        for feature in features:
            if feature.feature_type != 'gene':
                children[feature.attrs['gene_id']].append(feature)
        features = []
        for _, gene in sorted(genes.items(), key=lambda item: (item[1].chrom, item[0])):
            features.append(gene)
            features.extend(children[gene.attrs['gene_id']])
    return features


def feature_to_line(feature: Feature, fmt: str) -> str:
    attr_text = format_gtf_attrs(feature.attrs) if fmt == "gtf" else format_gff3_attrs_perl(feature.attrs, feature.feature_type)
    return "\t".join(
        [
            feature.chrom,
            feature.source,
            feature.feature_type,
            str(feature.start),
            str(feature.end),
            feature.score,
            feature.strand,
            feature.phase,
            attr_text,
        ]
    )


def filter_fixture_features(
    path: str | Path,
    chromosome: str | tuple[str, ...] | None,
    start: int | None,
    end: int | None,
    gene_ids: set[str] | None,
    transcript_ids: set[str] | None,
    biotypes: set[str] | None,
    source: str | None,
    add_chr: bool,
    no_version: bool,
) -> list[Feature]:
    features = []
    allowed_transcripts: set[str] | None = None
    allowed_genes: set[str] | None = None
    with open_text(path) as handle:
        raw_features = [parse_feature(line, None) for line in handle if line.strip() and not line.startswith("#")]
    requested = parse_chromosomes(chromosome)
    available = {feature.chrom for feature in raw_features}
    validate_requested_regions(requested, available, 'fixture input')
    validate_slice_selection(len(requested if requested is not None else available), start, end)
    if gene_ids or transcript_ids or biotypes or source:
        allowed_transcripts = set()
        allowed_genes = set()
        for feature in raw_features:
            attrs = feature.attrs
            if feature.feature_type == "gene":
                gid = attrs.get("gene_id") or attrs.get("ID")
                if gene_ids and gid not in gene_ids and strip_version(gid) not in gene_ids:
                    continue
                if biotypes and (attrs.get("gene_type") or attrs.get("gene_biotype")) not in biotypes:
                    continue
                if source and feature.source != source:
                    continue
                if gid:
                    allowed_genes.add(gid)
            elif feature.feature_type in {"transcript", "mRNA"}:
                tid = attrs.get("transcript_id") or attrs.get("ID")
                gid = attrs.get("gene_id") or attrs.get("Parent")
                if transcript_ids and tid not in transcript_ids and strip_version(tid) not in transcript_ids:
                    continue
                if gene_ids and gid not in gene_ids and strip_version(gid) not in gene_ids:
                    continue
                if biotypes and (attrs.get("transcript_type") or attrs.get("transcript_biotype")) not in biotypes:
                    continue
                if source and feature.source != source:
                    continue
                if tid:
                    allowed_transcripts.add(tid)
                if gid:
                    allowed_genes.add(gid)
    for feature in raw_features:
        if requested is not None and feature.chrom not in requested:
            continue
        if start is not None and feature.end < start:
            continue
        if end is not None and feature.start > end:
            continue
        attrs = dict(feature.attrs)
        if allowed_transcripts is not None:
            tid = attrs.get("transcript_id") or attrs.get("Parent") or attrs.get("ID")
            gid = attrs.get("gene_id") or attrs.get("Parent") or attrs.get("ID")
            if feature.feature_type == "gene":
                if gid not in allowed_genes:
                    continue
            elif tid not in allowed_transcripts:
                continue
        chrom = feature.chrom
        if add_chr and not chrom.startswith("chr") and not chrom.startswith(("HG", "HS", "KI", "GL")):
            chrom = f"chr{chrom}"
        if no_version:
            for key in ["ID", "Parent", "gene_id", "transcript_id", "protein_id", "exon_id"]:
                if key in attrs:
                    attrs[key] = strip_version(attrs[key])
        features.append(Feature(chrom, feature.source, feature.feature_type, feature.start, feature.end, feature.score, feature.strand, feature.phase, attrs))
    return features


def copy_fixture_features(path: str | Path, chromosome: str | None) -> list[Feature]:
    return filter_fixture_features(path, chromosome, None, None, None, None, None, None, False, False)


def strip_version(value: Any) -> str:
    text = str(value or "")
    if "." in text and text.rsplit(".", 1)[1].isdigit():
        return text.rsplit(".", 1)[0]
    return text


def write_features(features: Iterable[Feature], output: str | Path, fmt: str,
                   sequence_regions: Iterable[tuple[str, int, int]] = ()) -> None:
    with Path(output).open("w", encoding="utf-8") as handle:
        if fmt == "gff3":
            handle.write("##gff-version 3\n")
            for name, start, end in sequence_regions:
                handle.write(f'##sequence-region {name} {start} {end}\n')
        for feature in features:
            handle.write(feature_to_line(feature, fmt) + "\n")


def export_main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Export Ensembl-like DB annotations to GFF3/GTF.", allow_abbrev=False)
    parser.add_argument(
        "--db-secrets",
        metavar="JSON",
        help="Database connection JSON; required unless --fixture-gff is used.",
    )
    parser.add_argument(
        "--fixture-gff",
        metavar="FILE",
        help="Read an existing GFF3/GTF file instead of querying a live database.",
    )
    parser.add_argument(
        "--chromosome",
        "--chromosomes",
        "--chr",
        dest="chromosome",
        metavar="NAMES",
        help="Comma-separated exact default-assembly region names, e.g. 1,2,X; may include named patches/scaffolds.",
    )
    parser.add_argument(
        "--region-set", choices=['CHR', 'PRI', 'ALL'],
        help="CHR: reference chromosomes (live default); PRI: chromosomes and reference scaffolds; ALL: also patches/haplotypes. Cannot combine with --chromosome.",
    )
    parser.add_argument("--start", type=int, metavar="POS", help="Keep annotation overlapping this 1-based lower bound; requires one selected region.")
    parser.add_argument("--end", type=int, metavar="POS", help="Keep annotation overlapping this 1-based upper bound; requires one selected region.")
    parser.add_argument(
        "--transcript-list",
        "--tlist",
        dest="transcript_list",
        metavar="FILE",
        help="Restrict output to transcript stable IDs listed one per line.",
    )
    parser.add_argument(
        "--gene-list",
        "--list",
        dest="gene_list",
        metavar="FILE",
        help="Restrict output to gene stable IDs listed one per line.",
    )
    parser.add_argument(
        "--biotype",
        action="append",
        metavar="NAME",
        help="Keep a gene or transcript biotype; repeat to allow multiple values.",
    )
    parser.add_argument("--source", metavar="NAME", help="Keep records matching this exact gene or transcript source.")
    parser.add_argument("--analysis", metavar="LOGIC_NAME", help="Keep genes with this exact analysis logic name (database mode only).")
    parser.add_argument(
        "--format",
        choices=["gff3", "gtf"],
        default="gff3",
        help="Output format (default: gff3).",
    )
    parser.add_argument(
        "--no-version",
        "--no_version",
        dest="no_version",
        action="store_true",
        help="Remove numeric version suffixes from stable identifiers.",
    )
    parser.add_argument(
        "--add-chr",
        "--add_chr",
        dest="add_chr",
        action="store_true",
        help="Add a chr prefix to standard sequence-region names.",
    )
    parser.add_argument("--mane", metavar="TSV", help="Add MANE RefSeq cross-references from a MANE summary TSV.")
    parser.add_argument(
        "--artifacts",
        action="store_true",
        help="Include all artifact-biotype models; otherwise only parent-gene artef_dupl exceptions survive.",
    )
    parser.add_argument(
        "--vega",
        action="store_true",
        help="Exclude genes/transcripts carrying the remark 'not for VEGA'.",
    )
    parser.add_argument(
        "--include-hidden-remarks",
        action="store_true",
        help="Include full hidden remarks on their owning gene/transcript rows (database mode only; default: off). May expose internal curation notes.",
    )
    parser.add_argument(
        "--output",
        "--file",
        dest="output",
        metavar="FILE",
        required=True,
        help="Output GFF3/GTF path (required).",
    )
    argv = list(sys.argv[1:] if argv is None else argv)
    if any(arg.split('=', 1)[0] in {'--non-ref', '--non_ref'} for arg in argv):
        parser.error('--non-ref/--non_ref were removed; use --region-set ALL for all regions, or omit the old flag with an explicit --chromosome list')
    args = parser.parse_args(argv)
    try:
        args.chromosome = validate_region_selection(args.chromosome, args.region_set)
    except ValueError as error:
        parser.error(str(error))

    gene_ids = read_id_list(args.gene_list)
    transcript_ids = read_id_list(args.transcript_list)
    biotypes = set(args.biotype or []) or None
    mane_map = load_mane_map(args.mane)

    if args.fixture_gff:
        if args.region_set in {'CHR', 'PRI'}:
            parser.error('Fixture input lacks assembly classification metadata; use explicit --chromosome names or --region-set ALL')
        try:
            features = filter_fixture_features(
                args.fixture_gff,
                args.chromosome,
                args.start,
                args.end,
                gene_ids,
                transcript_ids,
                biotypes,
                args.source,
                args.add_chr,
                args.no_version,
            )
        except ValueError as error:
            parser.error(str(error))
        write_features(features, args.output, args.format)
        return 0
    if not args.db_secrets:
        fail("--db-secrets is required without --fixture-gff")
    secrets = load_secrets(args.db_secrets)
    from .regions import fetch_regions
    try:
        with sql_connection(secrets) as conn:
            regions = fetch_regions(conn, args.chromosome, args.region_set)
        validate_slice_selection(len(regions), args.start, args.end)
    except ValueError as error:
        parser.error(str(error))
    sequence_regions = []
    for region in regions:
        start = args.start if args.start is not None else 1
        end = args.end if args.end is not None else region['length']
        if start < 1 or end < start or end > region['length']:
            fail(f"Invalid slice bounds for {region['name']}: {start}-{end}")
        name = region['name']
        if args.add_chr and not name.startswith(('chr', 'HG', 'HS', 'KI', 'GL')):
            name = 'chr' + name
        sequence_regions.append((name, start, end))
    rows = fetch_annotation_rows(
        secrets,
        chromosome=args.chromosome,
        start=args.start,
        end=args.end,
        transcript_ids=transcript_ids,
        gene_ids=gene_ids,
        biotypes=biotypes,
        source=args.source,
        analysis=args.analysis,
        region_set=args.region_set,
        allow_artifacts=args.artifacts,
        vega_only=args.vega,
        regions=regions,
    )
    gene_db_ids = {int(row["gene_db_id"]) for row in rows if row.get("gene_db_id") is not None}
    transcript_db_ids = {int(row["transcript_db_id"]) for row in rows if row.get("transcript_db_id") is not None}
    gene_attrs, transcript_attrs = fetch_attribute_rows(secrets, gene_db_ids, transcript_db_ids)
    genes = gene_features_from_rows(rows, gene_attrs, args.format, args.no_version, args.add_chr,
                                   include_hidden_remarks=args.include_hidden_remarks)
    records = rows_to_records(rows, args.no_version, args.add_chr, mane_map, gene_attrs, transcript_attrs)
    exon_sequences = fetch_exon_sequences(secrets, records)
    attach_exon_sequences(records, exon_sequences)
    require_sequence_validation(records)
    features = records_to_features(records, args.format, mane_map, genes,
                                   include_hidden_remarks=args.include_hidden_remarks)
    from .par import add_y_par_features
    features = add_y_par_features(secrets, args, regions, features, export_main)
    write_features(features, args.output, args.format, sequence_regions)
    return 0


def chromosomes_main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Get chromosome names from a DB or fixture annotation file.")
    parser.add_argument(
        "--db-secrets",
        metavar="JSON",
        help="Database connection JSON; required unless --fixture-gff is used.",
    )
    parser.add_argument(
        "--fixture-gff",
        metavar="FILE",
        help="Discover names from an existing GFF3/GTF file instead of a database.",
    )
    parser.add_argument(
        "--include-chromosomes",
        metavar="BOOL",
        default="true",
        help="Include chromosome coordinate systems (default: true).",
    )
    parser.add_argument(
        "--include-scaffolds",
        metavar="BOOL",
        default="false",
        help="Include scaffold coordinate systems (default: false).",
    )
    parser.add_argument(
        "--include-patches",
        metavar="BOOL",
        default="false",
        help="Include patch/haplotype sequence regions (default: false).",
    )
    parser.add_argument(
        "--output",
        metavar="FILE",
        default="chromosomes.txt",
        help="Destination text file, one sequence-region name per line (default: chromosomes.txt).",
    )
    args = parser.parse_args(argv)

    if args.fixture_gff:
        chroms = sorted({feature.chrom for feature in filter_fixture_features(args.fixture_gff, None, None, None, None, None, None, None, False, False)})
    else:
        if not args.db_secrets:
            fail("--db-secrets is required without --fixture-gff")
        truthy = {"1", "true", "yes", "y"}
        chroms = fetch_chromosomes_from_db(
            load_secrets(args.db_secrets),
            str(args.include_chromosomes).lower() in truthy,
            str(args.include_scaffolds).lower() in truthy,
            str(args.include_patches).lower() in truthy,
        )
    with Path(args.output).open("w", encoding="utf-8") as handle:
        for chrom in chroms:
            handle.write(f"{chrom}\n")
    return 0


def secrets_main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Validate DB secrets JSON and write normalized TSV.")
    parser.add_argument("db_secrets")
    parser.add_argument("--output", default="db_secrets.tsv")
    args = parser.parse_args(argv)
    secrets = load_secrets(args.db_secrets)
    with Path(args.output).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["host", "port", "user", "password", "database", "dialect"], delimiter="\t")
        writer.writeheader()
        writer.writerow(secrets)
    return 0
