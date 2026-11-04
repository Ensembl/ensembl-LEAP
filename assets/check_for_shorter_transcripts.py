#!/usr/bin/env python3
"""Check that extended GFF3 transcripts retain their original boundaries."""

import argparse
import csv
import gzip
import sys
from collections import Counter
from pathlib import Path


FIELDS = ["transcript_id", "status"] + [
    "{}_{}".format(field, annotation)
    for annotation in ("current", "extended")
    for field in ("chromosome", "strand", "start", "end", "length")
] + ["shorter_left", "shorter_right", "chromosome_mismatch", "strand_mismatch"]
STATUSES = ("unchanged", "extended", "shortened", "location_mismatch",
            "missing_from_extended", "missing_from_raw")


def get_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("current_annotation", help="Raw annotation GFF3 file (optionally .gz)")
    parser.add_argument("extended_annotation", help="Extended annotation GFF3 file (optionally .gz)")
    parser.add_argument("output", help="Comparison CSV output path")
    return parser.parse_args()


def read_transcripts(path):
    transcripts = {}
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as handle:
        for number, line in enumerate(handle, 1):
            if line.strip() == "##FASTA":
                break
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            context = "{}:{}".format(path, number)
            if len(fields) != 9:
                raise ValueError("{}: expected 9 GFF3 columns".format(context))
            if fields[2] != "transcript":
                continue
            attributes = {}
            for item in fields[8].split(";"):
                if not item:
                    continue
                key, separator, value = item.partition("=")
                if not separator or not key or key in attributes:
                    raise ValueError("{}: malformed or duplicate attribute".format(context))
                attributes[key] = value
            transcript_id = attributes.get("transcript_id", "")
            if not transcript_id.strip() or transcript_id == ".":
                raise ValueError("{}: missing transcript_id".format(context))
            if transcript_id in transcripts:
                raise ValueError("{}: duplicate transcript_id {}".format(context, transcript_id))
            try:
                start, end = int(fields[3]), int(fields[4])
            except ValueError:
                raise ValueError("{}: transcript coordinates must be integers".format(context)) from None
            if start < 1 or end < start:
                raise ValueError("{}: invalid transcript coordinates".format(context))
            if not fields[0] or fields[0] == "." or fields[6] not in ("+", "-", ".", "?"):
                raise ValueError("{}: invalid chromosome or strand".format(context))
            transcripts[transcript_id] = dict(chromosome=fields[0], strand=fields[6],
                                               start=start, end=end, length=end - start + 1)
    if not transcripts:
        raise ValueError("{}: no transcript records found".format(path))
    return transcripts


def compare_transcripts(current, extended):
    for transcript_id in sorted(current.keys() | extended.keys()):
        raw, new = current.get(transcript_id), extended.get(transcript_id)
        row = {"transcript_id": transcript_id}
        for annotation, record in (("current", raw), ("extended", new)):
            if record is not None:
                row.update({"{}_{}".format(key, annotation): value for key, value in record.items()})
        if raw is None:
            row["status"] = "missing_from_raw"
        elif new is None:
            row["status"] = "missing_from_extended"
        else:
            row["chromosome_mismatch"] = raw["chromosome"] != new["chromosome"]
            row["strand_mismatch"] = raw["strand"] != new["strand"]
            if row["chromosome_mismatch"] or row["strand_mismatch"]:
                row["status"] = "location_mismatch"
            else:
                row["shorter_left"] = new["start"] > raw["start"]
                row["shorter_right"] = new["end"] < raw["end"]
                if row["shorter_left"] or row["shorter_right"]:
                    row["status"] = "shortened"
                elif (raw["start"], raw["end"]) == (new["start"], new["end"]):
                    row["status"] = "unchanged"
                else:
                    row["status"] = "extended"
        yield row


def main():
    args = get_args()
    try:
        output = Path(args.output)
        for source in (args.current_annotation, args.extended_annotation):
            if output.resolve() == Path(source).resolve() or (output.exists() and output.samefile(source)):
                raise ValueError("output must not overwrite an input annotation")
        # Read both inputs before opening the report, so invalid inputs leave no partial report.
        current = read_transcripts(args.current_annotation)
        extended = read_transcripts(args.extended_annotation)
        counts = Counter()
        with open(args.output, "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=FIELDS)
            writer.writeheader()
            for row in compare_transcripts(current, extended):
                writer.writerow(row)
                counts[row["status"]] += 1
        failed = counts["shortened"] + counts["location_mismatch"]
        print("{}: {} failing transcript(s)".format("FAIL" if failed else "PASS", failed))
        for status in STATUSES:
            print("{}: {}".format(status, counts[status]))
        return 1 if failed else 0
    except (OSError, EOFError, UnicodeError, ValueError, csv.Error) as error:
        print("ERROR: {}".format(error), file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
