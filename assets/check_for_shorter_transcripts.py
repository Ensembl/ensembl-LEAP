#!/usr/bin/env python 
import gzip
import pandas as pd
import argparse

def get_args():
    parser = argparse.ArgumentParser(description="This script checks for transcripts which have gotten smaller compared to the original model")
    parser.add_argument("current_annotation", help="Path to the current annotation GFF3 file")
    parser.add_argument("extended_annotation", help="Path to the extended annotation GFF3 file")
    parser.add_argument("output", help="Path to save the output CSV file")
    return parser.parse_args()

def open_file(file_path):
    if file_path.endswith('.gz'):
        return gzip.open(file_path, 'rt')
    else:
        return open(file_path, 'r')

def main():
    args = get_args()
    current_annotation = args.current_annotation
    extended_annotation = args.extended_annotation
    output_file = args.output
    current_annotation_start_ends = []
    with open_file(current_annotation) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            attr_dict = dict(item.split("=", 1) for item in fields[-1].strip().split(";") if item and "=" in item)
            if fields[2] != "transcript":
                continue
            start = int(fields[3])
            end = int(fields[4])
            current_annotation_start_ends.append({"chromosome": fields[0], "strand": fields[6], "transcript_id": attr_dict.get("transcript_id"), "start": start, "end": end})

    current_annotation_start_ends = pd.DataFrame(current_annotation_start_ends)

    extended_annotation_start_ends = []
    with open_file(extended_annotation) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            attr_dict = dict(item.split("=", 1) for item in fields[-1].strip().split(";") if item and "=" in item)
            if fields[2] != "transcript":
                continue
            start = int(fields[3])
            end = int(fields[4])
            extended_annotation_start_ends.append({"chromosome": fields[0], "strand": fields[6], "transcript_id": attr_dict.get("transcript_id"), "start": start, "end": end})

    extended_annotation_start_ends = pd.DataFrame(extended_annotation_start_ends)

    merged_annotation = pd.merge(current_annotation_start_ends, extended_annotation_start_ends, on=["chromosome", "strand", "transcript_id"], suffixes=("_current", "_extended"))

    merged_annotation['start_diff'] = merged_annotation['start_extended'] - merged_annotation['start_current']
    merged_annotation['end_diff'] = merged_annotation['end_extended'] - merged_annotation['end_current']
    
    merged_annotation['shorter_left'] = merged_annotation['start_current'] < merged_annotation['start_extended']
    merged_annotation['shorter_right'] = merged_annotation['end_current'] > merged_annotation['end_extended']

    merged_annotation['label'] = merged_annotation['chromosome'].astype(str) + ":" + merged_annotation['start_current'].astype(str) + "-" + merged_annotation['end_current'].astype(str)
    merged_annotation.to_csv(output_file, index=False)

# bad_difference = merged_annotation[(merged_annotation['start_diff'] > 0) | (merged_annotation['end_diff'] < 0)]
if __name__ == "__main__":
    main()