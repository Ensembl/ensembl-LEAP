#!/usr/bin/env python 
import argparse
import gzip
import csv

def get_args():
    parser = argparse.ArgumentParser(description="Process BED files to extract longest possible extension and add shorter extensions as an atribute")
    parser.add_argument("input_file", help="Input BED file to process")
    parser.add_argument("output_prefix", help="Output file prefix")
    return parser.parse_args()

def open_file(file_path):
    if file_path.endswith('.gz'):
        return gzip.open(file_path, 'rt')
    else:
        return open(file_path, 'r')

def main():
    # Parse command line arguments
    args = get_args()
    input_file = args.input_file
    output_prefix = args.output_prefix

    # Check for strand and site in the input files
    if "negative" in input_file:
        strand = "negative"
    else:
        strand = "positive"

    if "_FE_" in input_file:
        site = "FE"
    else:
        site = "LE"

    if (strand == "positive" and site == "FE") or (strand == "negative" and site == "LE"):
        splice_site_position = 2
        other_position = 1
    else:
        splice_site_position = 1
        other_position = 2

    # Read the input file and create a dictionary to store splice site information
    splice_site_dict = {}
    with open_file(input_file) as f:
        for line in f:
            fields = line.strip().split('\t')
            if fields[splice_site_position] not in splice_site_dict:
                if fields[3] == "current_annotation":
                    extension_source_index = 0
                else:
                    extension_source_index = -1
                splice_site_dict[fields[splice_site_position]] = {"chromosome":fields[0],
                                                "start":int(fields[1]),
                                                "end":int(fields[2]),
                                                "strand":fields[5],
                                                "extension_source": extension_source_index,
                                                "shorter_extensions": [],
                                                "read_count": 0}
            elif splice_site_dict[fields[splice_site_position]]["extension_source"] == -1:
                splice_site_dict[fields[splice_site_position]]["shorter_extensions"].append(int(fields[other_position]))
            # if fields[3] != "current_annotation" and splice_site_dict[fields[splice_site_position]]["extension_source"] == -1:
            if fields[3] != "current_annotation":
                splice_site_dict[fields[splice_site_position]]["read_count"] += 1
            # if fields[3] == "current_annotation":
            if fields[3] == "current_annotation" and not (splice_site_dict[fields[splice_site_position]]["extension_source"] > -1):
                splice_site_dict[fields[splice_site_position]]["extension_source"] = splice_site_dict[fields[splice_site_position]]["read_count"]

    # Write splice_site_dict to BED file
    with open(f'{output_prefix}_EXT.bed', 'w') as bedfile:
        for key, value in splice_site_dict.items():
            # Use shorter_extensions as the name (joined by '/')
            name = '/'.join(map(str, value['shorter_extensions']))
            if name == "":
                # name = "current_annotation"
                name = "."
            # Use extension_source as the score
            score = value['extension_source']
            # Use read_count as itemRgb (as a string, e.g., "read_count,0,0")
            itemrgb = value['read_count']
            bed_fields = [
                value['chromosome'],
                str(value['start']),
                str(value['end']),
                name,
                str(score),
                value['strand'],
                ".",
                ".",
                str(itemrgb)
            ]
            bedfile.write('\t'.join(bed_fields) + '\n')

if __name__ == "__main__":
    main()