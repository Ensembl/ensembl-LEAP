#!/usr/bin/env python 
import argparse
import gzip
def get_args():
    parser = argparse.ArgumentParser(description="Add gene and transcript IDs to splice sites in a BED file based on a GFF annotation file.")
    parser.add_argument("current_annotation", help="Path to the current annotation GFF file.")
    parser.add_argument("possible_extensions", help="Path to the possible extensions BED file.")
    parser.add_argument("output_prefix", help="Prefix for the output BED file.")
    return parser.parse_args()

def open_file(file_path):
    if file_path.endswith('.gz'):
        return gzip.open(file_path, 'rt')
    else:
        return open(file_path, 'r')

def write_bed_with_comma_list(output_path, data):
    with open(output_path, 'w') as out:
        for row in data:
            # Convert the last element (a list) to a comma-separated string
            row_out = row[:-1] + [",".join(row[-1])]
            out.write('\t'.join(map(str, row_out)) + '\n')

def main():
    # Parse command line arguments
    args = get_args()
    input_current_annotation = args.current_annotation
    input_possible_extensions = args.possible_extensions
    output_prefix = args.output_prefix
    
    # Determine strand and feature based on the input file names
    if "negative" in input_possible_extensions: 
        strand = "negative"
    else:
        strand = "positive"
    
    if "_FE_" in input_possible_extensions:
        feature = "FE"
    else:
        feature = "LE"

    # Determine the splice site position based on strand and feature
    if (strand == "positive" and feature == "FE") or (strand == "negative" and feature == "LE"):
        splice_site_position = 4
    else:
        splice_site_position = 3

    # Read the current annotation file and create a dictionary mapping splice sites to gene and transcript IDs
    id_dictionary = {}
    with open_file(input_current_annotation) as f:
        for line in f:
            fields = line.strip().split('\t')
            attributes = fields[8]
            label = f'{fields[0]}_{fields[splice_site_position]}_{fields[6]}'
            attr_dict = dict(item.split('=', 1) for item in attributes.split(';') if '=' in item)
            if label not in id_dictionary:
                id_dictionary[label] = {"gene_id":[attr_dict['gene_id']]}
            else:
                if attr_dict['gene_id'] not in id_dictionary[label]["gene_id"]:
                    id_dictionary[label]["gene_id"].append(attr_dict['gene_id'])

    # Read the possible extensions file and add gene and transcript IDs to each splice site
    splice_site_with_gene_id = []
    with open_file(input_possible_extensions) as f:
        for line in f:
            fields = line.strip().split("\t")
            label = f'{fields[0]}_{fields[splice_site_position - 2]}_{fields[5]}'
            if label in id_dictionary:
                fields.append(id_dictionary[label]["gene_id"])
                splice_site_with_gene_id.append(fields)

    # Write the output BED file with splice sites and associated gene IDs
    output_bed = f"{output_prefix}_GID.bed"
    write_bed_with_comma_list(output_bed, splice_site_with_gene_id)

if __name__ == "__main__":
    main()
