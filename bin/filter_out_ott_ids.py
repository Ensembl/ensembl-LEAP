#!/usr/bin/env python 
import argparse

def get_args():
    parser = argparse.ArgumentParser(description="Filter out OTT transcript and genes while leaving OTT exons")
    parser.add_argument("input_gff", type=str, help="Input GFF to be filtered")
    parser.add_argument("output_gff", type=str, help="Output GFF name")
    return parser.parse_args()

def main():
    args = get_args()
    input_path = args.input_gff
    output_path = args.output_gff
    ott_ids = set()
    with open(input_path, "r") as input_gff:
        for line in input_gff:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if fields[2] not in ['gene', 'transcript']:
                continue
            attr_dict = dict(item.split("=") for item in [x for x in fields[-1].strip().split(";") if x])
            if fields[2] == "gene" and attr_dict['gene_id'].startswith("OTT"):
                ott_ids.add(attr_dict['gene_id'])
            if fields[2] == "transcript" and attr_dict['transcript_id'].startswith("OTT"):
                ott_ids.add(attr_dict['transcript_id'])

    with open(input_path, "r") as input_gff, open(output_path, "w") as output_gff:
        for line in input_gff:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            attr_dict = dict(item.split("=") for item in [x for x in fields[-1].strip().split(";") if x])
            if fields[2] == "gene":
                if attr_dict['gene_id'] in ott_ids:
                    continue
                else:
                    output_gff.write(line)
                    continue
            # Transcript parents are gene IDs, so also check the transcript's
            # own ID to remove OTT transcripts belonging to Ensembl genes.
            if (attr_dict.get('transcript_id') not in ott_ids
                    and attr_dict['Parent'] not in ott_ids
                    and attr_dict['gene_id'] not in ott_ids):
                output_gff.write(line)

if __name__ == "__main__":
    main()

