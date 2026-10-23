#!/usr/bin/env python 
import argparse
import gzip

def get_args():
    parser = argparse.ArgumentParser(description="Extract first and last exon from GFF file for each transcript.")
    parser.add_argument("input_file", help="Input GFF file (can be gzipped).")
    parser.add_argument("output_prefix", help="Output file for first exons.")
    parser.add_argument("-f", "--format", choices=["gff", "gtf"], default="gff", help="Input file format (default: gff).")
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
    input_format = args.format
    output_first_exon = f"{output_prefix}_FE.gff"
    output_last_exon = f"{output_prefix}_LE.gff"
    if input_format == "gtf" and "negative" in input_file:
        output_last_exon = f"{output_prefix}_FE.gff"
        output_first_exon = f"{output_prefix}_LE.gff"
    
    # Extract exon ids for each transcript
    transcript_id = {}
    with open_file(input_file) as input_gff:
        for line in input_gff:
            if line.startswith("#"):
                continue
            fields = line.strip().split('\t')
            if fields[2] == "exon":
                attributes = fields[8]
                if input_format == "gff":
                    attr_dict = dict(item.split('=', 1) for item in attributes.split(';') if '=' in item)
                elif input_format == "gtf":
                    attr_dict = {k.strip(): v.strip().replace('"', '') for k, v in (item.split(' "', 1) for item in attributes.split(';') if ' "' in item)}
                if attr_dict['transcript_id'] not in transcript_id:
                    transcript_id[attr_dict['transcript_id']] = [attr_dict['exon_id']]
                else:
                    transcript_id[attr_dict['transcript_id']].append(attr_dict['exon_id'])

    # Extract first and last exon for each transcript
    result = {k: {"first": v[0], "last": v[-1]} for k, v in transcript_id.items() if v}

    # Remove single exon transcripts
    result = {k: val for k, val in result.items() if val["first"] != val["last"]}

    # Load the input file line by line and write the first and last exons to different files
    with open_file(input_file) as input_gff:
        for line in input_gff:
            if line.startswith("#"):
                continue
            fields = line.strip().split('\t')
            attributes = fields[8]
            if input_format == "gff":
                attr_dict = dict(item.split('=', 1) for item in attributes.split(';') if '=' in item)
            elif input_format == "gtf":
                attr_dict = {k.strip(): v.strip().replace('"', '') for k, v in (item.split(' "', 1) for item in attributes.split(';') if ' "' in item)}
            if fields[2] != "exon":
                continue
            if attr_dict['transcript_id'] in result:
                if attr_dict['exon_id'] == result[attr_dict['transcript_id']]['first']:
                    with open(output_first_exon, 'a') as first_exon_file:
                        first_exon_file.write(line)
                elif attr_dict['exon_id'] == result[attr_dict['transcript_id']]['last']:
                    with open(output_last_exon, 'a') as last_exon_file:
                        last_exon_file.write(line)

if __name__ == "__main__":
    main()
