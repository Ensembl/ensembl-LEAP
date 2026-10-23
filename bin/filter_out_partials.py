#!/usr/bin/env python 
import pandas as pd
import argparse
import gzip

def get_args():
    parser = argparse.ArgumentParser(description="Filter a GFF file to retain transcripts that have both five prime and three prime UTRs.")
    parser.add_argument("input_file", help="Path to the input GFF file")
    parser.add_argument("output_prefix", help="Path to the output GFF file", default="filtered")
    return parser.parse_args()

def open_file(file_path):
    if file_path.endswith('.gz'):
        return gzip.open(file_path, 'rt')
    else:
        return open(file_path, 'r')

def count_utr(input_file):
    # Prepare dictionary to store transcripts
    transcript_ids = {}
    non_protein_coding = set()
    # Open file line by line
    with open_file(input_file) as input_gff:
        for line in input_gff:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            attributes = fields[8]
            attr_dict = dict(item.split('=', 1) for item in attributes.split(';') if '=' in item)
            if 'transcript_type' in attr_dict:
                if attr_dict['transcript_type'] != 'protein_coding':
                    # if attr_dict['gene_type'] != 'protein_coding':
                    non_protein_coding.add(attr_dict['transcript_id'])
                    continue
                if attr_dict['transcript_id'] not in transcript_ids:
                    transcript_ids[attr_dict['transcript_id']] = {'five_prime_UTR': 0, 'three_prime_UTR': 0}
                if fields[2] == 'five_prime_UTR':
                    transcript_ids[attr_dict['transcript_id']]['five_prime_UTR'] += 1
                elif fields[2] == 'three_prime_UTR':
                    transcript_ids[attr_dict['transcript_id']]['three_prime_UTR'] += 1
    return transcript_ids, non_protein_coding

def filter_gff(input_file, output_file, full_transcript_ids):
    with open_file(input_file) as input_gff, open(output_file, 'w') as output_gff:
        for line in input_gff:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            attributes = fields[8]
            attr_dict = dict(item.split('=', 1) for item in attributes.split(';') if '=' in item)
            if 'transcript_id' in attr_dict:
                if attr_dict['transcript_id'] in full_transcript_ids:
                    output_gff.write(line)

def main():
    # Parse command line arguments
    args = get_args()
    input_file = args.input_file
    output_file = f'{args.output_prefix}.gff'
    
    # Count UTRs in the input GFF file
    # utr_counts, non_protein_coding_tids = pd.DataFrame.from_dict(count_utr(input_file), orient='index').reset_index().rename(columns={'index': 'transcript_id'})
    utr_counts, non_protein_coding_tids = count_utr(input_file)
    utr_counts = pd.DataFrame.from_dict(utr_counts, orient='index').reset_index().rename(columns={'index': 'transcript_id'})
    
    # Filter transcript id's based on UTR counts
    full_transcripts_ids = utr_counts[(utr_counts['five_prime_UTR'] > 0) & (utr_counts['three_prime_UTR'] > 0)].set_index('transcript_id').to_dict(orient='index')
    partial_transcripts_ids = utr_counts[(utr_counts['five_prime_UTR'] == 0) | (utr_counts['three_prime_UTR'] == 0)].set_index('transcript_id').to_dict(orient='index')

    # Filter the GFF file based on the full transcripts
    filter_gff(input_file, output_file.replace(".gff", "_FT.gff"), full_transcripts_ids)
    
    # Filter the GFF file for partial transcripts
    filter_gff(input_file, output_file.replace(".gff", "_PT.gff"), partial_transcripts_ids)
    
    # Filter the GFF file for non-protein coding transcripts
    filter_gff(input_file, output_file.replace(".gff", "_NPC.gff"), non_protein_coding_tids)

if __name__ == "__main__":
    main()