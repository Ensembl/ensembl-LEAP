#!/usr/bin/env python 
import pandas as pd
import gzip
import argparse

def get_args():
    # Parse command-line arguments for input/output file paths
    parser = argparse.ArgumentParser(description = "Filter spreadsheet based on final extensions (External vs Internal)")
    parser.add_argument("final_gff", help="final_gff_with_removed_extensions")
    parser.add_argument("spreadsheet", help="Spreadsheet to filter")
    parser.add_argument("output", help="Output file name")
    return parser.parse_args()

def open_file(file_path):
    # Open file normally or with gzip depending on extension
    if file_path.endswith('.gz'):
        return gzip.open(file_path, 'rt')
    else:
        return open(file_path, 'r')

def main():
    args = get_args()
    final_gff = args.final_gff
    spreadsheet = args.spreadsheet
    output = args.output

    # List to store transcript info from GFF
    extended_transcripts_ids = []
    with open_file(final_gff) as f:
        for line in f:
            if line.startswith('#'):
                # Skip comment/header lines in GFF
                continue
            fields = line.strip().split('\t')
            # Parse attributes column into a dictionary
            attr_dict = dict(item.split('=') for item in fields[8].split(';') if '=' in item)
            if fields[2] != "transcript":
                # Only process transcript features
                continue
            # Collect chromosome, transcript_id, and leap_tag for each transcript
            extended_transcripts_ids.append({
                "chromosome": fields[0], 
                "transcript_id": attr_dict.get("transcript_id", ""), 
                "leap_tag": attr_dict.get("leap_tag", "")
            })

    # Convert list of dicts to DataFrame
    extended_transcripts_ids_df = pd.DataFrame(extended_transcripts_ids)
    # Filter out transcripts with leap_tag == "LEAP_none"
    extended_transcripts_ids_df = extended_transcripts_ids_df[extended_transcripts_ids_df['leap_tag'] != "LEAP_none"]
    
    # Split leap_tag string into list (e.g., "LEAP_5,LEAP_3" -> ["LEAP_5", "LEAP_3"])
    extended_transcripts_ids_df['leap_tag'] = extended_transcripts_ids_df['leap_tag'].apply(lambda x: x.split(',') if x else [])

    # Read input spreadsheet (tab-separated)
    input_spreadsheet = pd.read_csv(spreadsheet, sep='\t')
    
    # Extract chromosome from transcript_location column (format: chr:start-end)
    input_spreadsheet['chromosome'] = input_spreadsheet['transcript_location'].str.split(":").str[0]

    # Merge spreadsheet with GFF transcript info on chromosome and transcript_id
    merged_dataframe = pd.merge(
        input_spreadsheet, 
        extended_transcripts_ids_df, 
        left_on=['chromosome', 'transcript_id'], 
        right_on=['chromosome', 'transcript_id'], 
        how='inner'
    )

    # For transcripts lacking LEAP_5, set five_prime_extension/support to 0/'none'
    merged_dataframe.loc[
        merged_dataframe['leap_tag'].apply(lambda tags: 'LEAP_5' not in tags), 
        'five_prime_extension'
    ] = 0
    merged_dataframe.loc[
        merged_dataframe['leap_tag'].apply(lambda tags: 'LEAP_5' not in tags), 
        'five_prime_support'
    ] = 'none'
    # For transcripts lacking LEAP_3, set three_prime_extension/support to 0/'none'
    merged_dataframe.loc[
        merged_dataframe['leap_tag'].apply(lambda tags: 'LEAP_3' not in tags), 
        'three_prime_extension'
    ] = 0
    merged_dataframe.loc[
        merged_dataframe['leap_tag'].apply(lambda tags: 'LEAP_3' not in tags), 
        'three_prime_support'
    ] = 'none'

    # Remove helper columns before saving
    merged_dataframe = merged_dataframe.drop(columns=['chromosome', 'leap_tag'])

    # Write filtered/modified dataframe to output file (tab-separated)
    merged_dataframe.to_csv(output, sep='\t', index=False)

if __name__ == "__main__":
    main()