#!/usr/bin/env python 
import argparse
import gzip
import pandas as pd

def get_args():
    parser = argparse.ArgumentParser(description="Extract unique splice sites from GFF file.")
    parser.add_argument("input_path", type=str, help="Path to the input GFF file.")
    parser.add_argument("output_prefix", type=str, help="Prefix for the output files.")
    return parser.parse_args()

def open_file(file_path):
    if file_path.endswith('.gz'):
        return gzip.open(file_path, 'rt')
    else:
        return open(file_path, 'r')

def main():
    # Parse command line arguments
    args = get_args()
    input_path = args.input_path
    output_prefix = args.output_prefix

    # Initialize a dictionary to store unique splice sites
    unique_splice_sites = {}
    
    # Determine the type of splice site based on the input file name and process accordingly
    if ("positive_FE" in input_path) or ("negative_LE" in input_path):
        with open_file(input_path) as input_file:
            for line in input_file:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                splice_site_id = f'{fields[0]}_{fields[4]}_{fields[6]}'
                if splice_site_id not in unique_splice_sites:
                    unique_splice_sites[splice_site_id] = {"non_splice_site": int(fields[3]), "category": "current_annotation"}
                elif unique_splice_sites[splice_site_id]["non_splice_site"] > int(fields[3]): 
                    unique_splice_sites[splice_site_id]["non_splice_site"] = int(fields[3])
    elif ("positive_LE" in input_path) or ("negative_FE" in input_path):
        with open_file(input_path) as input_file:
            for line in input_file:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                splice_site_id = f'{fields[0]}_{fields[3]}_{fields[6]}'
                if splice_site_id not in unique_splice_sites:
                    unique_splice_sites[splice_site_id] = {"non_splice_site": int(fields[4]), "category": "current_annotation"}
                elif unique_splice_sites[splice_site_id]["non_splice_site"] < int(fields[4]): 
                    unique_splice_sites[splice_site_id]["non_splice_site"] = int(fields[4])
    # Convert the dictionary to a DataFrame
    unique_splice_sites_df = pd.DataFrame.from_dict(unique_splice_sites, orient='index').reset_index().rename(columns={'index': 'splice_id'})
    # Save the DataFrame to a CSV file
    unique_splice_sites_df.to_csv(f'{output_prefix}_USS.csv', index=False)

if __name__ == "__main__":
    main()