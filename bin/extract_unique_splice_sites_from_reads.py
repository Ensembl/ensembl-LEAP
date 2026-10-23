#!/usr/bin/env python 
import argparse
import gzip
import pandas as pd

def get_args():
    parser = argparse.ArgumentParser(description="Extract unique splice sites from GFF file.")
    parser.add_argument("input_path", type=str, help="Path to the input GFF file.")
    parser.add_argument("output_prefix", type=str, help="Prefix for the output files.")
    parser.add_argument("min_reads", type=int, help="Minimum number of reads to consider an extension.")
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
    min_reads = args.min_reads

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
                    unique_splice_sites[splice_site_id] = {"non_splice_site": [fields[3]], "category": "reads"}
                else:
                    unique_splice_sites[splice_site_id]["non_splice_site"].append(fields[3])
        # Sort and filter out lists with fewer than N elements
        to_remove = []
        for key, site in unique_splice_sites.items():
            if isinstance(site["non_splice_site"], list):
                if len(site["non_splice_site"]) >= min_reads:
                    site["non_splice_site"].sort()
                else:
                    to_remove.append(key)
        for key in to_remove:
            del unique_splice_sites[key]
        # Limit the extension to regions covered with at least N reads
        for key, site in unique_splice_sites.items():
            if isinstance(site["non_splice_site"], list):
                site["non_splice_site"] = site["non_splice_site"][min_reads - 1]
    elif ("positive_LE" in input_path) or ("negative_FE", input_path):
        with open_file(input_path) as input_file:
            for line in input_file:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                splice_site_id = f'{fields[0]}_{fields[3]}_{fields[6]}'
                if splice_site_id not in unique_splice_sites:
                    unique_splice_sites[splice_site_id] = {"non_splice_site": [fields[4]], "category": "reads"}
                else:
                    unique_splice_sites[splice_site_id]["non_splice_site"].append(fields[4])
        # Sort and filter out lists with fewer than N elements
        to_remove = []
        for key, site in unique_splice_sites.items():
            if isinstance(site["non_splice_site"], list):
                if len(site["non_splice_site"]) >= min_reads:
                    site["non_splice_site"].sort()
                else:
                    to_remove.append(key)
        for key in to_remove:
            del unique_splice_sites[key]
        for key, site in unique_splice_sites.items():
            if isinstance(site["non_splice_site"], list):
                site["non_splice_site"] = site["non_splice_site"][-min_reads]
    # Convert the dictionary to a DataFrame
    unique_splice_sites_df = pd.DataFrame.from_dict(unique_splice_sites, orient='index').reset_index().rename(columns={'index': 'splice_id'})
    # Save the DataFrame to a CSV file
    unique_splice_sites_df.to_csv(f'{output_prefix}_USS.csv', index=False)

if __name__ == "__main__":
    main()

