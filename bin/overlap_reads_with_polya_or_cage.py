#!/usr/bin/env python 
import pandas as pd
import argparse
import gzip

def get_args():
    parser = argparse.ArgumentParser(description="Overlap unique splice sites in reasds with polyA or CAGE features to create BED files.")
    parser.add_argument("input_reads", type=str, help="Input file with unique splice sites in reads.")
    parser.add_argument("input_feature", type=str, help="Input BED file with polyA or CAGE features.")
    parser.add_argument("output_prefix", type=str, help="Prefix for output BED file.")
    return parser.parse_args()

def open_file(file_path):
    if file_path.endswith('.gz'):
        return gzip.open(file_path, 'rt')
    else:
        return open(file_path, 'r')

def main():
    # Parse command line arguments
    args = get_args()
    input_reads = args.input_reads
    input_feature = args.input_feature
    output_prefix = args.output_prefix

    # Determine strand and evidence type from input files and read the feature file into a DataFrame
    if "negative" in input_reads:
        strand = "negative"
    elif "positive" in input_reads:
        strand = "positive"

    if "polya" in input_feature.lower():
        if input_feature.endswith('.bed'):
            use_cols = [0, 1, 2]
        else:
            use_cols = [0, 3, 4]
        evidence = "polya"
        feature_df = pd.read_csv(input_feature, sep="\t", names=['chromosome','start','end'], usecols=use_cols)
        feature_df['start'] = feature_df['start'].astype(int)
        feature_df['end'] = feature_df['end'].astype(int)
    elif "cage" in input_feature.lower():
        if input_feature.endswith('.bed'):
            use_cols = [0, 1, 2]
        else:
            use_cols = [0, 3, 4]
        evidence = "cage"
        feature_df = pd.read_csv(input_feature, sep="\t", names=['chromosome','start','end'], usecols=use_cols)
        feature_df['start'] = feature_df['start'].astype(int)
        feature_df['end'] = feature_df['end'].astype(int)
        
    
    # Create empty output file
    open(f'{output_prefix}_all.bed', "w").close()

    # Process the input reads based on strand and evidence type
    if (strand == "positive" and evidence == "cage") or (strand == "negative" and evidence == "polya"):
        with open_file(input_reads) as unique_splice_sites:
            header = next(unique_splice_sites)  # Skip header line
            for idx, line in enumerate(unique_splice_sites, start=2):
                fields = line.strip().split(",")
                info = fields[0].split("_")
                read_row = {
                    "splice_id": fields[0],
                    'chromosome': info[0],
                    'non_splice_site': int(fields[1]),
                    'splice_site': int(info[1]),
                    'strand': info[2],
                    'category': fields[2]
                }
                matching_evidence = feature_df[(feature_df['start'] >= read_row['non_splice_site']) & ((feature_df['end']) <= read_row['splice_site'])]
                if not matching_evidence.empty:
                    if len(matching_evidence) > 1:
                        for site in matching_evidence['start'].sort_values().tolist():
                            bed_formatted_line = f'{read_row["chromosome"]}\t{site}\t{read_row["splice_site"]}\t{read_row["category"]}\t0\t{read_row["strand"]}\n'
                            with open(f'{output_prefix}_all.bed', "a") as out_f:
                                out_f.write(bed_formatted_line)
                    else:
                        site = matching_evidence['start'].sort_values().tolist()[0]
                        bed_formatted_line = f'{read_row["chromosome"]}\t{site}\t{read_row["splice_site"]}\t{read_row["category"]}\t0\t{read_row["strand"]}\n'
                        with open(f'{output_prefix}_all.bed', "a") as out_f:
                            out_f.write(bed_formatted_line)
    elif (strand == "positive" and evidence == "polya") or (strand == "negative" and evidence == "cage"):
        with open_file(input_reads) as unique_splice_sites:
            header = next(unique_splice_sites)  # Skip header line
            for idx, line in enumerate(unique_splice_sites, start=2):
                fields = line.strip().split(",")
                info = fields[0].split("_")
                read_row = {
                    "splice_id": fields[0],
                    'chromosome': info[0],
                    'non_splice_site': int(fields[1]),
                    'splice_site': int(info[1]),
                    'strand': info[2],
                    'category': fields[2]
                }
                matching_evidence = feature_df[(feature_df['start'] >= read_row['splice_site']) & ((feature_df['end']) <= read_row['non_splice_site'])]
                if not matching_evidence.empty:
                    if len(matching_evidence) > 1:
                        for site in matching_evidence['end'].sort_values().tolist():
                            bed_formatted_line = f'{read_row["chromosome"]}\t{read_row["splice_site"]}\t{site}\t{read_row["category"]}\t0\t{read_row["strand"]}\n'
                            with open(f'{output_prefix}_all.bed', "a") as out_f:
                                out_f.write(bed_formatted_line)
                    else:
                        site = matching_evidence['end'].sort_values().tolist()[0]
                        bed_formatted_line = f'{read_row["chromosome"]}\t{read_row["splice_site"]}\t{site}\t{read_row["category"]}\t0\t{read_row["strand"]}\n'
                        with open(f'{output_prefix}_all.bed', "a") as out_f:
                            out_f.write(bed_formatted_line)

if __name__ == "__main__":
    main()
    