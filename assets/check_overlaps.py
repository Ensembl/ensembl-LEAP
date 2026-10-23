#!/usr/bin/env python 
import os
import pandas as pd
import numpy as np
import argparse
import gzip

def get_args():
    parser = argparse.ArgumentParser(description="Compare overlaps in the current annotation to the extended annotation")
    parser.add_argument('temp_output_dir', type=str, help='Path to temporary output directory')
    parser.add_argument('output_dir', type=str, help='Path to output directory')
    parser.add_argument('current_annotation', type=str, help='Path to current annotation gff3')
    parser.add_argument('extended_annotation', type=str, help='Path to extended annotation gff3')
    return parser.parse_args()

def open_gff(gff_path):
    if gff_path.endswith('.gz'):
        return gzip.open(gff_path, 'rt')
    else:
        return open(gff_path, 'r')

def main():
    # Get command line arguments
    args = get_args()
    temp_output_dir = args.temp_output_dir
    output_dir = args.output_dir

    # Path to Temporary output directory
    os.makedirs(temp_output_dir, exist_ok=True)

    # Path to output directory
    os.makedirs(output_dir, exist_ok=True)

    # Provided gffs
    current_annotation_gff = args.current_annotation
    extended_annotation_gff = args.extended_annotation

    # Convert GFF to BED format for current annotation
    input_gff = current_annotation_gff

    with open_gff(input_gff) as f, open(os.path.join(temp_output_dir, 'current_annotation_transcripts.bed'), 'w') as of:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if fields[2] !='transcript':
                continue
            attr_dict = dict(item.split("=") for item in fields[-1].split(";") if "=" in item)
            if attr_dict.get("transcript_type") != "protein_coding":
                continue
            # Write to BED format
            of.write('\t'.join([fields[0],fields[3],fields[4],attr_dict['gene_id'],attr_dict['transcript_id'],fields[6]]) + '\n')

    # Convert GFF to BED format for extended annotation
    input_gff = extended_annotation_gff

    with open_gff(input_gff) as f, open(os.path.join(temp_output_dir, 'extended_annotation_transcripts.bed'), 'w') as of:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if fields[2] !='transcript':
                continue
            attr_dict = dict(item.split("=") for item in fields[-1].split(";") if "=" in item)
            if attr_dict.get("transcript_type") != "protein_coding":
                continue
            of.write('\t'.join([fields[0],fields[3],fields[4],attr_dict['gene_id'],attr_dict['transcript_id'],fields[6]]) + '\n')
    
    # Sort and find overlaps
    os.system(f"sort -k1,1V -k2,2n {os.path.join(temp_output_dir, 'current_annotation_transcripts.bed')} > {os.path.join(temp_output_dir, 'current_annotation_transcripts.sorted.bed')}")
    os.system(f"sort -k1,1V -k2,2n {os.path.join(temp_output_dir, 'extended_annotation_transcripts.bed')} > {os.path.join(temp_output_dir, 'extended_annotation_transcripts.sorted.bed')}")
    os.system(f"bedtools intersect -a {os.path.join(temp_output_dir, 'current_annotation_transcripts.sorted.bed')} -b {os.path.join(temp_output_dir, 'current_annotation_transcripts.sorted.bed')} -s -wa -wb > {os.path.join(temp_output_dir, 'current_annotation_transcripts.overlaps.bed')}")
    os.system(f"bedtools intersect -a {os.path.join(temp_output_dir, 'extended_annotation_transcripts.sorted.bed')} -b {os.path.join(temp_output_dir, 'extended_annotation_transcripts.sorted.bed')} -s -wa -wb > {os.path.join(temp_output_dir, 'extended_annotation_transcripts.overlaps.bed')}")

    # Process current annotation overlaps
    current_annotation_between_gene_overlaps = []
    with open(os.path.join(temp_output_dir, 'current_annotation_transcripts.overlaps.bed'), 'r') as f:
        for line in f:
            fields = line.strip().split("\t")
            if fields[3] != fields[9]:
                # Sort genes to ensure consistent ordering (avoid X,Y vs Y,X duplicates)
                gene_pair = sorted([fields[3], fields[9]])
                current_annotation_between_gene_overlaps.append({"CA_gene_1": gene_pair[0], "CA_gene_2": gene_pair[1]})

    current_annotation_between_gene_overlaps_df = pd.DataFrame(current_annotation_between_gene_overlaps).drop_duplicates()

    # Process extended annotation overlaps
    extended_annotation_between_gene_overlaps = []
    with open(os.path.join(temp_output_dir, 'extended_annotation_transcripts.overlaps.bed'), 'r') as f:
        for line in f:
            fields = line.strip().split("\t")
            if fields[3] != fields[9]:
                # Sort genes to ensure consistent ordering (avoid X,Y vs Y,X duplicates)
                gene_pair = sorted([fields[3], fields[9]])
                extended_annotation_between_gene_overlaps.append({"EA_gene_1": gene_pair[0], "EA_gene_2": gene_pair[1]})

    extended_annotation_between_gene_overlaps_df = pd.DataFrame(extended_annotation_between_gene_overlaps).drop_duplicates()

    # Merge current and extended annotations
    merged_annotations = pd.merge(current_annotation_between_gene_overlaps_df, extended_annotation_between_gene_overlaps_df, left_on=['CA_gene_1', 'CA_gene_2'], right_on=['EA_gene_1', 'EA_gene_2'], how='outer')
    
    # Add columns to indicate if the gene pairs are the same in both annotations
    merged_annotations["same_gene_pair"] = merged_annotations.apply(
        lambda row: (
            (pd.isna(row['EA_gene_1']) and pd.isna(row['EA_gene_2']))
            or (row['CA_gene_1'] == row['EA_gene_1'] and row['CA_gene_2'] == row['EA_gene_2'])
        ),
        axis=1
    )
    
    # Output the merged annotations
    merged_annotations[merged_annotations['same_gene_pair'] == False].to_csv(os.path.join(output_dir, 'overlaps_different_gene_pairs.tsv'), sep = "\t", index=False)
    merged_annotations[merged_annotations['same_gene_pair'] == True].to_csv(os.path.join(output_dir, 'overlaps_same_gene_pairs.tsv'), sep = "\t", index=False)
    
    # Clean up temporary files
    os.system(f'rm -r {temp_output_dir}')

if __name__ == "__main__":
    main()