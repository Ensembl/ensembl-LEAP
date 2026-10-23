#!/usr/bin/env python 
import pandas as pd
import argparse
import os

def get_args():
    params = argparse.ArgumentParser(description="Extract chromosome names from a BAM file using samtools idxstats.")
    params.add_argument("reads", type=str, help="Path to the BAM file.")
    params.add_argument("reads_index", type=str, help="Path to the BAM index file (.bai/.csi).")
    return params.parse_args()

def get_chromosome_from_reads(reads, reads_index):
    # Extract the chromosome names from the bam file using idxstats
    os.system(f'samtools idxstats {reads} > idx_stats.txt')
    # Read the idx_stats.txt file into a pandas DataFrame
    idx_stats = pd.read_csv(f'idx_stats.txt', sep="\t", header=None, names=["chromosome", "length", "mapped_reads", "unmapped_reads"])
    # Filter out empty chromosomes
    idx_stats = idx_stats[idx_stats['mapped_reads'] > 0]
    # Filter out the 'MT', "chrM", and "chrMT" chromosomes
    idx_stats = idx_stats[~idx_stats['chromosome'].isin(['MT', 'chrM', 'chrMT', 'J01415.2', 'NC_012920.1'])]
    # Add information of idxstats
    idx_stats['bam'] = os.path.join(os.getcwd(),reads)
    idx_stats['bai'] = os.path.join(os.getcwd(),reads_index)
    # Reorder and select relevant columns
    idx_stats = idx_stats[['bam', 'bai', 'chromosome']]
    # Save sample sheet
    idx_stats.to_csv('chromosome_sample_sheet.csv', index=False)

def main():
    # Get command line arguments
    args = get_args()
    reads = args.reads
    reads_index = args.reads_index
    # Extract chromosome names from the reads
    get_chromosome_from_reads(reads, reads_index)

if __name__ == "__main__":
    main()
