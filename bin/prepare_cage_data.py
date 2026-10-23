#!/usr/bin/env python 
import pandas as pd
import argparse
import os
import gzip

def get_args():
    parser = argparse.ArgumentParser(description="Prepare CAGE data from GFF file.")
    parser.add_argument('input_file', type=str, help='Path to the input GFF file.')
    parser.add_argument("read_sample_sheet", type=str, help="Path to the read sample sheet.")
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
    read_sample_sheet = args.read_sample_sheet
    # Extract chromosomes from read sample sheet
    read_chromosomes = pd.read_csv(read_sample_sheet)['chromosome'].astype(str).str.strip()
    # Process the GFF file and create separate files for each chromosome and strand
    file_names = []
    chromosomes = []
    strands = []
    exons = []
    with open_file(input_file) as input_gff:
        for line in input_gff:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            chromosome = fields[0]
            if chromosome.startswith('chr'):
                chromosome = chromosome[3:]
                fields[0] = chromosome  # Update the chromosome in the fields
            if input_file.endswith('.bed'):
                extension = 'bed'
                strand = fields[5]
                fields[1] = str(int(fields[1]) + 1)  # Adjust start position for BED format
            else:
                strand = fields[6]
                extension = 'gff'
            exon = "FE"
            if strand == '+':
                strand = "positive"
            else:
                strand = "negative"
            write_file = f'CAGE_{chromosome}_{strand}.{extension}'
            if write_file not in file_names:
                file_names.append(write_file)
                chromosomes.append(chromosome)
                strands.append(strand)
                exons.append(exon)
            with open(write_file, 'a') as output_gff:
                output_gff.write("\t".join(fields) + '\n')
    # Add the current working directory to the file names
    file_names = [os.path.join(os.getcwd(), fname) for fname in file_names]
    # Filter the DataFrame to only include rows where chromosome is in read_chromosomes
    df = pd.DataFrame({'chromosome': chromosomes, 'file_name': file_names, 'strand': strands, 'exon': exons})
    df['chromosome'] = df['chromosome'].astype(str).str.strip()
    df = df[df['chromosome'].isin(read_chromosomes.values)]
    # Create a sample sheet with the chromosome, file name, and strand
    df.to_csv('CAGE_sample_sheet.csv', index=False)

if __name__ == "__main__":
    main()
