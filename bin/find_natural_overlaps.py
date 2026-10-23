#!/usr/bin/env python 
import argparse
import os
import gzip

def get_args():
    parser = argparse.ArgumentParser(description="Convert GFF files to BED format and find overlapping gene IDs.")
    parser.add_argument('input_fe', type=str, help='Input GFF file for FE')
    parser.add_argument('input_le', type=str, help='Input GFF file for LE')
    parser.add_argument('output_prefix', type=str, help='Prefix for output files')
    return parser.parse_args()

def open_file(file_path):
    if file_path.endswith('.gz'):
        return gzip.open(file_path, 'rt')
    else:
        return open(file_path, 'r')

def main():
    # Parse command line arguments
    args = get_args()
    input_fe = args.input_fe
    input_le = args.input_le
    output_prefix = args.output_prefix

    # Convert GFF files to BED format
    with open_file(input_fe) as fe, open_file(input_le) as le:
        with open('FE.bed', 'w') as output_file:
            for line in fe:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                attr_dict = {k: v for k, v in (x.split('=') for x in fields[8].split(';') if '=' in x)}
                output_line = f"{fields[0]}\t{fields[3]}\t{fields[4]}\t{attr_dict['gene_id']}\t0\t{fields[6]}\n"
                output_file.write(output_line)
        with open('LE.bed', 'w') as output_file:
            for line in le:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                attr_dict = {k: v for k, v in (x.split('=') for x in fields[8].split(';') if '=' in x)}
                output_line = f"{fields[0]}\t{fields[3]}\t{fields[4]}\t{attr_dict['gene_id']}\t0\t{fields[6]}\n"
                output_file.write(output_line)

    # Overlap the two BED files using bedtools
    os.system('bedtools intersect -a FE.bed -b LE.bed -wa -wb > FE_LE.bed')
    
    # If the output file is empty create an empty file
    if not os.path.exists('FE_LE.bed') or os.path.getsize('FE_LE.bed') == 0:
        os.system(f'touch {output_prefix}_NO.txt')
    else:
        # Extract the gene IDs from the overlapping regions
        with open('FE_LE.bed', 'r') as input_file, open('FE_LE_gene_ids.txt', 'w') as output_file:
            for line in input_file:
                fields = line.strip().split('\t')
                fe_gene_id = fields[3]
                le_gene_id = fields[9]
                output_file.write(f"{fe_gene_id}\t{le_gene_id}\n")
        # Remove duplicate gene ids
        os.system(f'uniq FE_LE_gene_ids.txt > {output_prefix}_NO.txt')
        # Remove the temporary file 
        os.remove('FE_LE_gene_ids.txt')
    
    # Remove the temporary files
    os.remove('FE.bed')
    os.remove('LE.bed')
    os.remove('FE_LE.bed')

if __name__ == "__main__":
    main()
    
    
