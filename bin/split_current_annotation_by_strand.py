#!/usr/bin/env python 
import argparse

def get_args():
    parser = argparse.ArgumentParser(description="A script to split a gff file by strand")
    parser.add_argument('input_file', type=str, help='Input GFF')
    parser.add_argument('output_prefix', type=str, help='Output prefix for the split files')
    return parser.parse_args()


def split_by_strand(input_file, output_prefix):
    with open(input_file, 'r') as infile:
        with open(f"{output_prefix}_positive.gff", 'w') as plus_out, \
             open(f"{output_prefix}_negative.gff", 'w') as minus_out:
            for line in infile:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                strand = fields[6]
                if strand == '+':
                    plus_out.write(line)
                elif strand == '-':
                    minus_out.write(line)

def main():
    # Parse command line arguments
    args = get_args()
    
    # Split the GFF file by strand
    split_by_strand(args.input_file, args.output_prefix)

if __name__ == "__main__":
    main()