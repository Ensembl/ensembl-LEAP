#!/usr/bin/env python 
import argparse
import gzip

def get_args():
    parser = argparse.ArgumentParser(description="Split GFF file by chromosome.")
    parser.add_argument("gff", type=str, default="current_annotation.gff3",
                        help="Path to the GFF file to split")
    return parser.parse_args()

def open_file(file_path):
        if file_path.endswith('.gz'):
            return gzip.open(file_path, 'rt')
        else:
            return open(file_path, 'r')

def main():
    # Parse command line arguments
    args = get_args()
    gff_file = args.gff
    
    # Split the GFF file by chromosome
    with open_file(gff_file) as infile:
        files = {}
        for line in infile:
            if line.startswith("#") or not line.strip():
                continue
            chrom = line.split('\t', 1)[0]
            if chrom not in files:
                files[chrom] = open(f"{chrom}.split.gff", "w")
            files[chrom].write(line)
        for f in files.values():
            f.close()

if __name__ == "__main__":
    main()