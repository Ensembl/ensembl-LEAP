#!/usr/bin/env python 
import argparse
import gzip
import os

def get_args():
    parser = argparse.ArgumentParser(description="Match USS from reads with USS from current annotation.")
    parser.add_argument("input_current_annotation", help="Input current annotation file in CSV format.")
    parser.add_argument("input_reads_uss", help="Input reads USS file in BED format.")
    parser.add_argument("output_prefix", help="Prefix for output files.")
    parser.add_argument("allow_current_annotation", help="Allow current annotation to be included in the output.")
    return parser.parse_args()

def open_file(file_path):
    if file_path.endswith('.gz'):
        return gzip.open(file_path, 'rt')
    else:
        return open(file_path, 'r')

def main():
    # Parse command line arguments
    args = get_args()
    input_current_annotation = args.input_current_annotation
    input_reads_uss = args.input_reads_uss
    allow_current_annotation= args.allow_current_annotation
    if allow_current_annotation.lower() == 'false':
        allow_current_annotation = False
    else:
        allow_current_annotation = True
    output_prefix = args.output_prefix
    strand = ""
    site = ""

    # Check for strand and site in the input files
    if "negative" in input_reads_uss:
        strand = "negative"
    else:
        strand = "positive"

    if "_FE_" in input_current_annotation:
        site = "FE"
    else:
        site = "LE"

    # Convert current annotation to BED format
    with open_file(input_current_annotation) as current_annotation:
        header = next(current_annotation)
        for line in current_annotation:
            fields = line.strip().split(',')
            info = fields[0].split("_")
            positions = sorted([int(fields[1]), int(info[1])])
            read_row = {
                "splice_id": fields[0],
                'chromosome': info[0],
                'start': min(positions),
                'end': max(positions),
                'strand': info[2],
                'category': fields[2]
            }
            bed_formatted_line = f'{read_row["chromosome"]}\t{read_row["start"]}\t{read_row["end"]}\t{read_row["category"]}\t0\t{read_row["strand"]}\n'
            with open(f'{output_prefix}_CA.bed', "a") as out_f:
                out_f.write(bed_formatted_line)

    # Setting up splice site and other positions columns based on strand and site
    if (strand == "positive" and site == "FE") or (strand == "negative" and site == "LE"):
        splice_site_column = 2
        other_positions_column = 1
    else:
        splice_site_column = 1
        other_positions_column = 2

    # Combine current annotation with reads USS
    os.system(f'cat {output_prefix}_CA.bed {input_reads_uss} > {output_prefix}_combined.bed')

    # Match current annotation with reads USS
    exon_matcher = {}
    with open_file(f'{output_prefix}_combined.bed') as combined_bed:
        for line in combined_bed:
            fields = line.strip().split('\t')
            if fields[splice_site_column] not in exon_matcher:
                exon_matcher[fields[splice_site_column]] = [fields[3]]
            else:
                exon_matcher[fields[splice_site_column]].append(fields[3])

    # Filter splice sites that match between current annotation and reads USS
    current_annotation_has_reads = []
    for exon, reads in exon_matcher.items():
        if len(set(reads)) > 1 or ("current_annotation" in reads and allow_current_annotation):
            current_annotation_has_reads.append(exon)

    # Write the filtered current annotation to a new BED file
    with open_file(f'{output_prefix}_combined.bed') as combined_bed:
        with open(f'temp.bed', 'w') as final_bed:
            for line in combined_bed:
                fields = line.strip().split('\t')
                if fields[splice_site_column] in current_annotation_has_reads:
                    final_bed.write(line)

    # Clean up temporary files and sort the final output
    if (strand == "positive" and site == "FE") or (strand == "negative" and site == "LE"):
        os.system(f'sort -k{splice_site_column + 1},{splice_site_column + 1}n -k {other_positions_column + 1},{other_positions_column + 1}n temp.bed > {output_prefix}_M.bed')
    else:
        os.system(f'sort -k{splice_site_column + 1},{splice_site_column + 1}n -k {other_positions_column + 1},{other_positions_column + 1}nr temp.bed > {output_prefix}_M.bed')
    os.system(f'rm {output_prefix}_CA.bed {output_prefix}_combined.bed temp.bed')

if __name__ == "__main__":
    main()