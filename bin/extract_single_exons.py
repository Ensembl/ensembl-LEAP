#!/usr/bin/env python 
import argparse
import os

# input_current_annotation = "current_annotation_13.gff3"
# output_single_exon = "current_annotation_13_S"
# output_multi_exon = "current_annotation_13_M"
def get_args():
    parser = argparse.ArgumentParser(description="Extract single and multi-exon transcripts from a GFF3 file.")
    parser.add_argument("input", help="Input GFF3 file with current annotations.")
    parser.add_argument("chromosome", help="Chromosome ID")
    parser.add_argument("output_prefix", help="Output file prefix")
    return parser.parse_args()

def main():
    # Get command line arguments
    args = get_args()
    
    input_current_annotation = args.input
    output_prefix = args.output_prefix
    chromosome = args.chromosome
    output_single_exon = f'{output_prefix}_S'
    output_multi_exon = f'{output_prefix}_M'

    # Count the number of exons for each transcript
    transcript_exon_count = {}
    with open(input_current_annotation, "r") as input_file:
        for line in input_file:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if fields[2] == "exon":
                attr_dict = dict(item.split('=', 1) for item in fields[8].split(';') if '=' in item)
                if attr_dict['transcript_id'] not in transcript_exon_count:
                    transcript_exon_count[attr_dict['transcript_id']] = 1
                else:
                    transcript_exon_count[attr_dict['transcript_id']] += 1

    # Separate transcripts into single-exon and multi-exon lists
    single_exon_transcripts = []
    multi_exon_transcripts = []
    for transcript_id, info in transcript_exon_count.items():
        if info == 1:
            single_exon_transcripts.append(transcript_id)
        else:
            multi_exon_transcripts.append(transcript_id)

    # Write single-exon and multi-exon transcripts to separate files
    with open(input_current_annotation, "r") as input_file, open(f'{output_single_exon}.gff3', "a") as output_se, open(f'{output_single_exon}.bed', "a") as output_se_bed, open(f'{output_multi_exon}.gff3', "a") as output_me, open(f'{output_multi_exon}.bed',"a") as output_me_bed:
        for line in input_file:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            attr_dict = dict(item.split('=', 1) for item in fields[8].split(';') if '=' in item)
            if fields[2] != "exon":
                if attr_dict['transcript_id'] in multi_exon_transcripts:
                    output_me.write(line)
                if attr_dict['transcript_id'] in single_exon_transcripts:
                    output_se.write(line)
                continue
            if attr_dict['transcript_type'] != 'protein_coding':
                continue
            output_line = [fields[0], fields[3], fields[4], "single_exon", "0", fields[6], attr_dict['gene_id'], ".", ".", "0"]
            if attr_dict['transcript_id'] in single_exon_transcripts:
                output_line = [fields[0], fields[3], fields[4], "single_exon", "0", fields[6], attr_dict['gene_id'], ".", ".", "0"]
                output_se_bed.write("\t".join(output_line) + "\n")
                output_se.write(line)
            elif attr_dict['transcript_id'] in multi_exon_transcripts:
                output_line = [fields[0], fields[3], fields[4], "multi_exon", "0", fields[6], attr_dict['gene_id'], ".", ".", "0"]
                output_me_bed.write("\t".join(output_line) + "\n")
                output_me.write(line)

    # Overlap single-exon transcripts with multi-exon transcripts
    os.system(f'bedtools intersect -a {output_single_exon}.bed -b {output_multi_exon}.bed -s -wa -wb > single_exon_overlaps.bed')
    
    # Separate overlaps into positive and negative strands
    with open('single_exon_overlaps.bed', "r") as se_file, open(f'{chromosome}_SE_positive_NO_all.bed', "a") as se_positive_output, open(f'{chromosome}_SE_negative_NO_all.bed', "a") as se_negative_output:
        for line in se_file:
            fields = line.strip().split("\t")
            if fields[6] == fields[16]:
                output_line = "\t".join([fields[6], fields[16]]) + "\n"
                if fields[5] == "-":
                    se_negative_output.write(output_line)
                else:
                    se_positive_output.write(output_line)

    # Sort and remove duplicates from the output files
    if os.path.getsize(f'{chromosome}_SE_positive_NO_all.bed') > 0:
        os.system(f'sort {chromosome}_SE_positive_NO_all.bed | uniq > {chromosome}_SE_positive_NO.bed')
        os.system(f'rm {chromosome}_SE_positive_NO_all.bed')
    else:
        os.system(f'touch {chromosome}_SE_positive_NO.bed')

    if os.path.getsize(f'{chromosome}_SE_negative_NO_all.bed') > 0:
        os.system(f'sort {chromosome}_SE_negative_NO_all.bed | uniq > {chromosome}_SE_negative_NO.bed')
        os.system(f'rm {chromosome}_SE_negative_NO_all.bed')
    else:
        os.system(f'touch {chromosome}_SE_negative_NO.bed')
    # Split current annotation single exons by strand
    with open(f'{output_single_exon}.bed', "r") as se_file, open(f'{chromosome}_positive_SE.bed', "a") as se_positive_output, open(f'{chromosome}_negative_SE.bed', "a") as se_negative_output:
        for line in se_file:
            fields = line.strip().split("\t")
            if fields[5] == "-":
                se_negative_output.write("\t".join(fields) + "\n")
            else:
                se_positive_output.write("\t".join(fields) + "\n")
        
    os.system(f'rm single_exon_overlaps.bed {output_multi_exon}.bed {output_single_exon}.bed')

if __name__ == "__main__":
    main()
