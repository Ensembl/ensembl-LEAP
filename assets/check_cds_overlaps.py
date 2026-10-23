#!/usr/bin/env python 
import os
import pandas as pd
import numpy as np
import argparse
import gzip

def open_file(file_path):
    if file_path.endswith('.gz'):
        return gzip.open(file_path, 'rt')
    else:
        return open(file_path, 'r')

def get_cds_overlap_side(row):
    if (row['cds_start'] < row['CA_start']) and (row['cds_end'] > row["CA_start"] and row["cds_end"] <= row['CA_end']):
        return 'left'
    elif (row['cds_start'] < row['CA_end'] and row['cds_start'] >= row['CA_start']) and (row['cds_end'] > row['CA_end']):
        return 'right'
    elif (row['cds_start'] > row['CA_end']):
        return 'right_new_overlap'
    elif (row['cds_end'] < row['CA_start']):
        return 'left_new_overlap'

def overlap_direction_checker(row):
    if row['cds_overlap_location'] == "left_new_overlap":
        return "check"
    if row['cds_overlap_location'] == "right_new_overlap":
        return "check"
    if row['exon_strand'] == "+":
        if row['leap_tag'] == "LEAP_5":
            if row['cds_overlap_location'] == "left":
                return "check"
            elif row['cds_overlap_location'] == "right":
                return "safe"
        if row['leap_tag'] == "LEAP_3":
            if row['cds_overlap_location'] == "left":
                return "safe"
            elif row['cds_overlap_location'] == "right":
                return "check"
    if row['exon_strand'] == "-":
        if row['leap_tag'] == "LEAP_5":
            if row['cds_overlap_location'] == "left":
                return "safe"
            elif row['cds_overlap_location'] == "right":
                return "check"
        if row['leap_tag'] == "LEAP_3":
            if row['cds_overlap_location'] == "left":
                return "check"
            elif row['cds_overlap_location'] == "right":
                return "safe"
    return "unknown"

def get_args():
    parser = argparse.ArgumentParser(description="This script checks for new cds overlaps in the extended annotation")
    parser.add_argument("current_annotation", help="Path to the current annotation GFF3 file")
    parser.add_argument("extended_annotation", help="Path to the extended annotation GFF3 file")
    parser.add_argument("output", help="Path to save the output CSV file")
    return parser.parse_args()

def main():
    args = get_args()
    extended_annotation = args.extended_annotation
    current_annotation = args.current_annotation
    temp_dir_path = args.output
    
    if not os.path.exists(temp_dir_path):
        os.system(f'mkdir {temp_dir_path}')
    
    # First pass extract leap tags and last exon number
    transcript_info = {}
    with open_file(extended_annotation) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            attr_dict = dict(item.split("=", 1) for item in fields[-1].strip().split(";") if item and "=" in item)
            if attr_dict['transcript_id'] not in transcript_info:
                transcript_info[attr_dict['transcript_id']] = {"leap_tag":None, "min_exon":1, "max_exon": 1}
            if fields[2] == "transcript":
                transcript_info[attr_dict['transcript_id']]['leap_tag'] = attr_dict["leap_tag"].strip().split(",")
            if fields[2] == "exon":
                transcript_info[attr_dict['transcript_id']]['max_exon'] = max(transcript_info[attr_dict['transcript_id']]['max_exon'], int(attr_dict["exon_number"]))
    
    # Second pass output two files the first for all cds regions. The second for the first and last exon of each transcript. Attach leap tags to the exons.
    with open_file(extended_annotation) as f, open(os.path.join(temp_dir_path, "extended_regions_cds.bed"), "w") as cds_out, open(os.path.join(temp_dir_path, "extended_regions_exons.bed"), "w") as exons_out:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            attr_dict = dict(item.split("=", 1) for item in fields[-1].strip().split(";") if item and "=" in item)
            if fields[2] == "transcript":
                continue
            if fields[2] == "CDS":
                output_line = '\t'.join([fields[0], fields[3], fields[4], attr_dict['transcript_id'], attr_dict['exon_number'], fields[6]]) + "\n"
                cds_out.write(output_line)
            if fields[2] == "exon":
                if (int(attr_dict['exon_number']) == 1) and (fields[6] == "+"):
                    if "LEAP_5" in transcript_info[attr_dict['transcript_id']]['leap_tag']:
                        output_line = '\t'.join([fields[0], fields[3], fields[4], attr_dict['transcript_id'], attr_dict['exon_number'], fields[6], "LEAP_5"]) + "\n"
                    else:
                        output_line = '\t'.join([fields[0], fields[3], fields[4], attr_dict['transcript_id'], attr_dict['exon_number'], fields[6], "LEAP_none"]) + "\n"
                elif (int(attr_dict['exon_number']) == transcript_info[attr_dict['transcript_id']]['max_exon']) and (fields[6] == "+"):
                    if "LEAP_3" in transcript_info[attr_dict['transcript_id']]['leap_tag']:
                        output_line = '\t'.join([fields[0], fields[3], fields[4], attr_dict['transcript_id'], attr_dict['exon_number'], fields[6], "LEAP_3"]) + "\n"
                    else:
                        output_line = '\t'.join([fields[0], fields[3], fields[4], attr_dict['transcript_id'], attr_dict['exon_number'], fields[6], "LEAP_none"]) + "\n"
                elif (int(attr_dict['exon_number']) == 1) and (fields[6] == "-"):
                    if "LEAP_3" in transcript_info[attr_dict['transcript_id']]['leap_tag']:
                        output_line = '\t'.join([fields[0], fields[3], fields[4], attr_dict['transcript_id'], attr_dict['exon_number'], fields[6], "LEAP_3"]) + "\n"
                    else:
                        output_line = '\t'.join([fields[0], fields[3], fields[4], attr_dict['transcript_id'], attr_dict['exon_number'], fields[6], "LEAP_none"]) + "\n"
                elif (int(attr_dict['exon_number']) == transcript_info[attr_dict['transcript_id']]['max_exon']) and (fields[6] == "-"):
                    if "LEAP_5" in transcript_info[attr_dict['transcript_id']]['leap_tag']:
                        output_line = '\t'.join([fields[0], fields[3], fields[4], attr_dict['transcript_id'], attr_dict['exon_number'], fields[6], "LEAP_5"]) + "\n"
                    else:
                        output_line = '\t'.join([fields[0], fields[3], fields[4], attr_dict['transcript_id'], attr_dict['exon_number'], fields[6], "LEAP_none"]) + "\n"
                else:
                    continue
                exons_out.write(output_line)
    
    # Now I process the current annotation taking just the exons but I don't need to add the leap tag here
    current_annotation_data = []
    with open_file(current_annotation) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            attr_dict = dict(item.split("=", 1) for item in fields[-1].strip().split(";") if item and "=" in item)
            if fields[2] != "exon":
                continue
            if fields[2] == "exon" and attr_dict['transcript_type'] == "protein_coding":
                current_annotation_data.append({'chromosome': str(fields[0]), 'CA_start': int(fields[3]), 'CA_end': int(fields[4]), 'transcript_id': str(attr_dict['transcript_id']), 'exon_number': int(attr_dict['exon_number']), 'CA_strand': str(fields[6])})
                
    current_annotation_data_df = pd.DataFrame(current_annotation_data)
    
    # Now i need to sort the bed files from the extended annotation
    os.system(f'sort -k1,1V -k2,2n {os.path.join(temp_dir_path, "extended_regions_cds.bed")} > {os.path.join(temp_dir_path, "extended_regions_cds.sorted.bed")}')
    os.system(f'sort -k1,1V -k2,2n {os.path.join(temp_dir_path, "extended_regions_exons.bed")} > {os.path.join(temp_dir_path, "extended_regions_exons.sorted.bed")}')
    
    # Now i want to reduce the CDS regions to remove redundant regions
    # Reduce CDS regions to non-overlapping intervals using bedtools merge
    os.system(f'bedtools merge -i {os.path.join(temp_dir_path, "extended_regions_cds.sorted.bed")} -s -c 4,5,6 -o distinct,distinct,distinct > {os.path.join(temp_dir_path, "extended_regions_cds.merged.bed")}')
    
    # Now i need to overlap the CDS and the exons
    os.system(f'bedtools intersect -a {os.path.join(temp_dir_path, "extended_regions_cds.merged.bed")} -b {os.path.join(temp_dir_path, "extended_regions_exons.sorted.bed")} -s -wa -wb > {os.path.join(temp_dir_path, "extended_regions_overlaps.bed")}')
    
    # Now I need to read the file in and make it a dataframe
    dtypes_for_df = {'cds_chromosome':'str','start':'int','end':'int','cds_tid':'str','cds_number':'str','strand':'str','exon_chromosome':'str','exon_start':'int','exon_end':'int','exon_tid':'str','exon_number':'int','exon_strand':'str','leap_tag':'str'}
    extended_regions_overlaps = pd.read_csv(os.path.join(temp_dir_path, "extended_regions_overlaps.bed"), sep="\t", names=['cds_chromosome','cds_start','cds_end','cds_tid','cds_number','cds_strand','exon_chromosome','exon_start','exon_end','exon_tid','exon_number','exon_strand','leap_tag'],dtype=dtypes_for_df)
    
    # First, split the comma-separated columns into lists
    extended_regions_overlaps['cds_tid'] = extended_regions_overlaps['cds_tid'].str.split(',')
    extended_regions_overlaps['cds_number'] = extended_regions_overlaps['cds_number'].str.split(',')
    extended_regions_overlaps = extended_regions_overlaps.explode('cds_tid').explode('cds_number')
    extended_regions_overlaps['cds_number'] = extended_regions_overlaps['cds_number'].str.strip().astype(int)
    
    # Now I need to join the current annotation dataframe based on chromosome, transcript id and exon number
    merged_dataframe = pd.merge(extended_regions_overlaps, current_annotation_data_df, left_on=['exon_chromosome','exon_tid','exon_number','exon_strand'], right_on=['chromosome','transcript_id','exon_number','CA_strand'], how="left")
    
    # Remove rows that have no extension
    merged_dataframe = merged_dataframe[merged_dataframe['leap_tag'] != "LEAP_none"]
    
    # Now I need to filter out rows that the cds is completely inside the unextended region i just joined the the dataframe
    merged_dataframe = merged_dataframe[~((merged_dataframe['cds_start'] >= merged_dataframe['CA_start']) & (merged_dataframe['cds_end'] <= merged_dataframe['CA_end']))]
    
    # Now i need to remove cds that are part of the extended exon
    merged_dataframe = merged_dataframe[~((merged_dataframe['cds_tid'] == merged_dataframe['exon_tid']) & (merged_dataframe['cds_number'] == merged_dataframe['exon_number']))]
    
    # Cleaning up columns
    merged_dataframe = merged_dataframe[['chromosome', 'exon_strand', 'cds_start', 'cds_end', 'cds_tid', 'exon_start', 'exon_end', 'exon_tid', 'leap_tag', 'CA_start','CA_end']]
    
    # Adding side of CDS overlap
    merged_dataframe['cds_overlap_location'] = merged_dataframe.apply(get_cds_overlap_side, axis=1)
    
    # Now i want to remove the cds_tid column as it is not needed anymore and drop duplicates
    merged_dataframe = merged_dataframe.drop(columns=['cds_tid']).drop_duplicates()
    
    # Now i need to check if a cds_overlap_location overlaps the extending direction
    merged_dataframe['safety_check'] = merged_dataframe.apply(overlap_direction_checker, axis=1)
    
    merged_dataframe.to_csv(os.path.join(temp_dir_path, "cds_overlaps_final.csv"), index=False)
    
    merged_dataframe[merged_dataframe['cds_overlap_location'].isin(["left_new_overlap", "right_new_overlap"])].to_csv(os.path.join(temp_dir_path, "cds_overlaps_final_check.csv"), index=False)

if __name__ == "__main__":
    main()