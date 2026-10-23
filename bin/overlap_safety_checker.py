#!/usr/bin/env python 
import argparse
import os
import gzip

def get_args():
    parser = argparse.ArgumentParser(description="Check the safety of overlaps based on naturally occurring overlaps.")
    parser.add_argument("naturally_occurring_overlaps", help="Path to the file containing naturally occurring overlaps.")
    parser.add_argument("input_overlaps", help="Path to the input overlaps file.")
    parser.add_argument("safe_fe_prefix", help="Path to the output file for safe FE overlaps.")
    parser.add_argument("safe_le_prefix", help="Path to the output file for safe LE overlaps.")
    parser.add_argument("unsafe_overlaps_path", help="Path to the output file for unsafe overlaps.")
    
    return parser.parse_args()

def open_file(file_path):
    if file_path.endswith('.gz'):
        return gzip.open(file_path, 'rt')
    else:
        return open(file_path, 'r')

def main():
    # Parse command line arguments
    args = get_args()
    naturally_occurring_overlaps = args.naturally_occurring_overlaps
    input_overlaps = args.input_overlaps
    safe_fe_prefix = args.safe_fe_prefix
    safe_le_prefix = args.safe_le_prefix
    unsafe_overlaps_path = args.unsafe_overlaps_path

    # Check if the naturally occurring overlaps file is empty
    if not os.path.getsize(naturally_occurring_overlaps) == 0:
        fe_end_no = []
        le_end_no = []
        with open_file(naturally_occurring_overlaps) as no_file:
            for line in no_file:
                fields = line.strip().split("\t")
                fe_end_no.append(fields[0])
                le_end_no.append(fields[1])
    else:
        fe_end_no = []
        le_end_no = []

    # Read the input overlaps file and check for safety
    safe_overlaps = []
    unsafe_overlaps = []
    with open_file(input_overlaps) as overlaps_file:
        for line in overlaps_file:
            fields = line.strip().split("\t")
            if fields[3] == "single_exon":
                fe_end = fields[6].strip().split(",")
            else:
                fe_end = fields[9].strip().split(",")
            if fields[13] == "single_exon":
                le_end = fields[16].strip().split(",")
            else:
                le_end = fields[19].strip().split(",")
            is_safe = False
            if any(fe in le_end for fe in fe_end):
                is_safe = True
            else:
                for fe, le in zip(fe_end, le_end):
                    for idx, (fe_no, le_no) in enumerate(zip(fe_end_no, le_end_no)):
                        if fe == fe_no and le == le_no:
                            is_safe = True
                            break
                    if is_safe:
                        break
            if is_safe:
                safe_overlaps.append(line)
            else:
                unsafe_overlaps.append(line)
    if safe_overlaps:
        # Split safe overlaps up into FE and LE
        safe_fe_exon = []
        safe_le_exon = []
        for i in range(len(safe_overlaps)):
            fields = safe_overlaps[i].strip().split("\t")
            fe_fields = fields[0:10]
            le_fields = fields[10:20]
            safe_fe_exon.append("\t".join(fe_fields) + "\n")
            safe_le_exon.append("\t".join(le_fields) + "\n")

        # Write safe FE and LE overlaps to separate files
        with open(f"{safe_fe_prefix}.bed", "a") as safe_fe_file:
            for overlap in safe_fe_exon:
                safe_fe_file.write(overlap)
        with open(f"{safe_le_prefix}.bed", "a") as safe_le_file:
            for overlap in safe_le_exon:
                safe_le_file.write(overlap)
    os.rename(f"{safe_fe_prefix}.bed", f"{safe_fe_prefix}C.bed")
    os.rename(f"{safe_le_prefix}.bed", f"{safe_le_prefix}C.bed")
    if unsafe_overlaps:
        # Write unsafe overlaps to a file
        with open(unsafe_overlaps_path, "w") as unsafe_file:
            for overlap in unsafe_overlaps:
                unsafe_file.write(overlap)
    else:
        os.system(f'touch {unsafe_overlaps_path}')

if __name__ == "__main__":
    main()