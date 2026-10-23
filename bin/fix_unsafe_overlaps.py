#!/usr/bin/env python 
import argparse
import gzip
from itertools import product
import shutil
import os

def get_args():
    parser = argparse.ArgumentParser(description="Fix unsafe overlaps in BED files.")
    parser.add_argument("input_unsafe_extensions", help="Input BED file with unsafe extensions.")
    parser.add_argument("safe_fe_extensions", help="Output BED prefix for safe first exon extensions.")
    parser.add_argument("safe_le_extensions", help="Output BED prefix for safe last exon extensions.")
    parser.add_argument("unsafe_overlaps", help="Output BED prefix for unsafe overlaps.")
    parser.add_argument("--alpha", type=float, default=1.0, help="Weighting factor for the score calculation.")
    return parser.parse_args()

def open_file(file_path):
    if file_path.endswith('.gz'):
        return gzip.open(file_path, 'rt')
    else:
        return open(file_path, 'r')

def find_best_extension(X, Y, end_A, start_B, alpha=1.0):
    max_X = max(X)
    max_Y = max(Y)
    best_score = -float("inf")
    best_pair = (0, 0)
    best_sum = -float("inf")
    max_extension = max(0, abs(start_B - end_A) - 1)  # Ensure at least 1bp gap
    for x, y in product(X, Y):
        if x + y <= max_extension:  # Ensure regions do not touch or overlap
            x_norm = x / max_X if max_X > 0 else 0
            y_norm = y / max_Y if max_Y > 0 else 0
            score = (x_norm + y_norm) - alpha * abs(x_norm - y_norm)
            total_extension = x + y
            if score > best_score or (score == best_score and total_extension > best_sum):
                best_score = score
                best_pair = (x, y)
                best_sum = total_extension
    return best_pair, best_score

def main():
    # Parse command line arguments
    args = get_args()
    input_unsafe_extensions = args.input_unsafe_extensions
    safe_fe_extensions = args.safe_fe_extensions
    safe_le_extensions = args.safe_le_extensions
    unsafe_overlaps = args.unsafe_overlaps
    alpha = args.alpha
    
    # Creating files for output
    shutil.copy(f'{safe_fe_extensions}.bed', f'{safe_fe_extensions}F.bed')
    shutil.copy(f'{safe_le_extensions}.bed', f'{safe_le_extensions}F.bed')
    os.system(f'touch {unsafe_overlaps}_rejected.bed')
    os.system(f'touch potentially_safe_fe.bed')
    os.system(f'touch potentially_safe_le.bed')

    # Determine strand based on the input file name
    if "negative" in input_unsafe_extensions:
        strand = "negative"
    else:
        strand = "positive"

    # Open the input file and process each line
    if os.path.getsize(input_unsafe_extensions) != 0:
        with open_file(input_unsafe_extensions) as f:
            for line in f:
                fields = line.strip().split('\t')
                if strand == "positive":
                    fe_splice_site = int(fields[2])
                    # Treat "." and "single_exon" the same
                    fe_exts = fields[3].split("/") if fields[3] not in [".", "single_exon"] else []
                    fe_exts.append(fields[1])
                    fe_extension_options = [abs(int(x) - fe_splice_site) for x in fe_exts]
                    le_splice_site = int(fields[11])
                    le_exts = fields[13].split("/") if fields[13] not in [".", "single_exon"] else []
                    le_exts.append(fields[12])
                    le_extension_options = [abs(int(x) - le_splice_site) for x in le_exts]
                else:
                    fe_splice_site = int(fields[1])
                    fe_exts = fields[3].split("/") if fields[3] not in [".", "single_exon"] else []
                    fe_exts.append(fields[2])
                    fe_extension_options = [abs(int(x) - fe_splice_site) for x in fe_exts]
                    le_splice_site = int(fields[12])
                    le_exts = fields[13].split("/") if fields[13] not in [".", "single_exon"] else []
                    le_exts.append(fields[11])
                    le_extension_options = [abs(int(x) - le_splice_site) for x in le_exts]
                best_pair, score = find_best_extension(fe_extension_options, le_extension_options, fe_splice_site, le_splice_site, alpha)
                if score != float('-inf'):
                    new_fe_extension = fe_exts[fe_extension_options.index(best_pair[0])]
                    new_le_extension = le_exts[le_extension_options.index(best_pair[1])]
                    if strand == "positive":
                        fe_safe = [fields[0],new_fe_extension, fields[2], fields[3], fields[4], fields[5], fields[6], fields[7], fields[8], fields[9]]
                        le_safe = [fields[10], fields[11], new_le_extension, fields[13], fields[14], fields[15], fields[16], fields[17], fields[18], fields[19]]
                    else:
                        fe_safe = [fields[0], fields[1], new_fe_extension, fields[3], fields[4], fields[5], fields[6], fields[7], fields[8], fields[9]]
                        le_safe = [fields[10], new_le_extension, fields[12], fields[13], fields[14], fields[15], fields[16], fields[17], fields[18], fields[19]]
                    # with open(f'{safe_fe_extensions}F.bed', "a") as out_fe, open(f'{safe_le_extensions}F.bed', "a") as out_le:
                    with open(f'potentially_safe_fe.bed', "a") as out_fe, open(f'potentially_safe_le.bed', "a") as out_le:
                        out_fe.write('\t'.join(fe_safe) + '\n')
                        out_le.write('\t'.join(le_safe) + '\n')
                else:
                    with open(f"{unsafe_overlaps}_rejected.bed", "a") as out_unsafe:
                        out_unsafe.write(line.strip() + '\n')
        if strand == "positive":
            splice_site_column_fe, other_position_fe = 2, 1
            splice_site_column_le, other_position_le = 1, 2
        else:
            splice_site_column_fe, other_position_fe = 1, 2
            splice_site_column_le, other_position_le = 2, 1
        with open("potentially_safe_fe.bed") as fe_bed, open("potentially_safe_le.bed") as le_bed, open(f"{safe_fe_extensions}F.bed", "a") as out_fe, open(f"{safe_le_extensions}F.bed", "a") as out_le:
            safe_regions = {}
            for line in fe_bed:
                fields = line.strip().split('\t')
                if fields[3] == "single_exon":
                    continue
                if fields[splice_site_column_fe] not in safe_regions:
                    safe_regions[fields[splice_site_column_fe]] = fields
                else:
                    existing_region = safe_regions[fields[splice_site_column_fe]]
                    if strand == "positive":
                        if int(fields[other_position_fe]) > int(existing_region[other_position_fe]):
                            existing_region[other_position_fe] = fields[other_position_fe]
                    else:
                        if int(fields[other_position_fe]) < int(existing_region[other_position_fe]):
                            existing_region[other_position_fe] = fields[other_position_fe]
                    safe_regions[fields[splice_site_column_fe]] = existing_region
            for region in safe_regions.values():
                out_fe.write('\t'.join(region) + '\n')
            safe_regions.clear()
            for line in le_bed:
                fields = line.strip().split('\t')
                if fields[3] == "single_exon":
                    continue
                if fields[splice_site_column_le] not in safe_regions:
                    safe_regions[fields[splice_site_column_le]] = fields
                else:
                    existing_region = safe_regions[fields[splice_site_column_le]]
                    if strand == "positive":
                        if int(fields[other_position_le]) < int(existing_region[other_position_le]):
                            existing_region[other_position_le] = fields[other_position_le]
                    else:
                        if int(fields[other_position_le]) > int(existing_region[other_position_le]):
                            existing_region[other_position_le] = fields[other_position_le]
                    safe_regions[fields[splice_site_column_le]] = existing_region
            for region in safe_regions.values():
                out_le.write('\t'.join(region) + '\n')
        # os.remove("potentially_safe_fe.bed")
        # os.remove("potentially_safe_le.bed")    
                

if __name__ == "__main__":
    main()
