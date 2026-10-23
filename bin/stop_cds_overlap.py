#!/usr/bin/env python 
import argparse
import gzip
from itertools import chain
import pandas as pd

def get_args():
    parser = argparse.ArgumentParser(description="Overlap CDS with extensions to ensure no overlap with CDS regions.")
    parser.add_argument("gff_file", help="GFF file containing CDS regions.")
    parser.add_argument("input_extensions", help="Input file for extensions.")
    parser.add_argument("output_prefix", help="Prefix for output files.")
    parser.add_argument("single_exon_regions", help="File containing single exon transcripts.")
    parser.add_argument("non_protein_coding_regions", help="File containing non-protein coding regions.")
    parser.add_argument("natural_overlaps_file", help="File containing natural overlaps.")
    return parser.parse_args()

def open_file(file_path):
    if file_path.endswith('.gz'):
        return gzip.open(file_path, 'rt')
    else:
        return open(file_path, 'r')

def get_exon_coords(group):
    exon_rows = group[group['type'] == 'exon']
    if not exon_rows.empty:
        # If multiple exon rows exist, check if they have the same coordinates
        if len(exon_rows) > 1:
            # Check if all exon rows have the same start/end
            starts = exon_rows['start'].unique()
            ends = exon_rows['end'].unique()
            if len(starts) > 1 or len(ends) > 1:
                print(f"Warning: Multiple exon rows with different coordinates for transcript {group.iloc[0]['transcript_id']}, exon {group.iloc[0]['exon_number']}")
                print(f"Starts: {starts}, Ends: {ends}")
        
        # Use coordinates from first exon row
        return pd.Series({
            'chromosome': exon_rows.iloc[0]['chromosome'],
            'start': exon_rows.iloc[0]['start'],
            'end': exon_rows.iloc[0]['end'],
            'type': list(group['type'])
        })
    else:
        # Fallback to first row if no exon found
        return pd.Series({
            'chromosome': group.iloc[0]['chromosome'],
            'start': group.iloc[0]['start'],
            'end': group.iloc[0]['end'],
            'type': list(group['type'])
        })

def split_gff(input_gff, remove_utr_type, natural_overlaps_ids):
    # Extract all CDS features from the GFF file
    cds_features = []
    single_cds_exons = []
    with open(input_gff, 'r') as gff_file:
        for line in gff_file:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if fields[2] != "CDS":
                continue
            attr_dict = dict(item.split("=") for item in fields[-1].split(';') if '=' in item)
            if attr_dict["transcript_type"] != "protein_coding":
                continue
            cds_features.append({"chromosome": fields[0],
                                    "start": int(fields[3]),
                                    "end": int(fields[4]),
                                    "gene_id": attr_dict["gene_id"]})
    all_cds_features = pd.DataFrame(cds_features).drop_duplicates()
    # Extract gene ids for genes that only contain a single CDS but have multiple exons
    single_cds_gene_ids = all_cds_features.groupby("gene_id").size().reset_index(name='counts').query('counts == 1')['gene_id'].tolist()
    single_cds_gene_ids = [gene_id for gene_id in single_cds_gene_ids if gene_id in natural_overlaps_ids]
    # Extract all exons for genes with a single CDS but multiple exons except for the opposite end exon to the one being processed in this script
    with open(input_gff, 'r') as gff_file:
        for line in gff_file:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            attr_dict = dict(item.split("=") for item in fields[-1].split(';') if '=' in item)
            if attr_dict['gene_id'] in single_cds_gene_ids:
                if fields[2] in ["exon", "five_prime_UTR", "three_prime_UTR", "CDS"]:
                    single_cds_exons.append({"chromosome": fields[0],
                                                "type": fields[2],
                                                "start": int(fields[3]),
                                                "end": int(fields[4]),
                                                "transcript_id": attr_dict["transcript_id"],
                                                "gene_id": attr_dict["gene_id"],
                                                "exon_number": attr_dict.get("exon_number", None)
                                                })
    single_cds_exons = pd.DataFrame(single_cds_exons)
    if not single_cds_exons.empty:
        single_cds_exons_combined = single_cds_exons.groupby(['transcript_id', 'gene_id', 'exon_number']).apply(get_exon_coords, include_groups=False).reset_index()
        if remove_utr_type is not None:
            other_type = 'three_prime_UTR' if remove_utr_type == 'five_prime_UTR' else 'five_prime_UTR'
            # Only filter if remove_utr_type is present AND other_type is NOT present
            single_cds_exons_filtered = single_cds_exons_combined[
                ~(single_cds_exons_combined['type'].apply(lambda x: remove_utr_type in x) & 
                    ~single_cds_exons_combined['type'].apply(lambda x: other_type in x))
            ].reset_index(drop=True)
        else:
            single_cds_exons_filtered = single_cds_exons_combined.reset_index(drop=True)
    else:
        single_cds_exons_filtered = single_cds_exons.copy()
    return all_cds_features[['chromosome', 'start', 'end']].drop_duplicates(), single_cds_gene_ids, single_cds_exons_filtered

def filter_exons_by_side(exons_df, side, strand):
    """Extract all exons on one side of the CDS based on side and strand"""
    filtered_exons = []
    transcript_ids = exons_df['transcript_id'].unique()
    if strand == "negative":
            side = 'right' if side == 'left' else 'left'
    
    for transcript_id in transcript_ids:
        exons_by_transcript = exons_df[exons_df['transcript_id'] == transcript_id]
        cds_row = exons_by_transcript[exons_by_transcript['type'].apply(lambda x: 'CDS' in x)]
        
        if cds_row.empty:
            continue
            
        cds_exon_number = float(cds_row['exon_number'].iloc[0])
        
        # For left side: get all exons with number < cds_exon_number
        # For right side: get all exons with number > cds_exon_number
        if side == 'left':
            # Get all exons with exon number less than CDS
            mask = exons_by_transcript['exon_number'].astype(float) < cds_exon_number
            selected = exons_by_transcript[mask]
            # If no exons found on the left, include the CDS exon
            if selected.empty:
                selected = cds_row
        else:  # side == 'right'
            # Get all exons with exon number greater than CDS
            mask = exons_by_transcript['exon_number'].astype(float) > cds_exon_number
            selected = exons_by_transcript[mask]
            # If no exons found on the right, include the CDS exon
            if selected.empty:
                selected = cds_row
        
        if not selected.empty:
            filtered_exons.append(selected)
    
    if filtered_exons:
        return pd.concat(filtered_exons, ignore_index=True)
    else:
        return pd.DataFrame(columns=exons_df.columns)

def overlap_region_with_feature(npc, fields, exon_type, strand):
    # Process the fields to create all possible extensions with the longest first
    if fields[3].strip() != ".":
        fields[1], fields[2] = int(fields[1]), int(fields[2])
        if strand == "positive":
            if exon_type == "FE":
                splice_site_index = 2
                extension_index = 1
                original_region = fields.copy()
                original_region[1] = int(fields[3].strip().split("/")[-1])
            else:
                splice_site_index = 1
                extension_index = 2
                original_region = fields.copy()
                original_region[2] = int(fields[3].strip().split("/")[-1])
        else:
            if exon_type == "FE":
                splice_site_index = 1
                extension_index = 2
                original_region = fields.copy()
                original_region[2] = int(fields[3].strip().split("/")[-1])
            else:
                splice_site_index = 2
                extension_index = 1
                original_region = fields.copy()
                original_region[1] = int(fields[3].strip().split("/")[-1])
        
        possible_extensions = [fields[extension_index]] + fields[3].strip().split("/")[:-1]
        replace_pos = possible_extensions.copy()
        replace_pos = replace_pos + [fields[3].strip().split("/")[-1]]
        
        safe_extension_found = False
        safe_extended_region = None
        
        for pos in possible_extensions:
            replace_pos.remove(pos)
            if replace_pos == []:
                replace_pos = "."
            check_extension = fields.copy()
            check_extension[1] = min(int(pos), int(fields[splice_site_index]))
            check_extension[2] = max(int(pos), int(fields[splice_site_index]))
            check_extension[3] = "/".join(map(str, replace_pos))
            
            # Extract extended region
            if strand == "positive":
                if exon_type == "FE":
                    extended_region = check_extension.copy()
                    extended_region[0], extended_region[1], extended_region[2] = check_extension[0], check_extension[1], original_region[1]
                else:
                    extended_region = check_extension.copy()
                    extended_region[0], extended_region[1], extended_region[2] = check_extension[0], original_region[2], check_extension[2]
            else:
                if exon_type == "FE":
                    extended_region = check_extension.copy()
                    extended_region[0], extended_region[1], extended_region[2] = check_extension[0], original_region[2], check_extension[2]
                else:
                    extended_region = check_extension.copy()
                    extended_region[0], extended_region[1], extended_region[2] = check_extension[0], check_extension[1], original_region[1]
            
            # Check for overlap with CDS
            overlap_found = False
            not_overlapping_regions = npc[(npc['end'] < extended_region[1]) | (npc['start'] > extended_region[2])]
            if not len(not_overlapping_regions) == len(npc):
                overlap_found = True
            
            # If no overlap found, use this extension
            if not overlap_found:
                safe_extended_region = check_extension.copy()
                safe_extension_found = True
                break
        
        # If no safe extension found, use original region
        if not safe_extension_found:
            safe_extended_region = original_region.copy()
            safe_extended_region[3] = "."
    else:
        safe_extended_region = fields.copy()
        safe_extended_region[3] = "."
    
    return safe_extended_region

def main():
        # Parse command line arguments
    args = get_args()
    gff_file = args.gff_file
    input_extensions = args.input_extensions
    output_prefix = args.output_prefix
    single_exon_regions = args.single_exon_regions
    non_protein_coding_regions = args.non_protein_coding_regions
    natural_overlaps_file = args.natural_overlaps_file
    
    if "FE" in input_extensions:
        extension_category = "FE"
    elif "LE" in input_extensions:
        extension_category = "LE"
    else:
        raise ValueError("Input file name must contain 'FE' or 'LE' to determine extension type.")

    if "positive" in input_extensions:
        strand = "positive"
    elif "negative" in input_extensions:
        strand = "negative"
    else:
        raise ValueError("Input file name must contain 'positive' or 'negative' to determine strand type.")

    if (strand == "positive" and extension_category == "FE"):
        # remove_utr_type = 'three_prime_UTR'
        side = 'left'
    elif (strand == "negative" and extension_category == "FE"):
        # remove_utr_type = 'three_prime_UTR'
        side = 'right'
    elif (strand == "positive" and extension_category == "LE"):
        # remove_utr_type = 'five_prime_UTR'
        side = 'right'
    elif (strand == "negative" and extension_category == "LE"):
        # remove_utr_type = 'five_prime_UTR'
        side = 'left'

    # Read natural overlaps file and create combinations
    natural_overlaps_ids = set()
    with open(natural_overlaps_file, 'r') as no_file:
        for line in no_file:
            fields = line.strip().split('\t')
            if fields[0] != fields[1]:
                natural_overlaps_ids.add(fields[0])
                natural_overlaps_ids.add(fields[1])

    # Run the function with the specified GFF file and UTR type to remove
    cds_df, single_cds_genes, single_cds_exons = split_gff(gff_file, None, natural_overlaps_ids)

    # Load non-protein coding regions
    strand_symbol = "+" if strand == "positive" else "-"
    # pattern = '|'.join([f'transcript_type={t_type}' for t_type in [
    #     'snrna', 'snRNA', 'snorna', 'snoRNA', 'mirna', 'miRNA', 
    #     'misc_rna', 'misc_RNA', 'ribozyme', 'vault_rna', 
    #     'scarna', 'scaRNA', 'srna', 'sRNA', 'rRNA', 'rrna_pseudogene','rrna', 
    #     'rRNA_pseudogene'
    # ]])

    non_protein_coding_df = pd.read_csv(non_protein_coding_regions, sep="\t", header=None, names=["chromosome", "source", "type", "start", "end", "score", "strand", "phase", "attributes"])
    # non_protein_coding_df = non_protein_coding_df[~non_protein_coding_df['attributes'].str.contains(pattern, na=False)]
    non_protein_coding_exon_df = non_protein_coding_df[(non_protein_coding_df['type'] == "exon") & (non_protein_coding_df['strand'] == strand_symbol)][["chromosome", "start", "end"]].drop_duplicates()

    # Load single exon regions
    single_exon_regions_df = pd.read_csv(single_exon_regions, sep="\t", header=None, names=["chromosome", "source", "type", "start", "end", "score", "strand", "phase", "attributes"])
    # single_exon_regions_df = single_exon_regions_df[~single_exon_regions_df['attributes'].str.contains(pattern, na=False)]
    single_exon_regions_df = single_exon_regions_df[(single_exon_regions_df['type'] == "CDS") & (single_exon_regions_df['strand'] == strand_symbol)][["chromosome", "start", "end"]].drop_duplicates()

    # Combine features that cannot overlap extensions
    merged_df = pd.concat([cds_df, non_protein_coding_exon_df, single_exon_regions_df])
    
    with open_file(input_extensions) as fe_file, open(f'{output_prefix}C.bed', "w") as output_file:
        for line in fe_file:
            fields = line.strip().split("\t")
            # Find overlapping regions with CDS
            no_overlaps = overlap_region_with_feature(merged_df, fields, extension_category, strand)
            # Replace the updated gene IDs with the original ones
            if (fields[-1] in natural_overlaps_ids) and (fields[-1] in single_cds_genes):
                filtered_single_cds_exons = filter_exons_by_side(single_cds_exons, side, strand)
                no_overlaps = overlap_region_with_feature(filtered_single_cds_exons[['chromosome', 'start', 'end']], no_overlaps, extension_category, strand)
            no_overlaps[4] = str(len(no_overlaps[3].strip().split("/")) - 1)
            # Write the updated line to the output file
            output_file.write("\t".join(map(str, no_overlaps)) + "\n")

if __name__ == "__main__":
    main()