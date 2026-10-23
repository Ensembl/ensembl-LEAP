#!/usr/bin/env python 
import gzip
import os
import argparse


def get_args():
    parser = argparse.ArgumentParser(description="Filter out internal transcripts from LEAP output GFF files.")
    parser.add_argument("leap_output_directory", type=str, help="Directory containing LEAP output files")
    parser.add_argument("output_prefix", type=str, help="Prefix for output files")
    parser.add_argument("version", type=str, help="Version of the LEAP output")
    return parser.parse_args()

def open_file(file_path):
    if file_path.endswith('.gz'):
        return gzip.open(file_path, 'rt')
    else:
        return open(file_path, 'r')

def main():
    args = get_args()
    leap_output_directory = args.leap_output_directory
    output_prefix = args.output_prefix
    version = args.version
    # Combine extended transcript GFFs
    input_files = os.listdir(os.path.join(leap_output_directory, 'extended_transcripts'))
    output_file = os.path.join(leap_output_directory, f"{output_prefix}_{version}")
    current_iteration = 0
    for file in input_files:
        file = os.path.join(leap_output_directory, 'extended_transcripts', file)
        with open(file, "r") as f, open(f"{output_file}.unstranded.combined.gff", "a") as out_file:
            if current_iteration == 0:
                header = "##gff-version 3\n"
                current_iteration = 1
                out_file.write(header)
            for line in f:
                fields = line.strip().split("\t")
                if fields[2] == "gene":
                    continue
                attr_dict = dict(item.split('=', 1) for item in fields[8].split(';') if '=' in item)
                fields[8] = ';'.join([f"{k}={v}" for k, v in attr_dict.items()])
                out_file.write("\t".join(fields) + "\n")

    # Extract the most extreme CDS positions
    all_cds = {}
    with open_file(f"{output_file}.unstranded.combined.gff") as f:
        for line in f:
            if line.startswith("##"):
                continue
            fields = line.strip().split("\t")
            if fields[2] != 'CDS':
                continue
            attr_dict = dict(item.split("=") for item in fields[-1].split(";"))
            start_pos = int(fields[3])
            end_pos = int(fields[4])
            if f'{fields[0]}_{attr_dict["gene_id"]}' not in all_cds:
                all_cds[f'{fields[0]}_{attr_dict["gene_id"]}'] = {
                    'chromosome': fields[0], 
                    'strand': fields[6],
                    'earliest_start_cds': {'start': start_pos, 'end': end_pos},
                    'latest_end_cds': {'start': start_pos, 'end': end_pos}
                }
            else:
                # Check if this CDS has an earlier start
                if start_pos < all_cds[f'{fields[0]}_{attr_dict["gene_id"]}']['earliest_start_cds']['start']:
                    all_cds[f'{fields[0]}_{attr_dict["gene_id"]}']['earliest_start_cds'] = {'start': start_pos, 'end': end_pos}
                # Check if this CDS has a later end
                if end_pos > all_cds[f'{fields[0]}_{attr_dict["gene_id"]}']['latest_end_cds']['end']:
                    all_cds[f'{fields[0]}_{attr_dict["gene_id"]}']['latest_end_cds'] = {'start': start_pos, 'end': end_pos}

    # Match extreme cds positions with other transcripts
    transcript_matches_extreme_cds = {}
    with open_file(f"{output_file}.unstranded.combined.gff") as f:
        for line in f:
            if line.startswith("##"):
                continue
            fields = line.strip().split("\t")
            if fields[2] != 'CDS':
                continue
            attr_dict = dict(item.split("=") for item in fields[-1].split(";"))
            cds_extremes = all_cds[f"{fields[0]}_{attr_dict['gene_id']}"]
            if f"{fields[0]}_{attr_dict['transcript_id']}" not in transcript_matches_extreme_cds:
                transcript_matches_extreme_cds[f"{fields[0]}_{attr_dict['transcript_id']}"] = {'matches_start': False, 'matches_end': False}
            if int(fields[3]) == cds_extremes['earliest_start_cds']['start'] or int(fields[4]) == cds_extremes['earliest_start_cds']['end']:
                transcript_matches_extreme_cds[f"{fields[0]}_{attr_dict['transcript_id']}"]['matches_start'] = True
            if int(fields[3]) == cds_extremes['latest_end_cds']['start'] or int(fields[4]) == cds_extremes['latest_end_cds']['end']:
                transcript_matches_extreme_cds[f"{fields[0]}_{attr_dict['transcript_id']}"]['matches_end'] = True

    # Extract extreme exon numbers so they can be removed and replaced
    extreme_exon_numbers = {}
    with open_file(f"{output_file}.unstranded.combined.gff") as f, open(f"{output_file}.external.gff", "w") as external_out, open(f"{output_file}.potential_internal.gff", "w") as internal_out:
        for line in f:
            if line.startswith("##"):
                continue
            fields = line.strip().split("\t")
            attr_dict = dict(item.split("=") for item in fields[-1].split(";"))
            if transcript_matches_extreme_cds[f"{fields[0]}_{attr_dict['transcript_id']}"]['matches_start'] and transcript_matches_extreme_cds[f"{fields[0]}_{attr_dict['transcript_id']}"]['matches_end']:
                attr_dict["transcript_class"] = "external"
                fields[-1] = ';'.join([f"{k}={v}" for k, v in attr_dict.items()])
                external_out.write("\t".join(fields) + "\n")
                # external_out.write(line)
            else:
                internal_out.write(line)
                # Track exon numbers for potential internal transcripts
                if fields[2] == "exon":
                    transcript_id = f"{fields[0]}_{attr_dict['transcript_id']}"
                    strand = fields[6]
                    exon_number = int(attr_dict['exon_number'])

                    if transcript_id not in extreme_exon_numbers:
                        extreme_exon_numbers[transcript_id] = {
                            'first_exon_number': exon_number,
                            'last_exon_number': exon_number
                        }
                    else:
                        # Update first and last exon numbers based on strand
                        if strand == '+':
                            # Positive strand: first = smallest number, last = largest number
                            if exon_number < extreme_exon_numbers[transcript_id]['first_exon_number']:
                                extreme_exon_numbers[transcript_id]['first_exon_number'] = exon_number
                            if exon_number > extreme_exon_numbers[transcript_id]['last_exon_number']:
                                extreme_exon_numbers[transcript_id]['last_exon_number'] = exon_number
                        else:
                            # Negative strand: first = largest number, last = smallest number
                            if exon_number > extreme_exon_numbers[transcript_id]['first_exon_number']:
                                extreme_exon_numbers[transcript_id]['first_exon_number'] = exon_number
                            if exon_number < extreme_exon_numbers[transcript_id]['last_exon_number']:
                                extreme_exon_numbers[transcript_id]['last_exon_number'] = exon_number

    # Removing terminal exons as needed
    with open_file(f"{output_file}.potential_internal.gff") as internal_file, open(f"{output_file}.removed_exons.gff", "w") as internal_out:
        for line in internal_file:
            if line.startswith("##"):
                continue
            fields = line.strip().split("\t")
            attr_dict = dict(item.split("=") for item in fields[-1].split(";"))
            exon_to_remove = []
            if fields[2] == "transcript":
                internal_out.write(line)
                continue
            if not transcript_matches_extreme_cds[f"{fields[0]}_{attr_dict['transcript_id']}"]['matches_start']:
                exon_to_remove.append(str(extreme_exon_numbers[f"{fields[0]}_{attr_dict['transcript_id']}"]['first_exon_number']))
            if not transcript_matches_extreme_cds[f"{fields[0]}_{attr_dict['transcript_id']}"]['matches_end']:
                exon_to_remove.append(str(extreme_exon_numbers[f"{fields[0]}_{attr_dict['transcript_id']}"]['last_exon_number']))
            if attr_dict['exon_number'] in exon_to_remove:
                continue
            internal_out.write(line)

    # Loading unexteded annotation and extracting exons to replace the internal extensions
    altered_transcript_start_end = {}
    with open_file(f'{leap_output_directory}/db_annotation.combined.gff3.gz') as original_file, open(f"{output_file}.exons_to_add.gff","w") as internal_file:
        for line in original_file:
            if line.startswith("##"):
                continue
            fields = line.strip().split("\t")
            attr_dict = dict(item.split("=") for item in fields[-1].split(";"))
            if fields[2] == "transcript":
                continue
            if f"{fields[0]}_{attr_dict['transcript_id']}" not in transcript_matches_extreme_cds:
                continue
            exons_to_add = []
            if f"{fields[0]}_{attr_dict['transcript_id']}" not in altered_transcript_start_end:
                altered_transcript_start_end[f"{fields[0]}_{attr_dict['transcript_id']}"] = {"new_start": None, "new_end": None}
            if not transcript_matches_extreme_cds[f"{fields[0]}_{attr_dict['transcript_id']}"]['matches_start']:
                exons_to_add.append(str(extreme_exon_numbers[f"{fields[0]}_{attr_dict['transcript_id']}"]['first_exon_number']))
                if attr_dict['exon_number'] == str(extreme_exon_numbers[f"{fields[0]}_{attr_dict['transcript_id']}"]['first_exon_number']) and fields[2] == "exon":
                    altered_transcript_start_end[f"{fields[0]}_{attr_dict['transcript_id']}"]["new_start"] = fields[3]
            if not transcript_matches_extreme_cds[f"{fields[0]}_{attr_dict['transcript_id']}"]['matches_end']:
                exons_to_add.append(str(extreme_exon_numbers[f"{fields[0]}_{attr_dict['transcript_id']}"]['last_exon_number']))
                if attr_dict['exon_number'] == str(extreme_exon_numbers[f"{fields[0]}_{attr_dict['transcript_id']}"]['last_exon_number']) and fields[2] == "exon":
                    altered_transcript_start_end[f"{fields[0]}_{attr_dict['transcript_id']}"]["new_end"] = fields[4]
            if attr_dict['exon_number'] in exons_to_add:
                internal_file.write(line)

    # Combine the exons to add with the gff containing the remaining exons
    os.system(f'cat {output_file}.exons_to_add.gff {output_file}.removed_exons.gff > {output_file}.combined.gff')

    # Correcting transcript start and end positions in the combined GFF and altering leap tags as needed
    with open_file(f'{output_file}.combined.gff') as input_gff, open(f'{output_file}.combined.fixed.gff', "w") as output_gff:
        for line in input_gff:
            if line.startswith("##"):
                continue
            fields = line.strip().split("\t")
            attr_dict = dict(item.split("=") for item in fields[-1].split(";"))
            if fields[2] != "transcript":
                output_gff.write(line)
                continue
            leap_tags = attr_dict['leap_tag'].split(',')
            transcript_class = []
            if f"{fields[0]}_{attr_dict['transcript_id']}" in altered_transcript_start_end:
                if altered_transcript_start_end[f"{fields[0]}_{attr_dict['transcript_id']}"]['new_start'] is not None:
                    fields[3] = altered_transcript_start_end[f"{fields[0]}_{attr_dict['transcript_id']}"]['new_start']
                    transcript_class.append('internal_left')
                    if fields[6] == "+":
                        leap_tags = [item for item in leap_tags if item != "LEAP_5"]
                    else:
                        leap_tags = [item for item in leap_tags if item != "LEAP_3"]
                if altered_transcript_start_end[f"{fields[0]}_{attr_dict['transcript_id']}"]['new_end'] is not None:
                    fields[4] = altered_transcript_start_end[f"{fields[0]}_{attr_dict['transcript_id']}"]['new_end']
                    transcript_class.append('internal_right')
                    if fields[6] == "+":
                        leap_tags = [item for item in leap_tags if item != "LEAP_3"]
                    else:
                        leap_tags = [item for item in leap_tags if item != "LEAP_5"]
                if leap_tags == []:
                    leap_tags = ['LEAP_none']
                attr_dict['leap_tag'] = ','.join(leap_tags)
                attr_dict['transcript_class'] = ','.join(transcript_class)
                fields[-1] = ';'.join([f"{k}={v}" for k, v in attr_dict.items()])
                output_gff.write(('\t').join(fields) + "\n")

    # Combining the final output before sorting it
    os.system(f'(echo "##gff-version 3"; cat {output_file}.combined.fixed.gff {output_file}.external.gff) > {output_file}.final.gff')
    os.system(f"perl /hps/software/users/ensembl/havana/rmerritt/scripts/gff3sort-1.0.0/gff3sort.pl --precise {output_file}.final.gff > {output_file}.final.sorted.gff")


    # Making the final output viewable by removing references to genes in the gff
    with open(f"{output_file}.final.gff", "r") as f, open(f"{output_file}.finalRG.gff", "w") as out_file:
        for line in f:
            if line.startswith("##"):
                out_file.write(line)
                continue
            fields = line.strip().split("\t")
            if fields[2] == "gene":
                continue
            attr_dict = dict(item.split('=', 1) for item in fields[8].split(';') if '=' in item)
            keys_to_remove = ["gene_id", "gene_type", "gene_name"]
            if fields[2] == "transcript":
                keys_to_remove = ["gene_id", "gene_type", "gene_name", "Parent"]
            attr_dict = {k: v for k, v in attr_dict.items() if k not in keys_to_remove}
            fields[8] = ';'.join([f"{k}={v}" for k, v in attr_dict.items()])
            out_file.write("\t".join(fields) + "\n")

    # Sort the output with removed genes before bgziping and tabixing
    os.system(f"perl /hps/software/users/ensembl/havana/rmerritt/scripts/gff3sort-1.0.0/gff3sort.pl --precise {output_file}.finalRG.gff > {output_file}.finalRG.sorted.gff")
    os.system(f'bgzip {output_file}.finalRG.sorted.gff')
    os.system(f'tabix -p gff {output_file}.finalRG.sorted.gff.gz')

    # Some text files containing which transcripts have which class will do more work here as more stats are needed
    with open(f"{output_file}_internal_both.txt", "w") as leap_internal_both, \
        open(f"{output_file}_internal_left.txt", "w") as leap_internal_left, \
        open(f"{output_file}_internal_right.txt", "w") as leap_internal_right, \
        open(f"{output_file}_external.txt", "w") as leap_external:
        for key, value in transcript_matches_extreme_cds.items():
            key = key.split('_')[1].split('.')[0]  # Extract transcript ID
            if value['matches_start'] and value['matches_end']:
                leap_external.write(f"{key}\n")
            elif not value['matches_start'] and value['matches_end']:
                leap_internal_left.write(f"{key}\n")
            elif value['matches_start'] and not value['matches_end']:
                leap_internal_right.write(f"{key}\n")
            else:
                leap_internal_both.write(f"{key}\n")

    # Remove temporary files
    os.system(f'rm {output_file}.combined.fixed.gff {output_file}.combined.gff {output_file}.exons_to_add.gff {output_file}.external.gff {output_file}.final.gff {output_file}.potential_internal.gff {output_file}.removed_exons.gff {output_file}.unstranded.combined.gff')

if __name__ == "__main__":
    main()