#!/usr/bin/env python 
import argparse
import csv
import os
import gzip

def get_args():
    parser = argparse.ArgumentParser(description="Build GFF file with extended annotations based on splice site extensions.")
    parser.add_argument("input_fe_extensions", type=str, help="Input file for five prime extensions.")
    parser.add_argument("input_le_extensions", type=str, help="Input file for three prime extensions.")
    parser.add_argument("input_annotation", type=str, help="Input GFF annotation file.")
    parser.add_argument("black_list", type=str, help="File with blacklisted transcripts to exclude from the output.")
    parser.add_argument("output_gff", type=str, help="Output GFF file with extended annotations.")
    parser.add_argument("output_stats", type=str, help="Output file for extension statistics.")
    parser.add_argument("--partial-annotation", required=True,
                        help="Classified partial GFF; its models participate in overlap processing but are omitted from final GFF and stats. May be empty.")
    return parser.parse_args()

def open_file(file_path):
    if file_path.endswith('.gz'):
        return gzip.open(file_path, 'rt')
    else:
        return open(file_path, 'r')

def position_comparison(possible_extension, current_annotation, strand, end):
    output = False
    possible_extension = int(possible_extension)
    current_annotation = int(current_annotation)
    if (strand == "positive" and end == "five" and current_annotation > possible_extension) or (strand == "positive" and end == "three" and current_annotation < possible_extension):
        output = True
    elif (strand == "negative" and end == "five" and current_annotation < possible_extension) or (strand == "negative" and end == "three" and current_annotation > possible_extension):
        output = True
    return output


def partial_transcript_keys(path):
    """Use exact IDs and locations, including version suffixes and PAR copies."""
    keys = set()
    with open_file(path) as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith('#'):
                continue
            fields = line.rstrip('\r\n').split('\t')
            if len(fields) != 9:
                raise ValueError(f'{path}:{number}: expected 9 GFF columns')
            if fields[2] == 'gene':
                continue
            attrs = dict(item.split('=', 1) for item in fields[8].split(';') if '=' in item)
            identifier = attrs.get('transcript_id')
            if not identifier or identifier == '.' or fields[0] in ('', '.') or fields[6] not in ('+', '-'):
                raise ValueError(f'{path}:{number}: partial record needs transcript_id, chromosome and strand')
            keys.add((fields[0], fields[6], identifier))
    return keys

def main():
    # Parse command line arguments
    args = get_args()
    input_fe_extensions = args.input_fe_extensions
    input_le_extensions = args.input_le_extensions
    input_annotation = args.input_annotation
    black_list = args.black_list
    output_gff = args.output_gff
    output_stats = args.output_stats
    # Exclude whole partial models before applying either end, so their child
    # records, statistics and gene-bound contributions cannot leak into output.
    partials = partial_transcript_keys(args.partial_annotation)
    
    # Determine strand and options based on input file names
    if "positive" in input_fe_extensions:
        strand = "positive"
        five_prime_options = {"splice_site_column":2, "extension_position_column":1,"annotation_splice_site_column":4,"annotation_position_column" :3}
        three_prime_options = {"splice_site_column":1, "extension_position_column":2,"annotation_splice_site_column":3,"annotation_position_column":4}
    else:
        strand = "negative"
        five_prime_options = {"splice_site_column":1, "extension_position_column":2,"annotation_splice_site_column":3,"annotation_position_column" :4}
        three_prime_options = {"splice_site_column":2, "extension_position_column":1,"annotation_splice_site_column":4,"annotation_position_column":3}

    # Read input files and build dictionaries for extensions    
    possible_five_prime_extensions = {}
    with open_file(input_fe_extensions) as f:
        for line in f:
            fields = line.strip().split('\t')
            extension_support_category = "current_annotation" if fields[3] == "." else "reads"
            reads_supporting_splice_site = fields[8]
            possible_five_prime_extensions[fields[five_prime_options['splice_site_column']]] = {"extension": fields[five_prime_options['extension_position_column']], "support": extension_support_category, "read_count": reads_supporting_splice_site}

    possible_three_prime_extensions = {}
    with open_file(input_le_extensions) as f:
        for line in f:
            fields = line.strip().split('\t')
            extension_support_category = "current_annotation" if fields[3] == "." else "reads"
            reads_supporting_splice_site = fields[8]
            possible_three_prime_extensions[fields[three_prime_options['splice_site_column']]] = {"extension": fields[three_prime_options['extension_position_column']], "support": extension_support_category, "read_count": reads_supporting_splice_site}
    
    # Read blacklisted transcripts
    blacklisted_transcripts = set()
    with open_file(black_list) as f:
        for line in f:
            transcript_id = line.strip()
            if transcript_id == "stable_id":
                continue
            if transcript_id:
                blacklisted_transcripts.add(transcript_id)

    # Process the input annotation file and apply extensions 5' first
    current_transcript = {}
    with open_file(input_annotation) as f, open("temp_output.gff", "w") as output_file:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            attr_dict = dict(item.split('=', 1) for item in fields[8].split(';') if '=' in item)
            if (fields[0], fields[6], attr_dict.get('transcript_id')) in partials:
                continue
            if fields[2] == "transcript":
                transcript_fields = fields.copy()
                current_transcript[attr_dict['transcript_id']] = transcript_fields
                current_transcript[attr_dict['transcript_id']].append([])
                current_transcript[attr_dict['transcript_id']].append([])
                current_transcript[attr_dict['transcript_id']].append([])
                current_transcript[attr_dict['transcript_id']].append([])
                current_transcript[attr_dict['transcript_id']].append([])
                current_transcript[attr_dict['transcript_id']].append(transcript_fields[five_prime_options['annotation_position_column']])
                current_transcript[attr_dict['transcript_id']].append(transcript_fields[five_prime_options['annotation_splice_site_column']])
            if attr_dict['transcript_id'].split('.')[0] in blacklisted_transcripts:
                if "LEAP_blacklisted" not in current_transcript[attr_dict['transcript_id']][-3]:
                    current_transcript[attr_dict['transcript_id']][-3].append("LEAP_blacklisted")
            elif (fields[2] == "exon") and (fields[five_prime_options['annotation_position_column']] == current_transcript[attr_dict['transcript_id']][five_prime_options['annotation_position_column']] or fields[five_prime_options['annotation_position_column']] == current_transcript[attr_dict['transcript_id']][-2]):
                if fields[five_prime_options['annotation_splice_site_column']] in possible_five_prime_extensions:
                    if position_comparison(possible_five_prime_extensions[fields[five_prime_options['annotation_splice_site_column']]]['extension'], fields[five_prime_options['annotation_position_column']], strand, "five"):
                        fields[five_prime_options['annotation_position_column']] = possible_five_prime_extensions[fields[five_prime_options['annotation_splice_site_column']]]['extension']
                        current_transcript[attr_dict['transcript_id']][five_prime_options['annotation_position_column']] = possible_five_prime_extensions[fields[five_prime_options['annotation_splice_site_column']]]['extension']
                        current_transcript[attr_dict['transcript_id']][-3].append("LEAP_5")
                        current_transcript[attr_dict['transcript_id']][-7].append(possible_five_prime_extensions[fields[five_prime_options['annotation_splice_site_column']]]["support"])
                        current_transcript[attr_dict['transcript_id']][-6].append(possible_five_prime_extensions[fields[five_prime_options['annotation_splice_site_column']]]["read_count"])
            if fields[2] == "five_prime_UTR" and (fields[five_prime_options['annotation_position_column']] == current_transcript[attr_dict['transcript_id']][five_prime_options['annotation_position_column']] or fields[five_prime_options['annotation_position_column']] == current_transcript[attr_dict['transcript_id']][-2]):
                if position_comparison(current_transcript[attr_dict['transcript_id']][five_prime_options['annotation_position_column']], fields[five_prime_options['annotation_position_column']], strand, "five"):
                    fields[five_prime_options['annotation_position_column']] = current_transcript[attr_dict['transcript_id']][five_prime_options['annotation_position_column']]
            if fields[2] != "transcript":
                writer = csv.writer(output_file, delimiter="\t")
                writer.writerow(fields)


    # Process the input annotation file and apply extensions 3' next
    with open("temp_output.gff", "r") as f, open(output_gff, "w") as output_file:
        for line in f:
            fields = line.strip().split('\t')
            attr_dict = dict(item.split('=', 1) for item in fields[8].split(';') if '=' in item)
            if fields[2] == "exon" and (fields[three_prime_options['annotation_position_column']] == current_transcript[attr_dict['transcript_id']][three_prime_options['annotation_position_column']] or fields[three_prime_options['annotation_position_column']] == current_transcript[attr_dict['transcript_id']][-1]):
                if fields[three_prime_options['annotation_splice_site_column']] in possible_three_prime_extensions:
                    if position_comparison(possible_three_prime_extensions[fields[three_prime_options['annotation_splice_site_column']]]['extension'], fields[three_prime_options['annotation_position_column']], strand, "three"):
                            fields[three_prime_options['annotation_position_column']] = possible_three_prime_extensions[fields[three_prime_options['annotation_splice_site_column']]]['extension']
                            current_transcript[attr_dict['transcript_id']][three_prime_options['annotation_position_column']] = possible_three_prime_extensions[fields[three_prime_options['annotation_splice_site_column']]]['extension']
                            current_transcript[attr_dict['transcript_id']][-3].append('LEAP_3')
                            current_transcript[attr_dict['transcript_id']][-5].append(possible_three_prime_extensions[fields[three_prime_options['annotation_splice_site_column']]]["support"])
                            current_transcript[attr_dict['transcript_id']][-4].append(possible_three_prime_extensions[fields[three_prime_options['annotation_splice_site_column']]]["read_count"])
            if fields[2] == "three_prime_UTR" and (fields[three_prime_options['annotation_position_column']] == current_transcript[attr_dict['transcript_id']][three_prime_options['annotation_position_column']] or fields[three_prime_options['annotation_position_column']] == current_transcript[attr_dict['transcript_id']][-1]):
                if position_comparison(current_transcript[attr_dict['transcript_id']][three_prime_options['annotation_position_column']], fields[three_prime_options['annotation_position_column']], strand, "three"):
                    fields[three_prime_options['annotation_position_column']] = current_transcript[attr_dict['transcript_id']][three_prime_options['annotation_position_column']]
            writer = csv.writer(output_file, delimiter="\t")
            writer.writerow(fields)

    # Write the final output GFF file and extension statistics
    extension_stats = {}
    with open(output_gff, "a") as output_file, open(output_stats, "w") as stats_file:
        stats_writer = csv.writer(stats_file, delimiter="\t")
        # Write header for extension_stats.tsv
        stats_writer.writerow([
            "gene_id",
            "transcript_id",
            "transcript_location",
            "five_prime_extension",
            "three_prime_extension",
            "five_prime_support",
            "five_prime_read_count",
            "three_prime_support",
            "three_prime_read_count"
        ])
        for transcript_id, transcript_data in current_transcript.items():
            attr_dict = dict(item.split('=', 1) for item in transcript_data[8].split(';') if '=' in item)
            extension_stats[transcript_id] = {"transcript_id": transcript_id, "five_prime_extension": [], "three_prime_extension": [], "five_prime_support": [], "five_prime_read_count": [], "three_prime_support": [], "three_prime_read_count": []}
            leap_tag = set()
            if "LEAP_blacklisted" in transcript_data[-3]:
                extension_stats[transcript_id]["five_prime_extension"] = 0
                extension_stats[transcript_id]["five_prime_support"] = "blacklisted"
                extension_stats[transcript_id]["five_prime_read_count"] = "blacklisted"
                leap_tag.add("LEAP_blacklisted")
            elif "LEAP_5" in transcript_data[-3]:
                extension_stats[transcript_id]["five_prime_extension"] = abs(int(transcript_data[-2]) - int(transcript_data[five_prime_options['annotation_position_column']]))
                extension_stats[transcript_id]["five_prime_support"] = transcript_data[-7]
                extension_stats[transcript_id]["five_prime_read_count"] = transcript_data[-6]
                leap_tag.add("LEAP_5")
            else:
                extension_stats[transcript_id]["five_prime_extension"] = 0
                extension_stats[transcript_id]["five_prime_support"] = "none"
                extension_stats[transcript_id]["five_prime_read_count"] = 0
            if "LEAP_3" in transcript_data[-3]:
                extension_stats[transcript_id]["three_prime_extension"] = abs(int(transcript_data[-1]) - int(transcript_data[three_prime_options['annotation_position_column']]))
                extension_stats[transcript_id]["three_prime_support"] = transcript_data[-5]
                extension_stats[transcript_id]["three_prime_read_count"] = transcript_data[-4]
                leap_tag.add("LEAP_3")
            else:
                extension_stats[transcript_id]["three_prime_extension"] = 0
                extension_stats[transcript_id]["three_prime_support"] = "none"
                extension_stats[transcript_id]["three_prime_read_count"] = 0
            if len(leap_tag) == 0:
                leap_tag.add("LEAP_none")
            leap_tag = ','.join(leap_tag)
            transcript_data = transcript_data[0:9]
            transcript_data[8] = f'{transcript_data[8]};leap_tag={leap_tag}'
            # Write transcript line to final_output_with_transcripts.gff
            transcript_line = "\t".join(transcript_data[0:9])
            output_file.write(transcript_line + "\n")
            # Write extension stats to extension_stats.tsv
            stats_writer.writerow([
                attr_dict["gene_id"],
                extension_stats[transcript_id]["transcript_id"],
                f"{transcript_data[0]}:{transcript_data[3]}-{transcript_data[4]}",
                extension_stats[transcript_id]["five_prime_extension"],
                extension_stats[transcript_id]["three_prime_extension"],
                ",".join(extension_stats[transcript_id]["five_prime_support"]) if isinstance(extension_stats[transcript_id]["five_prime_support"], list) else extension_stats[transcript_id]["five_prime_support"],
                ",".join(extension_stats[transcript_id]["five_prime_read_count"]) if isinstance(extension_stats[transcript_id]["five_prime_read_count"], list) else extension_stats[transcript_id]["five_prime_read_count"],
                ",".join(extension_stats[transcript_id]["three_prime_support"]) if isinstance(extension_stats[transcript_id]["three_prime_support"], list) else extension_stats[transcript_id]["three_prime_support"],
                ",".join(extension_stats[transcript_id]["three_prime_read_count"]) if isinstance(extension_stats[transcript_id]["three_prime_read_count"], list) else extension_stats[transcript_id]["three_prime_read_count"]
            ])
    os.remove("temp_output.gff")
    # Write gene rows to output
    genes = {}
    with open(output_gff, "a") as output_file:
        for transcript_id, transcript_data in current_transcript.items():
            transcript_data = transcript_data[0:9]
            attr_dict = dict(item.split('=', 1) for item in transcript_data[8].split(';') if '=' in item)
            gene_line = f"ID={attr_dict['gene_id']};gene_id={attr_dict['gene_id']};gene_type={attr_dict['gene_type']};gene_name={attr_dict['gene_name']}"
            if attr_dict["gene_id"] not in genes:
                gene_row = transcript_data.copy()
                gene_row[2] = "gene"
                gene_row[8] = gene_line
                genes[attr_dict["gene_id"]] = gene_row
            else:
                # Compare as integers
                if int(transcript_data[3]) < int(genes[attr_dict["gene_id"]][3]):
                    genes[attr_dict["gene_id"]][3] = transcript_data[3]
                if int(transcript_data[4]) > int(genes[attr_dict["gene_id"]][4]):
                    genes[attr_dict["gene_id"]][4] = transcript_data[4]
        for gene, gene_info in genes.items():
            output_file.write("\t".join(gene_info) + "\n")

if __name__ == "__main__":
    main()
