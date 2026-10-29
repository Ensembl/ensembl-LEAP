#!/usr/bin/env bash 
# Download the current annotation from the provided database
db2gxf.pl \
-host "$1" \
-port "$2" \
-user "$3" \
-pass "$4" \
-dbname "$5" \
-chr "$6" \
-vega \
-format gff3 \
-file db_out.gff3

query="SELECT gene.stable_id
FROM transcript
JOIN transcript_attrib ON transcript.transcript_id = transcript_attrib.transcript_id
JOIN gene ON transcript.gene_id = gene.gene_id
JOIN seq_region ON transcript.seq_region_id = seq_region.seq_region_id
WHERE
    transcript_attrib.value = 'readthrough'
    AND transcript.is_current = 1
    AND LENGTH(seq_region.name) <= 2
    AND seq_region.name = '${6}';
"
mysql -h "$1" -u "$3" -p"$4" -P "$2" -D "$5" -e "$query" > readthroughs.txt

grep -v -F -f readthroughs.txt db_out.gff3 > db_out_RR.gff3

# grep -v -E "=OTTHUM|=OTTMUS" db_out_RR.gff3 > current_annotation_"$6".gff3
filter_out_ott_ids.py db_out_RR.gff3 current_annotation_"$6".gff3
query="SELECT DISTINCT transcript.stable_id
FROM transcript_attrib
JOIN transcript ON transcript.transcript_id = transcript_attrib.transcript_id
WHERE transcript.is_current = 1
    AND transcript.biotype = 'protein_coding'
    AND transcript_attrib.attrib_type_id = 54
    AND transcript_attrib.value IN ('upstream atg', 'cage supported tss', 'RP supported TIS')
    AND transcript.seq_region_id IN (
        SELECT seq_region_id FROM seq_region WHERE name = '${6}'
    );"
mysql -h "$1" -u "$3" -p"$4" -P "$2" -D "$5" -e "$query" > do_not_extend_transcripts_positive.txt
cp do_not_extend_transcripts_positive.txt do_not_extend_transcripts_negative.txt

if [ ! -s current_annotation_"$6".gff3 ]; then
    echo "Error: Output file current_annotation_${6}.gff3 is empty."
    exit 1
fi