process SPLIT_TRANSCRIPTS_BASED_ON_EXON_NUMBER {
    label "process_high"
    label "error_retry"
    input:
        tuple val(chromosome), path(full_length_transcripts)
    output:
        tuple val("${chromosome}_negative"), path("*SE_negative_NO.bed"), emit: single_exon_negative_no
        tuple val("${chromosome}_positive"), path("*SE_positive_NO.bed"), emit: single_exon_positive_no
        tuple val("${chromosome}_negative"), path("*negative_SE.bed"), emit: single_exon_negative
        tuple val("${chromosome}_positive"), path("*positive_SE.bed"), emit: single_exon_positive
        tuple val(chromosome), val("positive"), path("*_S.gff3"), emit: single_exon_positive_gff
        tuple val(chromosome), val("negative"), path("*_S.gff3"), emit: single_exon_negative_gff
        tuple val(chromosome), path("current_annotation_*_*FT_M.gff3"), emit: multi_exon_full_length_transcripts
    script:
    """
    extract_single_exons.py ${full_length_transcripts} ${chromosome} ${full_length_transcripts.simpleName}
    """
}