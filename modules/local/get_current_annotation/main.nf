process GET_CURRENT_ANNOTATION {
    label "process_high"
    label "error_retry"
    publishDir "${params.output_directory}/raw_annotation/", mode: 'copy', pattern: "current_annotation*.gff3"
    input:
        tuple val(reads), val(reads_index), val(chromosome)
    output:
        tuple val(chromosome), path("current_annotation*.gff3"), emit:current_annotation_by_chromosome
        tuple val("${chromosome}_positive"), path("do_not_extend_transcripts_positive.txt"), emit:do_not_extend_transcripts_positive
        tuple val("${chromosome}_negative"), path("do_not_extend_transcripts_negative.txt"), emit:do_not_extend_transcripts_negative
    script:
    """
    get_current_annotation.sh "${params.host}" "${params.port}" "${params.user}" "${params.pass}" "${params.dbname}" "${chromosome}"
    """
}