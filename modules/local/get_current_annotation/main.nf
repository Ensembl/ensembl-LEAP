process GET_CURRENT_ANNOTATION {
    label "process_high"
    label "error_retry"
    conda "${moduleDir}/environment.yml"
    publishDir "${params.output_directory}/raw_annotation/", mode: 'copy', pattern: "current_annotation*.gff3"
    input:
        tuple val(reads), val(reads_index), val(chromosome)
        path db_secrets, stageAs: 'connection/db_secrets.json'
        path export_helpers, stageAs: 'helpers/*'
    output:
        tuple val(chromosome), path("current_annotation*.gff3"), emit:current_annotation_by_chromosome
        tuple val("${chromosome}_positive"), path("do_not_extend_transcripts_positive.txt"), emit:do_not_extend_transcripts_positive
        tuple val("${chromosome}_negative"), path("do_not_extend_transcripts_negative.txt"), emit:do_not_extend_transcripts_negative
    script:
    def quote = { value -> "'" + value.toString().replace("'", "'\"'\"'") + "'" }
    """
    PYTHONDONTWRITEBYTECODE=1 python3 helpers/get_current_annotation.py ${quote.call(db_secrets)} ${quote.call(chromosome)}
    """
}
