process EXTEND_TRANSCRIPTS {
    label "process_high"
    label "error_retry"
    publishDir "${params.output_directory}/extended_transcripts", mode: 'copy', pattern: "*_extended.gff"
    publishDir "${params.output_directory}/extended_transcripts_stats", mode: 'copy', pattern: "*_stats.tsv"
    input:
        tuple val(id), path(fe), path(le), path(rejected), path(current_annotation), path(blacklist_transcripts), path(partial_annotation)
    output:
        path("${id}_extended.gff"), emit:extended_transcripts
        path("${id}_stats.tsv"), emit:extended_transcripts_stats
    script:
    """
    build_chromosome_gff_and_stats.py ${fe} ${le} ${current_annotation} ${blacklist_transcripts} ${id}_extended.gff ${id}_stats.tsv --partial-annotation ${partial_annotation}
    """
}
