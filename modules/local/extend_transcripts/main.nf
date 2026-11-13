process EXTEND_TRANSCRIPTS {
    label "process_high"
    label "error_retry"
    publishDir "${params.output_directory}/extended_transcripts", mode: 'copy', pattern: "*_extended.gff"
    publishDir "${params.output_directory}/extended_transcripts_stats", mode: 'copy', pattern: "*_stats.tsv"
    publishDir "${params.output_directory}/fixed_extensions", mode: 'copy', pattern: "*_protected.bed"
    publishDir "${params.output_directory}/partial_protection", mode: 'copy', pattern: "*_partial_protection.tsv"
    input:
        tuple val(id), path(fe), path(le), path(rejected), path(current_annotation), path(blacklist_transcripts), path(partial_annotation), path(natural_overlaps)
    output:
        path("${id}_extended.gff"), emit:extended_transcripts
        path("${id}_stats.tsv"), emit:extended_transcripts_stats
        tuple val(id), path("${id}_FE_protected.bed"), path("${id}_LE_protected.bed"), emit:protected_extensions
        path("${id}_partial_protection.tsv"), emit:partial_protection_report
        tuple val(id), path("${id}_extended.gff"), path("${id}_stats.tsv"),
            path("${id}_FE_protected.bed"), path("${id}_LE_protected.bed"),
            path("${id}_partial_protection.tsv"), emit: datacheck_outputs
    script:
    """
    protect_partial_extensions.py ${fe} ${le} ${current_annotation} ${partial_annotation} ${natural_overlaps} ${id}
    build_chromosome_gff_and_stats.py ${id}_FE_protected.bed ${id}_LE_protected.bed ${current_annotation} ${blacklist_transcripts} ${id}_extended.gff ${id}_stats.tsv --partial-annotation ${partial_annotation}
    """
}
