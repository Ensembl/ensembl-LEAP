process MATCH_POSSIBLE_EXTENSIONS_WITH_CURRENT_ANNOTATION {
    label "process_high"
    label "error_retry"
    input:
        tuple val(id), val(chromosome), val(strand), path(reads), path(current_annotation)
    output:
        tuple val(id), val(chromosome), val(strand), path("${reads.simpleName}_M.bed"), emit:matched_possible_extensions
    script:
    """
    overlap_reads_and_annotation.py ${current_annotation} ${reads} ${reads.simpleName} ${params.extend_on_current_annotation}
    """
}