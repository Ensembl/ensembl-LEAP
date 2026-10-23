process SPLIT_CURRENT_ANNOTATION_BY_STRAND {
    label "process_high"
    label "error_retry"
    input:
        tuple val(chromosome), path(current_annotation)
        val format
    output:
        tuple val(chromosome), path("*positive.gff"), val("positive"), val(format), emit:positive_annotation
        tuple val(chromosome), path("*negative.gff"), val("negative"), val(format), emit:negative_annotation
    script:
    """
    split_current_annotation_by_strand.py ${current_annotation} ${current_annotation.simpleName}
    """
}