process EXTRACT_FIRST_AND_LAST_EXONS {
    label "process_high"
    label "error_retry"
    input:
        tuple val(chromosome), path(split_annotation), val(strand), val(format)
    output:
        tuple val(chromosome), path("*FE.gff"), val(strand), val("FE"), emit:first_exon_annotation
        tuple val(chromosome), path("*LE.gff"), val(strand), val("LE"), emit:last_exon_annotation
    script:
    """
    extract_first_and_last_exon.py ${split_annotation} ${split_annotation.simpleName} -f ${format}
    """
}