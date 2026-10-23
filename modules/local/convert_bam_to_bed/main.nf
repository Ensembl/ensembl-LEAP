process CONVERT_BAM_TO_BED {
    label "process_high"
    label "error_retry"
    container "docker://quay.io/biocontainers/bedtools:2.24--1"
    input:
        tuple val(chromosome), path(reads), path(reads_index)
    output:
        tuple val(chromosome), path("*.bed"), emit: converted_bed
    script:
        """
        bedtools bamtobed -bed12 -i ${reads} > ${reads.simpleName}.bed
        """
}