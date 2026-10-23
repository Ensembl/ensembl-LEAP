process CONVERT_BED_TO_GENEPRED {
    label "process_high"
    label "error_retry"
    container "docker://quay.io/biocontainers/ucsc-bedtogenepred:377--ha8a8165_3"
    input:
        tuple val(chromosome), path(bed)
    output:
        tuple val(chromosome), path("*.txt"), emit: converted_genepred
    script:
        """
        bedToGenePred ${bed} ${bed.simpleName}.txt
        """
}