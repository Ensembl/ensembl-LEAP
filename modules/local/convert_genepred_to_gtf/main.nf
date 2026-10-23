process CONVERT_GENEPRED_TO_GTF {
    label "process_high"
    label "error_retry"
    container "docker://quay.io/biocontainers/ucsc-genepredtogtf:377--ha8a8165_5"
    input:
        tuple val(chromosome), path(genePred)
    output:
        tuple val(chromosome), path("*.gtf"), emit: converted_gtf
    script:
        """
        genePredToGtf file ${genePred} ${genePred.simpleName}.gtf
        """
}