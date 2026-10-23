process FIND_NATURAL_OVERLAPS_IN_CURRENT_ANNOTATION {
    label "process_high"
    label "error_retry"
    container "docker://quay.io/biocontainers/bedtools:2.27.1--h077b44d_9"
    publishDir "${params.output_directory}/natural_overlaps/", mode: 'copy', pattern: "*_NO.txt"
    input:
        tuple val(id), path(fe_gff), path(le_gff), path(single_exon_no)
    output:
        tuple val(id), path("${id}_NO.txt"), emit:natural_overlaps
    script:
    """
    source ~/venv/bin/activate
    find_natural_overlaps.py  ${fe_gff} ${le_gff} temp
    cat temp_NO.txt ${single_exon_no} > ${id}_NO.txt
    """
}