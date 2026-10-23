process CHECK_EXTENSION_SAFETY {
    label "process_high"
    label "error_retry"
    publishDir "${params.output_directory}/extension_safety/", mode: 'copy'
    input:
        tuple val(id), path(extension_overlaps), path(fe_safe), path(le_safe), path(natural_overlaps)
    output:
        tuple val(id), path("${fe_safe.simpleName}C.bed"), path("${le_safe.simpleName}C.bed"), path("${id}_unsafe_extensions.bed"), emit:safety_checked_extensions
    script:
    """
    overlap_safety_checker.py ${natural_overlaps} ${extension_overlaps} ${fe_safe.simpleName} ${le_safe.simpleName} ${id}_unsafe_extensions.bed
    """
}