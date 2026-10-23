process FIX_UNSAFE_OVERLAPS {
    label "process_high"
    label "error_retry"
    publishDir "${params.output_directory}/unsafe_overlaps/", mode: 'copy', pattern: "${id}_rejected.bed"
    publishDir "${params.output_directory}/fixed_extensions/", mode: 'copy', pattern: "*F.bed"
    input:
        tuple val(id), path(fe), path(le), path(unsafe_extensions)
    output:
        tuple val(id), path("${fe.simpleName}F.bed"), path("${le.simpleName}F.bed"), path("${id}_rejected.bed"), emit:fixed_extensions
    script:
    """
    fix_unsafe_overlaps.py ${unsafe_extensions} ${fe.simpleName} ${le.simpleName} ${id} --alpha ${params.overlap_penalty}
    """
}