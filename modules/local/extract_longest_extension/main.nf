process EXTRACT_LONGEST_EXTENSION {
    label "process_high"
    label "error_retry"
    input:
        tuple val(id), val(chromosome), val(strand), path(matched_possible_extensions)
    output:
        tuple val(id), val(chromosome), val(strand),path("*_EXT.bed"), emit:longest_extension
    script:
    """
    extract_longest_extension.py ${matched_possible_extensions} ${matched_possible_extensions.simpleName}
    """
}