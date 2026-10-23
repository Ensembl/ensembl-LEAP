process SPLIT_READS_BY_CHROMOSOME {
    label "process_high"
    label "error_retry"
    container "docker://quay.io/biocontainers/minimap2:2.26--he4a0461_2"
    input:
        path example
    output:
        path("example.txt"), emit:example_ch
    script:
    """
    code here
    """
}