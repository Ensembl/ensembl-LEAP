process SPLIT_CAGE_BY_CHROMOSOME_AND_STRAND {
    label "process_high"
    label "error_retry"
    input:
        path cage
        path chromosomes
    output:
        path("CAGE_sample_sheet.csv"), emit:cage_sample_sheet
    script:
    """
    prepare_cage_data.py ${cage} ${chromosomes}
    """
}