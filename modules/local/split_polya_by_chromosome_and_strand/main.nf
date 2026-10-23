process SPLIT_POLYA_BY_CHROMOSOME_AND_STRAND {
    label "process_high"
    label "error_retry"
    input:
        path polya
        path chromosomes
    output:
        path("polya_sample_sheet.csv"), emit:polya_sample_sheet
    script:
    """
    prepare_polya_data.py ${polya} ${chromosomes}
    """
}