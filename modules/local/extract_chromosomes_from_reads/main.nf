process EXTRACT_CHROMOSOMES_FROM_READS {
    label "process_high"
    label "error_retry"
    container "docker://quay.io/biocontainers/samtools:1.18--h50ea8bc_1"
    input:
        path reads 
        path reads_index
    output:
        path("chromosome_sample_sheet.csv"), emit:reads_chr_sample_sheet
    script:
    """
    source ~/venv/bin/activate
    extract_chromosomes_from_reads.py ${reads} ${reads_index}
    """
}