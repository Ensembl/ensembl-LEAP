process OVERLAP_READS_WITH_POLYA {
    label "process_high"
    label "error_retry"
    input:
        tuple val(id), path(reads), path(polya), val(chromosome), val(strand)
    output:
        tuple val(id), val(chromosome), val(strand), path("${reads.simpleName}_all.bed"), emit:overlapped_reads
    script:
    """
    overlap_reads_with_polya_or_cage.py ${reads} ${polya} ${reads.simpleName}
    """
}