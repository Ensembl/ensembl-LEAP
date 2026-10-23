process SPLIT_READS_BY_CHROMOSOME {
    label "process_high"
    label "error_retry"
    container "docker://quay.io/biocontainers/samtools:1.18--h50ea8bc_1"
    input:
        tuple path(reads), path(reads_index), val(chromosome)
    output:
        tuple val(chromosome), path("${reads.simpleName}_${chromosome}.bam"), path("${reads.simpleName}_${chromosome}.bam.csi"), emit:split_reads_by_chromosome
    script:
    """
    samtools view -b ${reads} ${chromosome} > temp.bam
    samtools sort temp.bam -o ${reads.simpleName}_${chromosome}.bam
    samtools index -c ${reads.simpleName}_${chromosome}.bam
    rm temp.bam
    """
}