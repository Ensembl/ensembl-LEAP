process EXTRACT_UNIQUE_SPLICE_SITES_FROM_READS {
    label "process_high"
    label "error_retry"
    input:
        tuple val(chromosome), path(gff), val(strand), val(exon_type)
    output:
        tuple val(chromosome), path("${gff.simpleName}_USS.csv"), val(strand), val(exon_type), emit:unique_splice_sites
    script:
    """
    extract_unique_splice_sites_from_reads.py ${gff} ${gff.simpleName} ${params.minimum_read_support}
    """
}