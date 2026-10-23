process EXTRACT_UNIQUE_SPLICE_SITES_FROM_CURRENT_ANNOTATION {
    label "process_high"
    label "error_retry"
    input:
        tuple val(chromosome), path(gff_file), val(strand), val(exon_type)
    output:
        tuple val(chromosome), path("${gff_file.simpleName}_USS.csv"), val(strand), val(exon_type), emit:unique_splice_sites
    script:
    """
    extract_unique_splice_sites_from_current_annotation.py ${gff_file} ${gff_file.simpleName}
    """
}