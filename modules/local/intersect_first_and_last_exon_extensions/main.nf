process INTERSECT_FIRST_AND_LAST_EXON_EXTENSIONS {
    label "process_high"
    label "error_retry"
    container "docker://quay.io/biocontainers/bedtools:2.27.1--h077b44d_9"
    input:
        tuple val(id), path(fe), path(le), path(single_exon_bed)
    output:
        tuple val(id), path("extension_overlaps.bed"), path("${fe.simpleName}_S.bed"), path("${le.simpleName}_S.bed"), emit:overlapped_extensions
    script:
    """
    bedtools intersect -a ${fe} -b ${single_exon_bed} -wa -wb > FE_SEO.bed
    bedtools intersect -a ${single_exon_bed} -b ${le} -wa -wb > LE_SEO.bed
    bedtools intersect -a ${fe} -b ${le} -wa -wb > FE_LE_O.bed
    cat FE_SEO.bed LE_SEO.bed FE_LE_O.bed > extension_overlaps.bed
    bedtools intersect -a ${fe} -b ${le} -v > ${fe.simpleName}_T.bed
    bedtools intersect -a ${fe} -b ${single_exon_bed} -v > ${fe.simpleName}_S.bed
    bedtools intersect -a ${le} -b ${fe} -v > ${le.simpleName}_T.bed
    bedtools intersect -a ${le} -b ${single_exon_bed} -v > ${le.simpleName}_S.bed
    """
}