process FILTER_OVERLAPPING_FEATURES {
    label "process_high"
    label "error_retry"
    input:
        tuple val(id), path(fe_extensions), path(le_extensions), path(current_annotation), path(single_exon_gff), path(non_protein_coding_regions), path(natural_overlaps)
    output:
        tuple val(id), path("${fe_extensions.simpleName}C.bed"), path("${le_extensions.simpleName}C.bed"), emit:filtered_extensions
    script:
    """
    source ~/venv/bin/activate
    stop_cds_overlap.py ${current_annotation} ${fe_extensions} ${fe_extensions.simpleName} ${single_exon_gff} ${non_protein_coding_regions} ${natural_overlaps}
    stop_cds_overlap.py ${current_annotation} ${le_extensions} ${le_extensions.simpleName} ${single_exon_gff} ${non_protein_coding_regions} ${natural_overlaps}
    """
}