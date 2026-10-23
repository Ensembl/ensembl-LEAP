process FILTER_PARTIAL_TRANSCRIPTS {
    label "process_high"
    label "error_retry"
    input:
        tuple val(chromosome), path(current_annotation)
    output:
        tuple val(chromosome), path("*_FT.gff"), emit:current_annotation_ft
        tuple val(chromosome), path("*_PT.gff"), emit:current_annotation_pt
        tuple val(chromosome), val("positive"), path("*NPC.gff"), emit:current_annotation_npc_positive
        tuple val(chromosome), val("negative"), path("*NPC.gff"), emit:current_annotation_npc_negative
    script:
    """
    filter_out_partials.py ${current_annotation} ${current_annotation.simpleName}
    """
}