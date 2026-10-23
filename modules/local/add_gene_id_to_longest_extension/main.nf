process ADD_GENE_ID_TO_LONGEST_EXTENSION {
    label "process_high"
    label "error_retry"
    input:
        tuple val(id), val(chromosome), val(strand), path(uss), path(current_annotation)
    output:
        tuple val(chromosome), val(strand), path("${uss.simpleName}_GID.bed"), emit:uss_with_gene_id
    script:
    """
    add_gene_and_transcript_id_to_USS.py ${current_annotation} ${uss} ${uss.simpleName} 
    """
}