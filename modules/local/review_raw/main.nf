process PREPARE_REVIEW_RAW {
    tag "${chromosome}"
    label 'process_medium'
    errorStrategy 'terminate'
    afterScript ''
    input:
    tuple val(chromosome), path(annotation)
    path helpers, stageAs: 'helpers/*'
    output:
    path("raw_${chromosome}"), emit: partition
    script:
    def q = { v -> "'" + v.toString().replace("'", "'\"'\"'") + "'" }
    """
    python3 helpers/prepare_manual_review.py raw ${q.call(chromosome)} ${q.call(annotation)} ${q.call("raw_${chromosome}")}
    """
}
