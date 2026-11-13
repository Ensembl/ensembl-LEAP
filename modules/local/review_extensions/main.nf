process PREPARE_REVIEW_EXTENSIONS {
    tag "${chromosome}_${strand}"
    label 'process_high'
    errorStrategy 'terminate'
    afterScript ''
    input:
    tuple val(chromosome), val(strand), path(raw), path(extended), path(stats)
    path helpers, stageAs: 'helpers/*'
    output:
    path("extended_${chromosome}_${strand}"), emit: partition
    script:
    def q = { v -> "'" + v.toString().replace("'", "'\"'\"'") + "'" }
    """
    python3 helpers/prepare_manual_review.py partition ${q.call(chromosome)} ${q.call(strand)} \\
      ${q.call(raw)} ${q.call(extended)} ${q.call(stats)} ${q.call("extended_${chromosome}_${strand}")}
    """
}
