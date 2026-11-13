process BUILD_REVIEW_BUNDLE {
    label 'process_high'
    errorStrategy 'terminate'
    afterScript ''
    publishDir "${params.output_directory}", mode: 'copy'
    input:
    path partitions, stageAs: 'partitions/*'
    val scopes
    val settings
    path helpers, stageAs: 'helpers/*'
    path reference_lengths, stageAs: 'reference/*'
    path read_evidence, stageAs: 'read_input'
    output:
    path('manual_review'), emit: bundle
    val(true), emit: ready
    script:
    // JSON is staged as content, rather than interpolated into shell commands.
    task.ext.reviewScopes = groovy.json.JsonOutput.toJson(scopes)
    task.ext.reviewSettings = groovy.json.JsonOutput.toJson(settings)
    def q = { v -> "'" + v.toString().replace("'", "'\"'\"'") + "'" }
    def lengths = reference_lengths ? "--reference-lengths ${q.call(reference_lengths)}" : ''
    """
    cat > scopes.json <<'LEAP_REVIEW_SCOPES'
    ${task.ext.reviewScopes}
    LEAP_REVIEW_SCOPES
    cat > settings.json <<'LEAP_REVIEW_SETTINGS'
    ${task.ext.reviewSettings}
    LEAP_REVIEW_SETTINGS
    python3 helpers/prepare_manual_review.py bundle partitions scopes.json settings.json manual_review ${lengths} --read-evidence read_input
    """
}
