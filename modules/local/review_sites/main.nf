process PREPARE_REVIEW_SITES {
    label 'process_high'
    errorStrategy 'terminate'
    afterScript ''
    input:
    val scopes
    path cage, stageAs: 'cage_input/*'
    path polya, stageAs: 'polya_input/*'
    path helpers, stageAs: 'helpers/*'
    path reference_lengths, stageAs: 'reference/*'
    output:
    path('site_evidence'), emit: evidence
    script:
    def manifest = groovy.json.JsonOutput.toJson([scopes:scopes,sources:[cage:cage.toString(),polya:polya.toString()]])
    def q = { v -> "'" + v.toString().replace("'", "'\"'\"'") + "'" }
    def lengths = reference_lengths ? "--reference-lengths ${q.call(reference_lengths)}" : ''
    """
    cat > site-inputs.json <<'LEAP_REVIEW_SITES'
    ${manifest}
    LEAP_REVIEW_SITES
    python3 helpers/prepare_review_sites.py site-inputs.json site_evidence ${lengths}
    """
}
