process SUMMARIZE_DATACHECKS {
    label 'process_very_low'
    errorStrategy 'terminate'
    afterScript ''
    publishDir "${params.output_directory}/datachecks", mode: 'copy', pattern: 'summary'

    input:
    path status_files, stageAs: 'statuses/*'
    val expected_manifest
    path runner, stageAs: 'helpers/*'

    output:
    path('summary'), emit: reports
    path('summary/verdict.json'), emit: verdict

    script:
    def manifest_json = groovy.json.JsonOutput.toJson(expected_manifest)
    """
    mkdir -p statuses
    cat > expected.json <<'MANIFEST'
    ${manifest_json}
    MANIFEST
    python3 '${runner}' summarize --expected expected.json --statuses statuses --output summary
    """
}
