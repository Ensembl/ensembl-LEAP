process BUILD_LEAP_IMPACT_REPORT {
    label 'process_high'
    errorStrategy 'terminate'
    afterScript ''
    publishDir "${params.output_directory}", mode: 'copy', pattern: 'impact'

    input:
    val manifest
    path raw_files, stageAs: 'raw/*'
    path final_files, stageAs: 'final/*'
    path stats_files, stageAs: 'stats/*'
    path check_summary, stageAs: 'checks/*'
    path helpers, stageAs: 'helpers/*'

    output:
    path('impact'), emit: reports
    path('impact/summary.json'), emit: summary
    val(true), emit: ready

    script:
    def q = { v -> "'" + v.toString().replace("'", "'\"'\"'") + "'" }
    def manifest_json = groovy.json.JsonOutput.toJson(manifest)
    def checks = check_summary ? "--datacheck-summary ${q.call(check_summary)}" : ''
    """
    cat > impact_inputs.json <<'LEAP_IMPACT_INPUTS'
    ${manifest_json}
    LEAP_IMPACT_INPUTS
    python3 helpers/build_leap_impact_report.py --manifest impact_inputs.json \\
      --output impact --defer-errors ${checks}
    """
}
