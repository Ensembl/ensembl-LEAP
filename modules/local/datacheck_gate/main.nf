process DATACHECK_GATE {
    label 'process_very_low'
    errorStrategy 'terminate'
    afterScript ''

    input:
    path verdict
    path runner, stageAs: 'helpers/*'
    val review_ready

    output:
    path('datachecks.passed')

    script:
    """
    python3 '${runner}' gate --verdict '${verdict}'
    """
}
