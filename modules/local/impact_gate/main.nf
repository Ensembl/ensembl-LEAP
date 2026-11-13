// Used only when the biological datachecks have been explicitly disabled.
process LEAP_IMPACT_GATE {
    label 'process_very_low'
    errorStrategy 'terminate'
    afterScript ''
    input:
    path summary
    val review_ready
    output:
    path('impact.complete')
    script:
    """
    python3 - <<'PY'
    import json
    from pathlib import Path
    data = json.loads(Path('${summary}').read_text())
    if data['status'] != 'COMPLETE' or data['exit_code'] != 0:
        raise SystemExit('Impact statistics are incomplete; see impact/diagnostics.csv. Datachecks were disabled.')
    Path('impact.complete').write_text('Impact statistics complete; biological datachecks were disabled.\\n')
    PY
    """
}
