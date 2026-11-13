include { BUILD_LEAP_IMPACT_REPORT } from '../../../modules/local/impact_report/main.nf'

workflow LEAP_IMPACT_REPORT {
    take:
    expected_scope
    raw_annotation
    final_outputs
    check_summary

    main:
    inputs = expected_scope.toList().map { [scopes: it] }
        .combine(raw_annotation.toList().map { [raw: it] }).map { a, b -> a + b }
        .combine(final_outputs.toList().map { [finals: it] }).map { a, b -> a + b }
    manifest = inputs.map { b ->
        def ids = b.scopes.collect { c, s -> "${c}_${s}".toString() }
        def chromosomes = b.scopes.collect { it[0].toString() }
        def issues = []
        [raw: b.raw, finals: b.finals].each { kind, rows ->
            rows.groupBy { it[0].toString() }.each { id, matches ->
                if (!(id in (kind == 'raw' ? chromosomes : ids)))
                    issues.add("Unexpected ${kind} scope: ${id}".toString())
                if (matches.size() > 1)
                    issues.add("Ambiguous ${kind} scope: ${id}".toString())
            }
        }
        [partitions: b.scopes.collect { c, s ->
            def finals = b.finals.findAll { it[0].toString() == "${c}_${s}" }
            [chromosome: c.toString(), strand: s.toString(),
             raw: b.raw.findAll { it[0].toString() == c.toString() }.collect { "raw/${it[1].name}".toString() },
             gff: finals.collect { "final/${it[1].name}".toString() },
             stats: finals.collect { "stats/${it[2].name}".toString() }]
        }, issues: issues]
    }
    helpers = Channel.value(['build_leap_impact_report.py', 'leap_impact_dashboard.html',
        'check_stats_gff_consistency.py', 'process_internal_extensions.py'].collect {
            file("${moduleDir}/../../../assets/${it}", checkIfExists: true)
        })
    BUILD_LEAP_IMPACT_REPORT(manifest,
        inputs.map { b -> b.raw.collect { it[1] }.unique() },
        inputs.map { b -> b.finals.collect { it[1] }.unique() },
        inputs.map { b -> b.finals.collect { it[2] }.unique() }, check_summary, helpers)

    emit:
    reports = BUILD_LEAP_IMPACT_REPORT.out.reports
    summary = BUILD_LEAP_IMPACT_REPORT.out.summary
    ready = BUILD_LEAP_IMPACT_REPORT.out.ready
}
