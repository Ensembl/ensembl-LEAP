// Run all checks from declared channel inputs, independently of publishDir copies.
include { RUN_DATACHECK as CHECK_BOUNDARIES } from '../../../modules/local/datacheck/main.nf'
include { RUN_DATACHECK as CHECK_INTERNAL_EXONS } from '../../../modules/local/datacheck/main.nf'
include { RUN_DATACHECK as CHECK_PRESERVATION } from '../../../modules/local/datacheck/main.nf'
include { RUN_DATACHECK as CHECK_STATS_GFF } from '../../../modules/local/datacheck/main.nf'
include { RUN_DATACHECK as CHECK_STRUCTURE } from '../../../modules/local/datacheck/main.nf'
include { RUN_DATACHECK as CHECK_COMPLETENESS } from '../../../modules/local/datacheck/main.nf'
include { RUN_DATACHECK as CHECK_PARTIAL_PROTECTION } from '../../../modules/local/datacheck/main.nf'
include { RUN_DATACHECK as CHECK_REFERENCE_COORDINATES } from '../../../modules/local/datacheck/main.nf'
include { RUN_DATACHECK as CHECK_CDS_OVERLAPS } from '../../../modules/local/datacheck/main.nf'
include { SUMMARIZE_DATACHECKS } from '../../../modules/local/datacheck_summary/main.nf'
include { DATACHECK_GATE } from '../../../modules/local/datacheck_gate/main.nf'

workflow DATACHECKS {
    take:
    expected_scope
    raw_annotation
    final_outputs
    selected_extensions
    natural_overlaps
    reference_lengths

    main:
    // toList emits [] even for an empty stream. Missing partitions must reach
    // completeness instead of disappearing through an inner join or collect.
    raw_context = raw_annotation.toList().map { rows -> [raw: rows] }
    final_context = final_outputs.toList().map { rows -> [finals: rows] }
    selected_context = selected_extensions.toList().map { rows -> [selected: rows] }
    natural_context = natural_overlaps.toList().map { rows -> [natural: rows] }
    context = raw_context.combine(final_context).map { a, b -> a + b }
        .combine(selected_context).map { a, b -> a + b }
        .combine(natural_context).map { a, b -> a + b }
    bundle = expected_scope.toList().map { rows -> [scopes: rows] }
        .combine(context).map { expected, inputs -> [scopes: expected.scopes, inputs: inputs] }

    jobs = bundle.flatMap { b ->
        b.scopes.collect { chromosome, strand ->
            if (!(chromosome.toString() ==~ /[A-Za-z0-9_.-]+/) || chromosome.toString() in ['.', '..'] || !(strand in ['positive', 'negative']))
                error("Invalid datacheck scope: ${chromosome} ${strand}")
            def id = "${chromosome}_${strand}".toString()
            def raw = b.inputs.raw.findAll { it[0].toString() == chromosome.toString() }
            def final_rows = b.inputs.finals.findAll { it[0].toString() == id }
            def selected = b.inputs.selected.findAll { it[0].toString() == id }
            def natural = b.inputs.natural.findAll { it[0].toString() == id }
            [chromosome.toString(), strand.toString(),
             raw.collect { it[1] }.unique(),
             final_rows.collect { it[1] }.unique(),
             final_rows.collect { it[2] }.unique(),
             selected.collectMany { [it[1], it[2]] }.unique(),
             final_rows.collectMany { [it[3], it[4]] }.unique(),
             final_rows.collect { it[5] }.unique(),
             natural.collect { it[1] }.unique()]
        }
    }
    expected_manifest = bundle.map { b ->
        def ids = b.scopes.collect { chromosome, strand -> "${chromosome}_${strand}".toString() }.toSet()
        def chromosomes = b.scopes.collect { it[0].toString() }.toSet()
        def issues = []
        [raw: b.inputs.raw, finals: b.inputs.finals, selected: b.inputs.selected, natural: b.inputs.natural].each { kind, rows ->
            rows.groupBy { it[0].toString() }.each { id, matches ->
                if (!(kind == 'raw' ? chromosomes.contains(id) : ids.contains(id)))
                    issues.add("Unexpected ${kind} scope: ${id}".toString())
                if (matches.size() > 1)
                    issues.add("Ambiguous ${kind} scope: ${id} (${matches.size()} records)".toString())
            }
        }
        [scopes: b.scopes, issues: issues]
    }

    check_scripts = Channel.value([
        'check_for_shorter_transcripts.py', 'check_internal_exon_positions.py',
        'check_transcript_preservation.py', 'check_stats_gff_consistency.py',
        'check_gff_structural_integrity.py', 'check_output_completeness.py',
        'check_partial_protection.py', 'check_reference_coordinates.py', 'check_cds_overlaps.py'
    ].collect { name -> file("${moduleDir}/../../../assets/${name}", checkIfExists: true) })
    runner = Channel.value(file("${moduleDir}/../../../bin/run_datacheck.py", checkIfExists: true))
    CHECK_BOUNDARIES(jobs, 'boundary', check_scripts, runner, reference_lengths)
    CHECK_INTERNAL_EXONS(jobs, 'internal_exons', check_scripts, runner, reference_lengths)
    CHECK_PRESERVATION(jobs, 'preservation', check_scripts, runner, reference_lengths)
    CHECK_STATS_GFF(jobs, 'stats_gff', check_scripts, runner, reference_lengths)
    CHECK_STRUCTURE(jobs, 'structure', check_scripts, runner, reference_lengths)
    CHECK_COMPLETENESS(jobs, 'completeness', check_scripts, runner, reference_lengths)
    CHECK_PARTIAL_PROTECTION(jobs, 'partial_protection', check_scripts, runner, reference_lengths)
    CHECK_REFERENCE_COORDINATES(jobs, 'reference_coordinates', check_scripts, runner, reference_lengths)
    CHECK_CDS_OVERLAPS(jobs, 'cds_overlaps', check_scripts, runner, reference_lengths)

    statuses = CHECK_BOUNDARIES.out.status.mix(
        CHECK_INTERNAL_EXONS.out.status, CHECK_PRESERVATION.out.status, CHECK_STATS_GFF.out.status,
        CHECK_STRUCTURE.out.status, CHECK_COMPLETENESS.out.status, CHECK_PARTIAL_PROTECTION.out.status,
        CHECK_REFERENCE_COORDINATES.out.status, CHECK_CDS_OVERLAPS.out.status)
        .map { chromosome, strand, check_name, result -> result }.toList()
    SUMMARIZE_DATACHECKS(statuses, expected_manifest, runner)
    DATACHECK_GATE(SUMMARIZE_DATACHECKS.out.verdict, runner)

    emit:
    reports = CHECK_BOUNDARIES.out.reports.mix(
        CHECK_INTERNAL_EXONS.out.reports, CHECK_PRESERVATION.out.reports, CHECK_STATS_GFF.out.reports,
        CHECK_STRUCTURE.out.reports, CHECK_COMPLETENESS.out.reports, CHECK_PARTIAL_PROTECTION.out.reports,
        CHECK_REFERENCE_COORDINATES.out.reports, CHECK_CDS_OVERLAPS.out.reports)
    summary = SUMMARIZE_DATACHECKS.out.reports
    passed = DATACHECK_GATE.out
}
