include { PREPARE_REVIEW_RAW } from '../../../modules/local/review_raw/main.nf'
include { PREPARE_REVIEW_EXTENSIONS } from '../../../modules/local/review_extensions/main.nf'
include { BUILD_REVIEW_BUNDLE } from '../../../modules/local/review_bundle/main.nf'
include { PREPARE_REVIEW_READS } from '../../../modules/local/review_reads/main.nf'

workflow MANUAL_REVIEW {
    take:
    expected_scope
    raw_annotation
    final_outputs
    reference_lengths
    settings
    read_alignments

    main:
    helpers = Channel.value([
        file("${moduleDir}/../../../bin/prepare_manual_review.py", checkIfExists: true),
        file("${moduleDir}/../../../bin/prepare_review_reads.py", checkIfExists: true),
        file("${moduleDir}/../../../assets/process_internal_extensions.py", checkIfExists: true),
        file("${moduleDir}/../../../assets/check_reference_coordinates.py", checkIfExists: true)
    ])
    inputs = expected_scope.toList().map { [scopes: it] }
        .combine(raw_annotation.toList().map { [raw: it] }).map { a, b -> a + b }
        .combine(final_outputs.toList().map { [finals: it] }).map { a, b -> a + b }
        .combine(read_alignments.toList().map { [reads: it] }).map { a, b -> a + b }
        .map { b ->
            def scopes = b.scopes
            def raw = b.raw
            def finals = b.finals
            def unique = scopes.collect { c, s -> [c.toString(), s.toString()] }
            if (!unique || unique.unique(false).size() != unique.size())
                error 'Manual review requires unique nonempty chromosome/strand scopes.'
            if (unique.any { c, s -> !(c ==~ /[A-Za-z0-9_.-]+/) || c in ['.', '..'] || !(s in ['positive', 'negative']) })
                error 'Invalid manual-review chromosome/strand scope.'
            def chromosomes = unique.collect { it[0] }.unique()
            def ids = unique.collect { "${it[0]}_${it[1]}".toString() }
            if (raw.size() != chromosomes.size() || raw.any { !(it[0].toString() in chromosomes) } ||
                finals.size() != ids.size() || finals.any { !(it[0].toString() in ids) })
                error 'Missing, duplicate or unexpected manual-review annotation partitions.'
            chromosomes.each { c ->
                if (raw.count { it[0].toString() == c } != 1)
                    error "Missing or duplicate raw annotation for manual review: ${c}"
            }
            ids.each { id ->
                if (finals.count { it[0].toString() == id } != 1)
                    error "Missing or duplicate final annotation for manual review: ${id}"
            }
            if (b.reads.size() != chromosomes.size() || b.reads.any { !(it[0].toString() in chromosomes) } ||
                chromosomes.any { c -> b.reads.count { it[0].toString() == c } != 1 })
                error 'Missing, duplicate or unexpected manual-review chromosome BAM/index pairs.'
            [scopes: unique, raw: raw, finals: finals, reads: b.reads]
        }
    raw_jobs = inputs.flatMap { b -> b.raw.collect { c, raw -> [c.toString(), raw] } }
    extension_jobs = inputs.flatMap { b ->
        b.scopes.collect { c, s ->
            def row = b.finals.find { it[0].toString() == "${c}_${s}" }
            [c, s, b.raw.find { it[0].toString() == c }[1], row[1], row[2]]
        }
    }
    PREPARE_REVIEW_RAW(raw_jobs, helpers)
    PREPARE_REVIEW_EXTENSIONS(extension_jobs, helpers)
    read_jobs = inputs.map { b -> [b.scopes, b.reads.collect { it[0].toString() },
        b.reads.collect { it[1] }, b.reads.collect { it[2] }] }
    read_helper = Channel.value(file("${moduleDir}/../../../bin/prepare_review_reads.py", checkIfExists: true))
    PREPARE_REVIEW_READS(read_jobs, read_helper)
    parts = PREPARE_REVIEW_RAW.out.partition.mix(PREPARE_REVIEW_EXTENSIONS.out.partition).toList()
    BUILD_REVIEW_BUNDLE(parts, inputs.map { it.scopes }, settings, helpers, reference_lengths,
        PREPARE_REVIEW_READS.out.evidence)

    emit:
    bundle = BUILD_REVIEW_BUNDLE.out.bundle
    ready = BUILD_REVIEW_BUNDLE.out.ready
}
