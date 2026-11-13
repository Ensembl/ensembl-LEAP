include { MANUAL_REVIEW } from '../../subworkflows/local/manual_review/main.nf'
include { DATACHECKS } from '../../subworkflows/local/datachecks/main.nf'
params.fixture_dir = ''
params.output_directory = ''
params.missing_negative = false
params.with_checks = false
params.no_reference = false
params.missing_reads = false
workflow {
    root = params.fixture_dir
    scopes = Channel.of(['1', 'positive'], ['1', 'negative'])
    raw = Channel.of(['1', file("${root}/raw_annotation/current_annotation_1.gff3", checkIfExists: true)])
    finals = Channel.of('positive', 'negative').filter { !(params.missing_negative && it == 'negative') }
        .map { s -> ["1_${s}".toString(), file("${root}/extended_transcripts/1_${s}_extended.gff"),
            file("${root}/extended_transcripts_stats/1_${s}_stats.tsv"),
            file("${root}/fixed_extensions/1_${s}_FE_protected.bed"),
            file("${root}/fixed_extensions/1_${s}_LE_protected.bed"),
            file("${root}/partial_protection/1_${s}_partial_protection.tsv")] }
    lengths = params.no_reference ? Channel.value([]) : Channel.value(file("${root}/reference.sizes"))
    settings = Channel.value([genome_id: 'hg38', species: 'homo_sapiens', publication: 'public',
        published_root: file("${params.output_directory}/manual_review").toAbsolutePath().toString()])
    reads = params.missing_reads ? Channel.empty() : Channel.of(['1',
        file("${root}/reads/1.bam", checkIfExists: true), file("${root}/reads/1.bam.csi", checkIfExists: true)])
    MANUAL_REVIEW(scopes, raw, finals, lengths, settings, reads)
    if (params.with_checks) {
        selected = Channel.of('positive', 'negative').map { s -> ["1_${s}".toString(),
            file("${root}/fixed_extensions_before_partial_protection/test_1_${s}_FE_SCF.bed"),
            file("${root}/fixed_extensions_before_partial_protection/test_1_${s}_LE_SCF.bed")] }
        natural = Channel.of('positive', 'negative').map { s -> ["1_${s}".toString(), file("${root}/natural_overlaps/1_${s}_NO.txt")] }
        DATACHECKS(scopes, raw, finals, selected, natural, lengths, MANUAL_REVIEW.out.ready)
    }
}
