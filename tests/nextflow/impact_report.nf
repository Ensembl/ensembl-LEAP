include { LEAP_IMPACT_REPORT } from '../../subworkflows/local/impact_report/main.nf'
include { LEAP_IMPACT_GATE } from '../../modules/local/impact_gate/main.nf'
params.fixture_dir = ''
params.output_directory = ''
workflow {
    root = params.fixture_dir
    scope = Channel.of(['1', 'positive'], ['1', 'negative'])
    raw = Channel.of(['1', file("${root}/raw_annotation/current_annotation_1.gff3")])
    finals = Channel.of('positive', 'negative').map { s -> ["1_${s}".toString(),
        file("${root}/extended_transcripts/1_${s}_extended.gff"),
        file("${root}/extended_transcripts_stats/1_${s}_stats.tsv")] }
    LEAP_IMPACT_REPORT(scope, raw, finals, Channel.value([]))
    LEAP_IMPACT_GATE(LEAP_IMPACT_REPORT.out.summary, Channel.value(true))
}
