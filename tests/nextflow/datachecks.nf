// A real local DSL2 smoke workflow, with no database, BAM processing or containers.
include { DATACHECKS } from '../../subworkflows/local/datachecks/main.nf'

params.fixture_dir = ''
params.output_directory = ''
params.missing_negative = false
params.reference_name = 'reference.sizes'

workflow {
    def root = params.fixture_dir
    scope = Channel.of(['1', 'positive'], ['1', 'negative'])
    raw = Channel.of(['1', file("${root}/raw_annotation/current_annotation_1.gff3", checkIfExists: true)])
    final_outputs = Channel.of('positive', 'negative')
        .filter { strand -> !(params.missing_negative && strand == 'negative') }
        .map { strand ->
            def id = "1_${strand}"
            [id.toString(), file("${root}/extended_transcripts/${id}_extended.gff", checkIfExists: true),
             file("${root}/extended_transcripts_stats/${id}_stats.tsv", checkIfExists: true),
             file("${root}/fixed_extensions/${id}_FE_protected.bed", checkIfExists: true),
             file("${root}/fixed_extensions/${id}_LE_protected.bed", checkIfExists: true),
             file("${root}/partial_protection/${id}_partial_protection.tsv", checkIfExists: true)]
        }
    selected = Channel.of('positive', 'negative').map { strand ->
        def id = "1_${strand}"
        [id.toString(), file("${root}/fixed_extensions_before_partial_protection/test_${id}_FE_SCF.bed", checkIfExists: true),
         file("${root}/fixed_extensions_before_partial_protection/test_${id}_LE_SCF.bed", checkIfExists: true)]
    }
    natural = Channel.of('positive', 'negative').map { strand ->
        def id = "1_${strand}"
        [id.toString(), file("${root}/natural_overlaps/${id}_NO.txt", checkIfExists: true)]
    }
    reference = Channel.value(file("${root}/${params.reference_name}", checkIfExists: true))
    DATACHECKS(scope, raw, final_outputs, selected, natural, reference, Channel.value(true))
}
