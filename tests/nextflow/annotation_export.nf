include { GET_CURRENT_ANNOTATION } from '../../modules/local/get_current_annotation/main.nf'
params.db_secrets = ''
params.helpers = ''
params.output_directory = ''
workflow {
    helpers = Channel.value(file("${params.helpers}/*").toList())
    GET_CURRENT_ANNOTATION(Channel.of(['unused', 'unused', '21']),
        Channel.value(file(params.db_secrets, checkIfExists: true)), helpers)
}
