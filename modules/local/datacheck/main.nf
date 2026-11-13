process RUN_DATACHECK {
    tag "${chromosome}_${strand}_${check_name}"
    label 'process_medium'
    maxForks 4
    errorStrategy 'terminate'
    afterScript ''
    publishDir "${params.output_directory}/datachecks", mode: 'copy',
        saveAs: { filename -> filename.endsWith('.status.json') ? null : filename }

    input:
    tuple val(chromosome), val(strand),
        path(raw_files, stageAs: 'context/raw_annotation/*'),
        path(gff_files, stageAs: 'context/extended_transcripts/*'),
        path(stats_files, stageAs: 'context/extended_transcripts_stats/*'),
        path(before_files, stageAs: 'context/fixed_extensions_before_partial_protection/*'),
        path(after_files, stageAs: 'context/fixed_extensions/*'),
        path(audit_files, stageAs: 'context/partial_protection/*'),
        path(natural_files, stageAs: 'context/natural_overlaps/*')
    val check_name
    path check_scripts, stageAs: 'checkers/*'
    path runner, stageAs: 'helpers/*'
    path reference_lengths, stageAs: 'reference/*'

    output:
    path("${chromosome}_${strand}_${check_name}"), emit: reports
    tuple val(chromosome), val(strand), val(check_name),
        path("${chromosome}_${strand}_${check_name}/${chromosome}_${strand}_${check_name}.status.json"), emit: status

    script:
    def quote = { value -> "'" + value.toString().replace("'", "'\"'\"'") + "'" }
    """
    python3 ${quote.call(runner)} check \\
      --check ${quote.call(check_name)} --chromosome ${quote.call(chromosome)} --strand ${quote.call(strand)} \\
      --context context --checkers checkers --reference-lengths ${quote.call(reference_lengths)} \\
      --output ${quote.call("${chromosome}_${strand}_${check_name}")}
    """
}
