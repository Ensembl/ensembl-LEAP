process PREPARE_REVIEW_READS {
    label 'process_high'
    conda "${moduleDir}/environment.yml"
    errorStrategy 'terminate'
    afterScript ''
    input:
    tuple val(scopes), val(chromosomes), path(bams, stageAs: 'input_bams/*'), path(indexes, stageAs: 'input_indexes/*')
    path helper, stageAs: 'helpers/*'
    output:
    path('read_evidence'), emit: evidence
    script:
    def bamList = bams instanceof List ? bams : [bams]
    def indexList = indexes instanceof List ? indexes : [indexes]
    def inputs = chromosomes.withIndex().collect { c, n ->
        [chromosome: c, bam: bamList[n].toString(), index: indexList[n].toString()]
    }
    def manifest = groovy.json.JsonOutput.toJson([scopes: scopes, partitions: inputs])
    """
    cat > read-inputs.json <<'LEAP_REVIEW_READS'
    ${manifest}
    LEAP_REVIEW_READS
    python3 helpers/prepare_review_reads.py read-inputs.json read_evidence --threads ${task.cpus}
    """
}
