#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LEAP: LEAP Analysis pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : <TBA>
    Author: Ryan Merritt
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PROCESS_ANNOTATION } from './subworkflows/local/process_annotation/main.nf'
include { READ_PROCESSING } from './subworkflows/local/read_processing/main.nf'
include { PROCESS_CAGE_AND_POLYA } from './subworkflows/local/process_cage_and_polya/main.nf'
include { OVERLAP_READS } from './subworkflows/local/overlap_reads/main.nf'
include { DATACHECKS } from './subworkflows/local/datachecks/main.nf'
include { MANUAL_REVIEW } from './subworkflows/local/manual_review/main.nf'
include { LEAP_IMPACT_REPORT } from './subworkflows/local/impact_report/main.nf'
include { LEAP_IMPACT_GATE } from './modules/local/impact_gate/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    legacy_database_params = ['host', 'port', 'user', 'pass', 'dbname'].findAll { params.containsKey(it) }
    if (legacy_database_params)
        error 'Legacy database parameters are no longer supported. Use --db_secrets with a database connection JSON file.'
    // CLI/config values may be strings: Groovy treats the string "false" as true.
    datacheck_switch = params.run_datachecks.toString().toLowerCase()
    if (!(datacheck_switch in ['true', 'false']))
        error '--run_datachecks must be true or false.'
    datachecks_enabled = datacheck_switch == 'true'
    if (datachecks_enabled && !params.datacheck_reference_lengths) {
        error 'Automatic datachecks require --datacheck_reference_lengths pointing to a matching FASTA .fai or chromosome-length file. Use --run_datachecks false only to explicitly disable validation.'
    }
    if (datachecks_enabled) {
        reference_path = file(params.datacheck_reference_lengths, checkIfExists: true)
        if (!java.nio.file.Files.isRegularFile(reference_path) || !java.nio.file.Files.isReadable(reference_path))
            error 'Automatic datachecks require a readable reference-length file, not a directory.'
        datacheck_reference = Channel.value(reference_path)
    }

    if (!params.db_secrets)
        error 'Annotation export requires --db_secrets pointing to a readable database connection JSON file.'
    db_secrets_path = file(params.db_secrets, checkIfExists: true)
    if (!java.nio.file.Files.isRegularFile(db_secrets_path) || !java.nio.file.Files.isReadable(db_secrets_path))
        error '--db_secrets must be a readable regular JSON file.'
    // Report schema errors without echoing JSON fragments or connection values.
    def connection_valid = false
    try {
        def connection = new groovy.json.JsonSlurper().parseText(db_secrets_path.text)
        if (connection instanceof Map) {
            def database = connection.database ?: connection.dbname
            def password = connection.containsKey('password') ? connection.password : connection.get('pass', '')
            connection_valid = connection.host instanceof String && connection.host.trim() &&
                connection.user instanceof String && connection.user.trim() &&
                database instanceof String && database.trim() && password instanceof String &&
                (connection.port instanceof Integer || connection.port instanceof Long || connection.port instanceof String) && (connection.port.toString() ==~ /[0-9]+/) &&
                connection.port.toString().toBigInteger() >= 1 && connection.port.toString().toBigInteger() <= 65535 &&
                connection.get('dialect', 'mysql') == 'mysql' &&
                !(connection.containsKey('database') && connection.containsKey('dbname') && connection.database != connection.dbname) &&
                !(connection.containsKey('password') && connection.containsKey('pass') && connection.password != connection.pass)
        }
    } catch (Exception ignored) { connection_valid = false }
    if (!connection_valid)
        error 'Invalid database secrets JSON: require host, user, database (or dbname), port 1..65535, a string password (or pass), and mysql dialect.'
    db_secrets = Channel.value(db_secrets_path)

    review_switch = params.prepare_manual_review.toString().toLowerCase()
    if (!(review_switch in ['true', 'false']))
        error '--prepare_manual_review must be true or false.'
    review_enabled = review_switch == 'true'
    if (review_enabled && params.review_reference) {
        review_reference_path = file(params.review_reference, checkIfExists: true)
        if (!java.nio.file.Files.isRegularFile(review_reference_path) || !java.nio.file.Files.isReadable(review_reference_path))
            error '--review_reference must point to a readable FASTA or 2bit file.'
    }

    // Process provided bam file containing the reads
    READ_PROCESSING("${params.reads}", "${params.reads_index}")

    // Process the CAGE and polyA data
    PROCESS_CAGE_AND_POLYA("${params.cage}", "${params.polya}", READ_PROCESSING.out.reads_sample_sheet)

    // Grab the current annotation
    PROCESS_ANNOTATION(READ_PROCESSING.out.reads_sample_sheet, db_secrets)

    // Overlap read unique splice sites with cage and polyA data
    OVERLAP_READS(READ_PROCESSING.out.unique_splice_sites, PROCESS_CAGE_AND_POLYA.out.cage_ch, PROCESS_CAGE_AND_POLYA.out.polya_ch, PROCESS_ANNOTATION.out.unique_splice_sites, PROCESS_ANNOTATION.out.exon_annotation, PROCESS_ANNOTATION.out.natural_overlaps, PROCESS_ANNOTATION.out.stranded_annotation, PROCESS_ANNOTATION.out.single_exon_beds, PROCESS_ANNOTATION.out.do_not_extend_transcripts, PROCESS_ANNOTATION.out.non_protein_coding_regions, PROCESS_ANNOTATION.out.single_exon_gff, PROCESS_ANNOTATION.out.partial_annotation)

    // Derive scope independently of completed outputs for checks and review prep.
    expected_scope = READ_PROCESSING.out.reads_sample_sheet.splitCsv(header: true)
        .map { row -> row.chromosome.toString() }.unique()
        .flatMap { chromosome -> [[chromosome, 'positive'], [chromosome, 'negative']] }
    review_ready = Channel.value(true)
    if (review_enabled) {
        review_lengths = params.datacheck_reference_lengths ? Channel.value(file(params.datacheck_reference_lengths, checkIfExists: true)) : Channel.value([])
        review_settings = Channel.value([
            published_root: file("${params.output_directory}/manual_review").toAbsolutePath().toString(),
            reference_genome: params.review_reference ? review_reference_path.toAbsolutePath().toString() : '',
            genome_id: params.review_genome_id, species: params.review_species,
            jbrowser_pipeline: params.review_jbrowser_pipeline, ftp_project: params.review_ftp_project,
            public_base_url: params.review_public_base_url, ftp_target_dir: params.review_ftp_target_dir,
            hub_email: params.review_hub_email, publication: 'public',
            read_input_bam: params.reads ? file(params.reads).toAbsolutePath().toString() : '',
            read_input_index: params.reads_index ? file(params.reads_index).toAbsolutePath().toString() : ''
        ])
        MANUAL_REVIEW(expected_scope, PROCESS_ANNOTATION.out.raw_annotation,
            OVERLAP_READS.out.datacheck_outputs, review_lengths, review_settings,
            READ_PROCESSING.out.chromosome_alignments)
        review_ready = MANUAL_REVIEW.out.ready
    }
    if (datachecks_enabled) {
        DATACHECKS(expected_scope, PROCESS_ANNOTATION.out.raw_annotation,
            OVERLAP_READS.out.datacheck_outputs, OVERLAP_READS.out.selected_extensions,
            PROCESS_ANNOTATION.out.natural_overlaps, datacheck_reference, review_ready)
    } else {
        LEAP_IMPACT_REPORT(expected_scope, PROCESS_ANNOTATION.out.raw_annotation,
            OVERLAP_READS.out.datacheck_outputs, Channel.value([]))
        LEAP_IMPACT_GATE(LEAP_IMPACT_REPORT.out.summary, review_ready)
    }

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
