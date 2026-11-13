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

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
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

    // Process provided bam file containing the reads
    READ_PROCESSING("${params.reads}", "${params.reads_index}")

    // Process the CAGE and polyA data
    PROCESS_CAGE_AND_POLYA("${params.cage}", "${params.polya}", READ_PROCESSING.out.reads_sample_sheet)

    // Grab the current annotation
    PROCESS_ANNOTATION(READ_PROCESSING.out.reads_sample_sheet)

    // Overlap read unique splice sites with cage and polyA data
    OVERLAP_READS(READ_PROCESSING.out.unique_splice_sites, PROCESS_CAGE_AND_POLYA.out.cage_ch, PROCESS_CAGE_AND_POLYA.out.polya_ch, PROCESS_ANNOTATION.out.unique_splice_sites, PROCESS_ANNOTATION.out.exon_annotation, PROCESS_ANNOTATION.out.natural_overlaps, PROCESS_ANNOTATION.out.stranded_annotation, PROCESS_ANNOTATION.out.single_exon_beds, PROCESS_ANNOTATION.out.do_not_extend_transcripts, PROCESS_ANNOTATION.out.non_protein_coding_regions, PROCESS_ANNOTATION.out.single_exon_gff, PROCESS_ANNOTATION.out.partial_annotation)

    if (datachecks_enabled) {
        // Derive expected partitions from the selected read chromosomes, not from
        // completed GFF outputs: absent annotation/output channels must be detected.
        expected_scope = READ_PROCESSING.out.reads_sample_sheet.splitCsv(header: true)
            .map { row -> row.chromosome.toString() }.unique()
            .flatMap { chromosome -> [[chromosome, 'positive'], [chromosome, 'negative']] }
        DATACHECKS(expected_scope, PROCESS_ANNOTATION.out.raw_annotation,
            OVERLAP_READS.out.datacheck_outputs, OVERLAP_READS.out.selected_extensions,
            PROCESS_ANNOTATION.out.natural_overlaps, datacheck_reference)
    }

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
