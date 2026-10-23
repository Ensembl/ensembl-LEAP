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

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    // Process provided bam file containing the reads
    READ_PROCESSING("${params.reads}", "${params.reads_index}")

    // Process the CAGE and polyA data
    PROCESS_CAGE_AND_POLYA("${params.cage}", "${params.polya}", READ_PROCESSING.out.reads_sample_sheet)

    // Grab the current annotation
    PROCESS_ANNOTATION(READ_PROCESSING.out.reads_sample_sheet)

    // Overlap read unique splice sites with cage and polyA data
    OVERLAP_READS(READ_PROCESSING.out.unique_splice_sites, PROCESS_CAGE_AND_POLYA.out.cage_ch, PROCESS_CAGE_AND_POLYA.out.polya_ch, PROCESS_ANNOTATION.out.unique_splice_sites, PROCESS_ANNOTATION.out.exon_annotation, PROCESS_ANNOTATION.out.natural_overlaps, PROCESS_ANNOTATION.out.stranded_annotation, PROCESS_ANNOTATION.out.single_exon_beds, PROCESS_ANNOTATION.out.do_not_extend_transcripts, PROCESS_ANNOTATION.out.non_protein_coding_regions, PROCESS_ANNOTATION.out.single_exon_gff)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/