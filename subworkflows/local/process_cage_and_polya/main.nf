//
// Subworkflow for processing CAGE and polyA data
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SPLIT_CAGE_BY_CHROMOSOME_AND_STRAND } from "../../../modules/local/split_cage_by_chromosome_and_strand/main.nf"
include { SPLIT_POLYA_BY_CHROMOSOME_AND_STRAND } from "../../../modules/local/split_polya_by_chromosome_and_strand/main.nf"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PROCESS_CAGE_AND_POLYA {
    take:
    cage
    polya
    reads_sample_sheet

    main:
    // Split the cage data by chromosome and strand
    SPLIT_CAGE_BY_CHROMOSOME_AND_STRAND(cage, reads_sample_sheet)
    
    // Split the polyA data by chromosome
    SPLIT_POLYA_BY_CHROMOSOME_AND_STRAND(polya, reads_sample_sheet)

    emit:
    cage_ch = SPLIT_CAGE_BY_CHROMOSOME_AND_STRAND.out.cage_sample_sheet.splitCsv(skip: 1).map { chromosome, file, strand, exon_type -> ["${chromosome}_${strand}_${exon_type}", file, chromosome, strand] }
    polya_ch = SPLIT_POLYA_BY_CHROMOSOME_AND_STRAND.out.polya_sample_sheet.splitCsv(skip: 1).map { chromosome, file, strand, exon_type -> ["${chromosome}_${strand}_${exon_type}", file, chromosome, strand] }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
