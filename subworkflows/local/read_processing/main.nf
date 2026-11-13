//
// Subworkflow for the processing of reads
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { EXTRACT_CHROMOSOMES_FROM_READS } from "../../../modules/local/extract_chromosomes_from_reads/main.nf"
include { SPLIT_READS_BY_CHROMOSOME } from "../../../modules/local/split_reads_by_chromosome/main.nf"
include { CONVERT_BAM_TO_BED } from "../../../modules/local/convert_bam_to_bed/main.nf"
include { CONVERT_BED_TO_GENEPRED } from "../../../modules/local/convert_bed_to_genepred/main.nf"
include { CONVERT_GENEPRED_TO_GTF } from "../../../modules/local/convert_genepred_to_gtf/main.nf"
include { SPLIT_CURRENT_ANNOTATION_BY_STRAND as SPLIT_READS_BY_STRAND } from "../../../modules/local/split_current_annotation_by_strand/main.nf"
include { EXTRACT_FIRST_AND_LAST_EXONS } from "../../../modules/local/extract_first_and_last_exons/main.nf"
include { EXTRACT_UNIQUE_SPLICE_SITES_FROM_READS } from "../../../modules/local/extract_unique_splice_sites_from_reads/main.nf"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow READ_PROCESSING {
    take:
    reads
    reads_index

    main:
    // Extract chromosomes from the provided bam file
    EXTRACT_CHROMOSOMES_FROM_READS(reads, reads_index)

    // Split bam file by chromosome
    SPLIT_READS_BY_CHROMOSOME(EXTRACT_CHROMOSOMES_FROM_READS.out.reads_chr_sample_sheet.splitCsv(header: true))

    // Convert bam to bed
    CONVERT_BAM_TO_BED(SPLIT_READS_BY_CHROMOSOME.out.split_reads_by_chromosome)

    // Convert bed to genepred
    CONVERT_BED_TO_GENEPRED(CONVERT_BAM_TO_BED.out.converted_bed)

    // Convery genepred to gtf
    CONVERT_GENEPRED_TO_GTF(CONVERT_BED_TO_GENEPRED.out.converted_genepred)

    // Split gtf by strand
    SPLIT_READS_BY_STRAND(CONVERT_GENEPRED_TO_GTF.out.converted_gtf, "gtf")

    // Extract the first and last exon
    EXTRACT_FIRST_AND_LAST_EXONS(SPLIT_READS_BY_STRAND.out.positive_annotation.mix(SPLIT_READS_BY_STRAND.out.negative_annotation))

    // Extract unique splice sites
    EXTRACT_UNIQUE_SPLICE_SITES_FROM_READS(EXTRACT_FIRST_AND_LAST_EXONS.out.first_exon_annotation.mix(EXTRACT_FIRST_AND_LAST_EXONS.out.last_exon_annotation))

    emit:
    reads_sample_sheet = EXTRACT_CHROMOSOMES_FROM_READS.out.reads_chr_sample_sheet
    chromosome_alignments = SPLIT_READS_BY_CHROMOSOME.out.split_reads_by_chromosome
    unique_splice_sites = EXTRACT_UNIQUE_SPLICE_SITES_FROM_READS.out.unique_splice_sites.map { chromosome, file, strand, exon_type -> ["${chromosome}_${strand}_${exon_type}", file] }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
