//
// Subworkflow for the example pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GET_CURRENT_ANNOTATION } from '../../../modules/local/get_current_annotation/main.nf'
include { FILTER_PARTIAL_TRANSCRIPTS } from '../../../modules/local/filter_partial_transcripts/main.nf'
include { SPLIT_CURRENT_ANNOTATION_BY_STRAND } from '../../../modules/local/split_current_annotation_by_strand/main.nf'
include { EXTRACT_FIRST_AND_LAST_EXONS} from '../../../modules/local/extract_first_and_last_exons/main.nf'
include { EXTRACT_UNIQUE_SPLICE_SITES_FROM_CURRENT_ANNOTATION } from '../../../modules/local/extract_unique_splice_sites_from_current_annotation/main.nf'
include { FIND_NATURAL_OVERLAPS_IN_CURRENT_ANNOTATION } from '../../../modules/local/find_natural_overlaps_in_current_annotation/main.nf'
include { SPLIT_TRANSCRIPTS_BASED_ON_EXON_NUMBER } from '../../../modules/local/split_transcripts_based_on_exon_number/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PROCESS_ANNOTATION {
    take:
    reads_sample_sheet
    
    main:
    // Get the new anotation from the provided database and filter to remove partial transcripts
    GET_CURRENT_ANNOTATION(reads_sample_sheet.splitCsv(header: true))

    // Filter the annotation to remove partial transcripts
    FILTER_PARTIAL_TRANSCRIPTS(GET_CURRENT_ANNOTATION.out.current_annotation_by_chromosome)

    // Split the current annotation based on multi-exon and single-exon transcripts
    SPLIT_TRANSCRIPTS_BASED_ON_EXON_NUMBER(FILTER_PARTIAL_TRANSCRIPTS.out.current_annotation_ft)

    // Split by strand
    SPLIT_CURRENT_ANNOTATION_BY_STRAND(SPLIT_TRANSCRIPTS_BASED_ON_EXON_NUMBER.out.multi_exon_full_length_transcripts, "gff")

    // Extract first and last exon 
    EXTRACT_FIRST_AND_LAST_EXONS(SPLIT_CURRENT_ANNOTATION_BY_STRAND.out.positive_annotation.mix(SPLIT_CURRENT_ANNOTATION_BY_STRAND.out.negative_annotation))

    // Find naturally occuring overlaps between the first and last exons
    FIND_NATURAL_OVERLAPS_IN_CURRENT_ANNOTATION(EXTRACT_FIRST_AND_LAST_EXONS.out.first_exon_annotation.map { chromosome, file, strand, _exon_type -> ["${chromosome}_${strand}", file] }.join(EXTRACT_FIRST_AND_LAST_EXONS.out.last_exon_annotation.map { chromosome, file, strand, _exon_type -> ["${chromosome}_${strand}", file] }).join(SPLIT_TRANSCRIPTS_BASED_ON_EXON_NUMBER.out.single_exon_positive_no.mix(SPLIT_TRANSCRIPTS_BASED_ON_EXON_NUMBER.out.single_exon_negative_no)))

    // Extract unique splice sites
    EXTRACT_UNIQUE_SPLICE_SITES_FROM_CURRENT_ANNOTATION(EXTRACT_FIRST_AND_LAST_EXONS.out.first_exon_annotation.mix(EXTRACT_FIRST_AND_LAST_EXONS.out.last_exon_annotation))

    emit:
    stranded_annotation = SPLIT_CURRENT_ANNOTATION_BY_STRAND.out.positive_annotation.map { chromosome, file, strand, _format -> ["${chromosome}_${strand}", file] }.mix(SPLIT_CURRENT_ANNOTATION_BY_STRAND.out.negative_annotation.map { chromosome, file, strand, _format -> ["${chromosome}_${strand}", file] })
    natural_overlaps = FIND_NATURAL_OVERLAPS_IN_CURRENT_ANNOTATION.out.natural_overlaps
    exon_annotation = EXTRACT_FIRST_AND_LAST_EXONS.out.first_exon_annotation.map { chromosome, file, strand, exon_type -> ["${chromosome}_${strand}_${exon_type}", file] }.mix(EXTRACT_FIRST_AND_LAST_EXONS.out.last_exon_annotation.map { chromosome, file, strand, exon_type -> ["${chromosome}_${strand}_${exon_type}", file] })
    unique_splice_sites = EXTRACT_UNIQUE_SPLICE_SITES_FROM_CURRENT_ANNOTATION.out.unique_splice_sites.map { chromosome, file, strand, exon_type -> ["${chromosome}_${strand}_${exon_type}", file] }
    single_exon_beds = SPLIT_TRANSCRIPTS_BASED_ON_EXON_NUMBER.out.single_exon_positive.mix(SPLIT_TRANSCRIPTS_BASED_ON_EXON_NUMBER.out.single_exon_negative)
    do_not_extend_transcripts = GET_CURRENT_ANNOTATION.out.do_not_extend_transcripts_positive.mix(GET_CURRENT_ANNOTATION.out.do_not_extend_transcripts_negative)
    non_protein_coding_regions = FILTER_PARTIAL_TRANSCRIPTS.out.current_annotation_npc_negative.mix(FILTER_PARTIAL_TRANSCRIPTS.out.current_annotation_npc_positive).map { chromosome, strand, file -> ["${chromosome}_${strand}", file] }
    single_exon_gff = SPLIT_TRANSCRIPTS_BASED_ON_EXON_NUMBER.out.single_exon_positive_gff.mix(SPLIT_TRANSCRIPTS_BASED_ON_EXON_NUMBER.out.single_exon_negative_gff).map { chromosome, strand, file -> ["${chromosome}_${strand}", file] }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
