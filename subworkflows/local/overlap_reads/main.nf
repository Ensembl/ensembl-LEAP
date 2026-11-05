//
// Subworkflow for the overlap of reads with CAGE and polyA data
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { OVERLAP_READS_WITH_CAGE } from "../../../modules/local/overlap_reads_with_cage/main.nf"
include { OVERLAP_READS_WITH_POLYA } from "../../../modules/local/overlap_reads_with_polya/main.nf"
include { MATCH_POSSIBLE_EXTENSIONS_WITH_CURRENT_ANNOTATION } from "../../../modules/local/match_possible_extensions_with_current_annotation/main.nf"
include { EXTRACT_LONGEST_EXTENSION } from "../../../modules/local/extract_longest_extension/main.nf"
include { ADD_GENE_ID_TO_LONGEST_EXTENSION } from "../../../modules/local/add_gene_id_to_longest_extension/main.nf"
include { INTERSECT_FIRST_AND_LAST_EXON_EXTENSIONS } from "../../../modules/local/intersect_first_and_last_exon_extensions/main.nf"
include { CHECK_EXTENSION_SAFETY } from "../../../modules/local/check_extension_safety/main.nf"
include { FIX_UNSAFE_OVERLAPS } from "../../../modules/local/fix_unsafe_overlaps/main.nf"
include { EXTEND_TRANSCRIPTS } from "../../../modules/local/extend_transcripts/main.nf"
include { FILTER_OVERLAPPING_FEATURES } from "../../../modules/local/filter_overlapping_features/main.nf"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow OVERLAP_READS {
    take:
    reads
    cage
    polya
    current_annotation_uss
    exon_annotation
    natural_overlaps
    stranded_annotation
    single_exon_beds
    blacklist_transcripts
    non_protein_coding_regions
    single_exon_gff
    partial_annotation

    main:
    // Overlap reads with CAGE data
    OVERLAP_READS_WITH_CAGE(reads.join(cage))

    // Overlap reads with polyA data
    OVERLAP_READS_WITH_POLYA(reads.join(polya))

    // Overlap possible extenstions with current annotation
    MATCH_POSSIBLE_EXTENSIONS_WITH_CURRENT_ANNOTATION(OVERLAP_READS_WITH_CAGE.out.overlapped_reads.join(current_annotation_uss).mix(OVERLAP_READS_WITH_POLYA.out.overlapped_reads.join(current_annotation_uss)))

    // Extract the longest possible extension for each exon
    EXTRACT_LONGEST_EXTENSION(MATCH_POSSIBLE_EXTENSIONS_WITH_CURRENT_ANNOTATION.out.matched_possible_extensions)

    // Add gene id to the longest extension for overlap
    ADD_GENE_ID_TO_LONGEST_EXTENSION(EXTRACT_LONGEST_EXTENSION.out.longest_extension.join(exon_annotation))

    // Handling internal overlaps in the extensions
    FILTER_OVERLAPPING_FEATURES(ADD_GENE_ID_TO_LONGEST_EXTENSION.out.uss_with_gene_id.map { chromosome, strand, file -> ["${chromosome}_${strand}", file] }.groupTuple(size: 2).map { unique_id, files -> 
        def fe_file = files.find { it.name.contains('_FE_') }
        def le_file = files.find { it.name.contains('_LE_') } 
        [unique_id, fe_file, le_file] }.join(stranded_annotation).join(single_exon_gff).join(non_protein_coding_regions).join(natural_overlaps))

    // Intersect first and last exon to find safe and unsafe extensions
    INTERSECT_FIRST_AND_LAST_EXON_EXTENSIONS(FILTER_OVERLAPPING_FEATURES.out.filtered_extensions.join(single_exon_beds))

    // Check for safe and unsafe extensions
    CHECK_EXTENSION_SAFETY(INTERSECT_FIRST_AND_LAST_EXON_EXTENSIONS.out.overlapped_extensions.join(natural_overlaps))
    
    // Process unsafe extensions returning non-overlapping extensions
    FIX_UNSAFE_OVERLAPS(CHECK_EXTENSION_SAFETY.out.safety_checked_extensions)

    // Create the extensions gff per chromosome with stats
    EXTEND_TRANSCRIPTS(FIX_UNSAFE_OVERLAPS.out.fixed_extensions.join(stranded_annotation).join(blacklist_transcripts).join(partial_annotation).join(natural_overlaps))

    emit:
    EXTEND_TRANSCRIPTS.out.extended_transcripts
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
