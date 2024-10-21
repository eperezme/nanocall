/*
 * Check input samplesheet and get read channels
 */

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv
    input_path

    main:
    /*
     * Check samplesheet is valid
     */
    SAMPLESHEET_CHECK ( samplesheet, input_path )
        .csv
        .splitCsv ( header:true, sep:',' )
        // .map { meta}
        .map { get_sample_info(it, params.genomes) }
        .set { ch_sample }



    // ch_sample.

    emit:
    ch_sample // [ sample, barcode, fasta, gtf, is_transcripts, annotation_str ]
}

// Function to resolve fasta and gtf file if using iGenomes
// Returns [ sample, input_file, barcode, fasta, gtf, is_transcripts, annotation_str, nanopolish_fast5 ]
def get_sample_info(LinkedHashMap sample, LinkedHashMap genomeMap) {
    def meta = [:]
    meta.id = sample.barcode
    meta.sample  = sample.sample
    meta.genome = sample.genome
    meta.flowcell_id = sample.flowcell_id
    def reads = [""]
    return [ meta, reads ]
}

// def add_unclassified(channel1) {
//     Channel.of(
//         meta.id = "unclassified",
//         meta.sample = "unclassified",
//         meta.flowcell_id = channel1[0].flowcell_id,
//         meta.genome = channel1[0].genome
//     ).view()
//     .set { added_unclassified }

//     return channel1.concat
// }


    //
    // THE GOOD ONE
    //
    // ch_sample.merge(ch_demuxed).map{ it ->
    //     def meta = it[0]
    //     def reads = it[3]
    //     return [meta, reads]
    // }.set { ch_demuxed }
