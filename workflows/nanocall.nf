/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_nanocall_pipeline'

// Initialize summary parameters
def summary_params = paramsSummaryMap(params, workflow)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC            } from '../modules/nf-core/fastqc/main'
include { MULTIQC           } from '../modules/nf-core/multiqc/main'
include { TOULLIGQC         } from '../modules/nf-core/toulligqc/main'
include { FAST5_TO_POD5     } from '../modules/local/pod5/fast5_to_pod5/main'
include { DORADO_BASECALLER } from '../modules/local/dorado/dorado_basecaller/main'
include { DORADO_DEMUX      } from '../modules/local/dorado/dorado_demux/main'
include { DORADO_TRIM       } from '../modules/local/dorado/dorado_trim/main'
include { DORADO_CORRECT    } from '../modules/local/dorado/dorado_correct/main'
include { DORADO_ALIGNER      } from '../modules/local/dorado/dorado_aligner/main'

include { PIGZ_COMPRESS     } from '../modules/nf-core/pigz/compress/main'
include { INPUT_CHECK       } from '../subworkflows/local/input_check'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow NANOCALL {
    // Declare the input channel: tuple of (id, fast5_file)
    // Removed 'input:' block to make the workflow self-contained

    main:
    // Define the fast5 channel by reading all FAST5 files in the input directory
    Channel.fromPath("${params.input_path}/*.fast5", checkIfExists: true)
        .map { path ->
            def id = path.getBaseName()
            tuple(id, path)
        }
        .set { ch_fast5_files }

    ch_input = Channel.fromPath(params.input, checkIfExists: true)
    ch_input_path = Channel.fromPath(params.input_path, checkIfExists: true)
    // Initialize channels for versions and MultiQC
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // Check samplesheet
    INPUT_CHECK (ch_input, ch_input_path)
        .set { ch_sample }

    // Define the barcode channel
    ch_barcodes = ch_sample.map{ sample ->
        def meta = sample[0]
        def barcode = sample[0].id
        tuple(barcode, meta)
        }
    ch_indexes = ch_sample.map{ sample ->
        def meta = sample[0]
        def index = sample[0].genome
        tuple(meta, index)
        }



//////////////////////////////////////////////////////////////////////
//          MODULE: Run FAST5 to POD5 conversion
//////////////////////////////////////////////////////////////////////

    FAST5_TO_POD5 (
        ch_fast5_files
    )

    // Collect all POD5 files into a single directory
    FAST5_TO_POD5.out.pod5.collect().set { ch_pod5_folder }
    ch_versions = ch_versions.mix(FAST5_TO_POD5.out.versions)
    // ch_multiqc_files = ch_multiqc_files.mix(FAST5_TO_POD5.out.pod5)




//////////////////////////////////////////////////////////////////////
//              MODULE: Run Dorado basecaller
//////////////////////////////////////////////////////////////////////

    DORADO_BASECALLER (
        ch_pod5_folder
    )

    ch_bam_files = DORADO_BASECALLER.out.bam.collectFile(name: 'basecall.bam')
    ch_fastq_files = DORADO_BASECALLER.out.fastq.collectFile(name: 'basecall.fastq.gz')
    ch_summary_files = DORADO_BASECALLER.out.summary.collectFile(name: 'summary.tsv')
    ch_versions = ch_versions.mix(DORADO_BASECALLER.out.versions)

    if (params.fastq == true ) {
        ch_basecalled = ch_fastq_files
    } else {
        ch_basecalled = ch_bam_files
    }


//////////////////////////////////////////////////////////////////////
//              MODULE: Dorado demultiplexing
//////////////////////////////////////////////////////////////////////
    if (params.skip_demux == false) {
    DORADO_DEMUX (
        ch_basecalled
    )

    if (params.fastq) {ch_demuxed_files = DORADO_DEMUX.out.fastq}
    else {ch_demuxed_files = DORADO_DEMUX.out.bam}
    ch_demuxed_summary = DORADO_DEMUX.out.summary
    ch_versions = ch_versions.mix(DORADO_DEMUX.out.versions)


    ch_barcodes.join(
        ch_demuxed_files.flatten()
        .map{ path ->
            def baseName = path.getBaseName()
            def matcher = baseName =~ /_(barcode\d{2}|unclassified)$/
            def barcode = matcher ? matcher[0][1] : baseName
            tuple(barcode, path)
        }, remainder: true)
        .set { ch_named_demuxed}

    ch_sample.map{ sample ->
        def barcode = sample[0].id
        def meta = sample[0]
        tuple (barcode, meta)
        }
        .join(ch_named_demuxed)
        .map{ sample ->
            def meta = sample[1]
            def reads = sample[3]
            return [meta, reads]
            }
        .set { ch_demuxed }


//////////////////////////////////////////////////////////////////////
//              MODULE: PIGZ_COMPRESS
//////////////////////////////////////////////////////////////////////
    if (params.fastq == true) {
        PIGZ_COMPRESS (ch_demuxed)
    ch_demuxed = PIGZ_COMPRESS.out.archive
    ch_versions = ch_versions.mix(PIGZ_COMPRESS.out.versions)
    }

    } else {
        ch_demuxed = ch_basecalled.map{ read ->
        def meta = [:]
        meta.id = 'unclassified'
        meta.genome = 'no-genome'
        meta.sample = 'unclassified'
        def reads = read
        return [meta, reads]}
        // ch_demuxed.view()
    }

//////////////////////////////////////////////////////////////////////
//               MODULE: ToulligQC
//////////////////////////////////////////////////////////////////////
    ch_summary = ch_summary_files.mix(DORADO_DEMUX.out.summary)
    ch_summary.view()
    TOULLIGQC (
        ch_summary,
        ch_pod5_folder,
        ch_basecalled
    )
    ch_versions = ch_versions.mix(TOULLIGQC.out.versions)


//////////////////////////////////////////////////////////////////////
//                  MODULE: Dorado trimming
//////////////////////////////////////////////////////////////////////

    if (params.skip_trimming == false) {
    DORADO_TRIM (
        ch_demuxed
    )
    ch_trimmed = DORADO_TRIM.out.bam
    ch_versions = ch_versions.mix(DORADO_TRIM.out.versions)
    if (params.fastq) {
        ch_trimmed = DORADO_TRIM.out.fastq
    }} else {
        ch_trimmed = ch_demuxed
    }


//////////////////////////////////////////////////////////////////////
//              MODULE: FASTQC
//////////////////////////////////////////////////////////////////////

    FASTQC (
        ch_trimmed
    )

    ch_fastqc = FASTQC.out.html
    ch_fastqc_zip = FASTQC.out.zip
    ch_versions = ch_versions.mix(FASTQC.out.versions)



//////////////////////////////////////////////////////////////////////
//           MODULE: Dorado Correct
//////////////////////////////////////////////////////////////////////


    if (params.error_correction) {
    DORADO_CORRECT (
        ch_trimmed
    )
    ch_corrected = DORADO_CORRECT.out.fasta
    ch_versions = ch_versions.mix(DORADO_CORRECT.out.versions)
    // Add index to the input channel for alignment
    ch_prealign = ch_corrected.map{ sample ->
        def meta = sample[0]
        def reads = sample[1]
        def index = sample[0].genome
        return [meta, reads, index]
        }
    }

// Add index to the input channel for alignment
    if (!params.error_correction) {
    ch_prealign = ch_trimmed.map{ sample ->
        def meta = sample[0]
        def reads = sample[1]
        def index = sample[0].genome
        return [meta, reads, index]
        }
    }

//////////////////////////////////////////////////////////////////////
//              MODULE: Dorado Aligner
//////////////////////////////////////////////////////////////////////
    if (!params.skip_align) {
        DORADO_ALIGNER (
            ch_prealign
            )
        ch_aligned = DORADO_ALIGNER.out.bam
        ch_versions = ch_versions.mix(DORADO_ALIGNER.out.versions)
    }



    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_pipeline_software_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }




//////////////////////////////////////////////////////////////////////
//               MODULE: MultiQC
//////////////////////////////////////////////////////////////////////

    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(params, workflow) // Corrected function call
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))

    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
