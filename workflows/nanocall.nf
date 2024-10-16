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

    // Initialize channels for versions and MultiQC
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // MODULE: Run FAST5 to POD5 conversion
    //
    FAST5_TO_POD5 (
        ch_fast5_files
    )

    // Collect all POD5 files into a single directory
    FAST5_TO_POD5.out.pod5.collect().set { ch_pod5_folder }
    ch_versions = ch_versions.mix(FAST5_TO_POD5.out.versions)
    // ch_multiqc_files = ch_multiqc_files.mix(FAST5_TO_POD5.out.pod5)

    //
    // MODULE: Run Dorado basecaller
    //
    DORADO_BASECALLER (ch_pod5_folder)
    DORADO_BASECALLER.out.bam.collectFile(name: 'basecall.bam').set { ch_bam_files }
    DORADO_BASECALLER.out.fastq.collectFile(name: 'basecall.fastq.gz').set { ch_fastq_files }
    DORADO_BASECALLER.out.summary.collectFile(name: 'summary.tsv').set { ch_summary_files }
    ch_versions = ch_versions.mix(DORADO_BASECALLER.out.versions)

    //
    // MODULE: ToulligQC
    //
    TOULLIGQC (
        ch_summary_files,
        ch_pod5_folder,
        ch_bam_files
    )
    //
    // MODULE: Run FastQC
    //
    // Uncomment and adjust if FastQC is needed
    /*
    FASTQC(id, FAST5_TO_POD5.out.pod5)
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{ it[1] })
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())
    */

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

    //
    // MODULE: MultiQC
    //
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
