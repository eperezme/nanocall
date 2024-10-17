process DORADO_DEMUX {
    tag "dorado_demux"
    label 'process_medium'

    // Define the environment using Conda
    conda "${moduleDir}/environment.yml"

    // Specify the container image based on the workflow's container engine
        container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE' :
        'docker.io/eperezme/dorado:full' }"

    //Define input channels: all BAM files
    input:
    path(bamfile)
    output:
    path "demuxed/*.bam", emit: bam, optional: true
    path "demuxed/*.fastq", emit: fastq, optional: true
    path "demuxed/barcoding_summary.txt", emit: summary
    path "versions.yml", emit: versions



    script:
    def args = task.ext.args ?: ''

    // Initialize additional_args based on parameters
    def additional_args = ""
    // Sample_sheet to restrict the barcode classifications to only those present, and apply aliases to the detected classifications
    if (params.sample_sheet) {
        additional_args += "--sample-sheet ${params.sample_sheet} "
    }

    // Kit name to use for demultiplexing
    if (params.kit) {
        additional_args += "--kit-name ${params.kit} "
    }
    // Require both ends of a read to be barcoded for a double ended barcode.
    if (params.barcode_both_ends) {
        additional_args += "--barcode-both-ends "
    }
    //
    if (params.barcode_arrangement) {
        additional_args += "--barcode-arrangement ${params.barcode_arrangement} "
    }
    if (params.barcode_sequences) {
        additional_args += "--barcode-sequences ${params.barcode_sequences} "
    }
    //  A file with a newline-delimited list of reads to demux.
    if (params.read_ids) {
        additional_args += "--read-ids ${params.read_ids} "
    }
    // Emit a FASTQ file for each barcode.
    if (params.fastq) {
        additional_args += "--emit-fastq "
    }

    """
    #!/bin/bash
    dorado demux --no-trim --emit-summary ${additional_args} -o demuxed/ ${bamfile}

    # Create a versions.yml file with the dorado version information
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dorado: \$(echo \$(dorado --version 2>&1) | sed -r 's/.{81}//')
    END_VERSIONS
    """
}
