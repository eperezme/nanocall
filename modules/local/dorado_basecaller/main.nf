process DORADO_BASECALLER {
    tag "dorado_basecaller"
    label 'process_high'

    // Define the environment using Conda
    conda "${moduleDir}/environment.yml"

    // Specify the container image based on the workflow's container engine
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE' :
        'docker.io/eperezme/dorado:full' }"

    // Define input channels: all POD5 files
    input:
    path pod5_dir

    // Define output channels: BAM files, FASTQ files, summary, and versions
    output:
    path "*.bam", emit: bam
    path "*.fastq.gz", emit: fastq
    path "summary.tsv", emit: summary
    path "versions.yml", emit: versions

    // Conditional execution based on task.ext.when
    when:
    task.ext.when == null || task.ext.when

    // Script section to perform the basecalling
    script:
    // Handle parameters and arguments
    def args = task.ext.args ?: ''
    def mode = (params.duplex == true) ? "duplex" : "basecaller"
    def dorado_model = params.modified_bases ? "${params.model},${params.modified_bases}" : "${params.model}"

    // Initialize emit_args based on parameters
    def emit_args = ""
    if (params.error_correction == true || (params.emit_fastq == true && params.emit_bam == false)) {
        emit_args = "> basecall.fastq && gzip basecall.fastq"
    } 
    else if (params.emit_bam == true || params.modified_bases || (params.demultiplexing && params.kit)) {
        emit_args = "> basecall.bam"
    }

    // Initialize additional_args based on parameters
    def additional_args = ""
    if (params.align) {
       additional_args += " --reference ${params.ref_genome} --mm2-opt '-k ${params.kmer_size} -w ${params.win_size}'" 
    }
    // Handle trimming options
    if (params.demultiplexing && params.kit) {
        additional_args += " --no-trim"
    } 
    else if (params.trim) {
        additional_args += " --trim ${params.trim}"
    }
    if (params.kit) {
        additional_args += " --kit-name ${params.kit}"
    }

    """
    # Create output directory if it doesn't exist
    mkdir -p dorado_output

    # Run the dorado basecaller with the specified mode and arguments
    dorado ${mode} ${dorado_model} ${additional_args} ${pod5_dir} ${emit_args}

    # Move output files to the main directory
    mkdir -p dorado_output .
    mv dorado_output/*.bam .
    mv dorado_output/*.fastq.gz .
    mv dorado_output/summary.tsv .
    mv dorado_output/versions.yml .

    # Clean up
    rm -rf dorado_output
    """
}
