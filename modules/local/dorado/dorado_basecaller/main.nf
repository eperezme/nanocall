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
    path (pod5_dir, stageAs: 'pod5_dir/*')

    // Define output channels: BAM files, FASTQ files, summary, and versions
    output:
    path "*.bam", emit: bam, optional: true
    path "*.fastq.gz", emit: fastq, optional: true
    path "summary.tsv", emit: summary
    path "versions.yml", emit: versions

    // Conditional execution based on task.ext.when
    when:
    task.ext.when == null || task.ext.when

    // Script section to perform the basecalling
    script:
    // Handle parameters and arguments
    def args = task.ext.args ?: ''

    // Define the basecalling method
    def mode = (params.duplex == true) ? "duplex" : "basecaller"

    // Define the model structure to be passed to dorado
    def dorado_model = params.modified_bases ? "${params.model},${params.modified_bases}" : "${params.model}"

    // Initialize emit_args and additional_args based on parameters
    def emit_args = "> basecall.bam"
    def additional_args = ""
    def outfile = "basecall.bam"

    if (params.error_correction == true || (params.emit_fastq == true)) {
        emit_args = "> basecall.fastq"
        outfile = "basecall.fastq"
        additional_args += " --emit-fastq"
    }
    """
    # Run the dorado basecaller with the specified mode and arguments
    dorado ${mode} ${dorado_model} ${additional_args} --no-trim pod5_dir/ ${emit_args}

    # Create the summary file
    dorado summary $outfile > summary.tsv

    # gzip the fastq file if it exists
    if [ -f basecall.fastq ]; then
        gzip -k basecall.fastq
    fi

    # Create a versions.yml file with the dorado version information
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dorado: \$(echo \$(dorado --version 2>&1) | sed -r 's/.{81}//')
    END_VERSIONS
    """
}

