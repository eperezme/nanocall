process DORADO_ALIGNER {
    tag "$meta.id"
    label 'process_medium'
//TODO: This should be a process that is runed for each sample.
    // Define the environment using Conda
    conda "${moduleDir}/environment.yml"

    // Specify the container image based on the workflow's container engine
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE' :
        'docker.io/eperezme/dorado:full' }"

    // Define input channels: all POD5 files
    input:
    tuple val(meta), path(reads), path(index)
    // tuple val(imeta), path(index)

    // Define output channels: BAM files, FASTQ files, summary, and versions
    output:
    tuple val(meta), path("aligned/${meta.id}/*.bam"), emit: bam
    tuple val(meta), path("aligned/${meta.id}/*.bai"), emit: bai
    tuple val(meta), path("aligned/${meta.id}/alignment_summary.txt"), emit: summary
    path "versions.yml", emit: versions

    // Conditional execution based on task.ext.when
    when:
    task.ext.when == null || task.ext.when

    // Script section to perform the basecalling
    script:
    // Handle parameters and arguments
    def args = task.ext.args ?: ''

    """
    # Run the dorado aligner with the specified mode and arguments
    dorado aligner --output-dir aligned/${meta.id}/ --emit-summary ${index} ${reads}

    # Create a versions.yml file with the dorado version information
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dorado: \$(echo \$(dorado --version 2>&1) | sed -r 's/.{81}//')
    END_VERSIONS
    """
}

