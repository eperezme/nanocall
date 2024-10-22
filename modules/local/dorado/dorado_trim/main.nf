process DORADO_TRIM {
    tag "$meta.id"
    label 'process_medium'

    // Define the environment using Conda
    conda "${moduleDir}/environment.yml"

    // Specify the container image based on the workflow's container engine
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE' :
        'docker.io/eperezme/dorado:full' }"


    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("trimmed/*.bam"), emit: bam, optional: true
    tuple val(meta), path("trimmed/*.fastq.gz"), emit: fastq, optional: true
    path "versions.yml", emit: versions

    // Conditional execution based on task.ext.when
    when:
    task.ext.when == null || task.ext.when

    script:
    // Handle parameters and arguments
    def args = task.ext.args ?: ''
    def additional_args = ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ftype = (params.fastq) ? "fastq" : "bam" // TODO: Check if this is correct
    if (params.fastq) {
        additional_args += "--emit-fastq "
    }

    """
    mkdir -p trimmed

    dorado trim ${additional_args} ${reads} > trimmed/${prefix}.$ftype

        # Create a versions.yml file with the dorado version information
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dorado: \$(echo \$(dorado --version 2>&1) | sed -r 's/.{81}//')
    END_VERSIONS
    """

}
