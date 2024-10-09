process FAST5_TO_POD5 {
    tag "$meta.id"
    label 'process_medium'

    // Define the environment using Conda
    conda "${moduleDir}/environment.yml"

    // Specify the container image based on the workflow's container engine
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pod5:0.3.15--pyhdfd78af_0' :
        'docker.io/eperezme/pod5:latest' }"

    // Define input channels: metadata and FAST5 reads
    input:
    tuple val(meta), path(reads)

    // Define output channels: POD5 files and versions file
    output:
    tuple val(meta), path("pod5/*_converted.pod5") , emit: pod5
    path "versions.yml"             , emit: versions

    // Conditional execution based on task.ext.when
    when:
    task.ext.when == null || task.ext.when

    // Script section to perform the conversion
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """

    
    # Convert FAST5 to POD5 using pod5 convert
    pod5 convert fast5 ${args} \
        -o pod5/${prefix}_converted.pod5 \
        ${reads}

    # Generate versions.yml if pod5 supports version reporting
    pod5 --version > versions.yml
    """
}
