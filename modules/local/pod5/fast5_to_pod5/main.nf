process FAST5_TO_POD5 {
    tag "$meta" // Use 'meta' directly as the tag
    label 'process_medium'

    // Define the environment using Conda
    conda "${moduleDir}/environment.yml"

    // Specify the container image based on the workflow's container engine
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pod5:0.3.15--pyhdfd78af_0' :
        'docker.io/eperezme/pod5:latest' }"

    // Define input channels: metadata (id) and FAST5 file
    input:
    tuple val(meta), path(fast5_file)

    // Define output channels: POD5 file and versions file
    output:
    path "pod5/${meta}_converted.pod5", emit: pod5
    path "versions_${meta}.yml", emit: versions

    // Conditional execution based on task.ext.when
    when:
    task.ext.when == null || task.ext.when

    // Script section to perform the conversion
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta}"

    """
    # Create output directory if it doesn't exist
    mkdir -p pod5

    # Convert FAST5 to POD5 using pod5 convert
    pod5 convert fast5 ${args} \
        -o pod5/${prefix}_converted.pod5 \
        ${fast5_file}

    # Generate versions.yml if pod5 supports version reporting
    pod5 --version > versions_${prefix}.yml
    """
}
