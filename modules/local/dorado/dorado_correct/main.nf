process DORADO_CORRECT {
    tag "$meta.id"
    label 'process_high'

    // Define the environment using Conda
    conda "${moduleDir}/environment.yml"

    // Specify the container image based on the workflow's container engine
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE' :
        'docker.io/eperezme/dorado:full' }"

    input:
    tuple val(meta), path(reads)

    output:
    path "*.fasta", emit: fasta
    path "versions.yml", emit: versions

    // Conditional execution based on task.ext.when
    when:
    task.ext.when == null || task.ext.when


    script:
    // Handle parameters and arguments
    def args = task.ext.args ?: ''

    """
    # Run the dorado basecaller with the specified mode and arguments
    dorado correct ${reads} > ${meta.id}.fasta

    # Create a versions.yml file with the dorado version information
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dorado: \$(echo \$(dorado --version 2>&1) | sed -r 's/.{81}//')
    END_VERSIONS
    """
}


}
