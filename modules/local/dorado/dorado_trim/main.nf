process DORADO_TRIM {
    tag "$meta.id"
    label 'process_medium'

    // Define the environment using Conda
    conda "${moduleDir}/environment.yml"

    // Specify the container image based on the workflow's container engine
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE' :
        'docker.io/eperezme/dorado:full' }"


    // Define input channels: all POD5 files
    input:
    tuple val(meta), path(reads)

    // Define output channels: BAM files, FASTQ files, summary, and versions
    output:
    tuple val(meta), path("*.bam"), emit: bam, optional: true
    tuple val(meta), path("*.fastq.gz"), emit: fastq, optional: true
    path "summary.tsv", emit: summary
    path "versions.yml", emit: versions

    // Conditional execution based on task.ext.when
    when:
    task.ext.when == null || task.ext.when

    // Script section to perform the basecalling
    script:
    // Handle parameters and arguments
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

}
