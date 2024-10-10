process DORADO_BASECALLER {
    tag "$meta.id"
    label 'process_medium'

    // Specify the container to use for this process
    container "docker.io/eperezme/dorado:full"

    input:
    tuple val(meta), path(pod5_path)  // Input tuple containing metadata and path to pod5 file
    val dorado_model                  // Input value for the dorado model

    output:
    tuple val(meta), path("basecall*")   , emit: dorado_out  // Output tuple containing metadata and path to basecall files
    path "versions.yml"                  , emit: versions    // Output path for versions.yml file

    script:
    // Determine the mode based on the duplex parameter
    def mode = (params.duplex == true) ? "duplex" : "basecaller"

    // Initialize emit_args based on parameters
    def emit_args = ""
    if (params.error_correction == true || params.emit_fastq == true) {
        emit_args = "> basecall.fastq && gzip basecall.fastq"
    } 
    // Emit bam if emit_bam is true or modified_bases is true
    // Demultiplexing with kit name, output will be in bam
    elif (params.emit_bam == true || params.modified_bases || (params.demultiplexing && params.kit_name)) {
        emit_args = "> basecall.bam"
    }

    // Initialize additional_args based on parameters
    def additional_args = ""
    if (params.align) {
       additional_args += " --reference $params.ref_genome --mm2-opt '-k $params.kmer_size -w $params.win_size'" 
    }

    // Handle trimming options
    if (params.demultiplexing && params.kit) {
        additional_args += " --no-trim"
    } 
    elif (params.trim) {
        additional_args += " --trim $params.trim"
    }
    if (params.kit) {
        additional_args += " --kit-name $params.kit"
    }

    """
    // Run the dorado basecaller with the specified mode and arguments
    dorado $mode $dorado_model $additional_args $pod5_path $emit_args 

    // Create a versions.yml file with the dorado version information
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dorado: \$(echo \$(dorado --version 2>&1) | sed -r 's/.{81}//')
    END_VERSIONS
    """
}