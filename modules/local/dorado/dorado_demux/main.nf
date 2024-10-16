process DORADO_DEMUX {
    input:
    path "*.bam"
    output:
    path "*.bam", emit: bam, optional: true
    path "*.fastq.gz", emit: fastq, optional: true
    path "summary.tsv", emit: summary
    path "versions.yml", emit: versions



    script:
    def args = task.ext.args ?: ''

    def additional_args = ""
    if (params.sample_sheet) {
        additional_args += "--sample-sheet ${params.sample_sheet} "
    }
    if (params.kit) {
        additional_args += "--kit-name ${params.kit} "
    }
    if (params.barcode_both_ends) {
        additional_args += "--barcode-both-ends "
    }
    if (params.barcode_arrangement) {
        additional_args += "--barcode-arrangement "
    }
    if (params.barcode_sequences) {
        additional_args += "--barcode-sequences "
    }
    if (params.fastq) {
        additional_args += "--emit-fastq "
    }



    """
    #!/bin/bash
    dorado demux --no-trim $additional_args-o output_dir/ $input
    """"
}
