The idea of this workflow is to process the data from the MinION sequencer. The data will be processed in the following steps:

1. Basecall the pod5 files
2. Demultiplex the basecalled data
3. Trim the adapters, barcodes and primers
4. Correct the FASTQ files to obtain FASTA (optional)
5. Align the reads to a reference genome

# Workflow

If you choose modified basecalling, the workflow will be with `BAM` files instead of `FASTQ` files.

If you choose to correct the FASTQ files, the workflow will be with `FASTA` files instead of `FASTQ` files.

For the diferents steps, we will use diferent modules:

### Basecalling

```bash
# POD5 IN
dorado basecaller model,modifications --no-trim <bam/fastq> pod5_dir/
# BAM/FASTQ OUT
```

### Demultiplexing

```bash
# BAM/FASTQ IN
dorado demux --no-trim --emit-summary --sort-bam --emit-fastq <fastq> --kit-name <kit> --sample-sheet <sample_sheet> --barcode-both-ends --barcode-arrangement --barcode-sequences -o output_dir/ <INPUT>
# Multiple BAM/FASTQ OUT
```

Here we should assign the files to the corresponding sample. The output will be `name_kit_barcodeXX.bam/fastq`. We can use the barcode XX to assign the files to the corresponding sample.

### Trimming

For each sample we do a worker.

```bash
# SINGLE BAM/FASTQ IN
dorado trim <INPUT> --emit-fastq --no-trim-primers --primer-sequences
# SINGLE BAM/FASTQ OUT
```

### Correcting (optional)

```bash
# SINGLE FASTQ IN
dorado correct $sample.fastq > $sample.fasta
# SINGLE FASTA OUT
```

### Aligning

```bash
# Multiple BED/FASTA IN
dorado aligner -o aligned/ --bed-file --emit-summary --mm2-opts "" index
# Multiple BAM OUT
```
