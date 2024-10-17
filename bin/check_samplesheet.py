#!/usr/bin/env python

import os
import sys
import errno
import argparse


def parse_args(args=None):
    Description = "Reformat nf-core/nanoseq samplesheet file and check its contents."
    Epilog = "Example usage: python check_samplesheet.py <FILE_IN> <FILE_OUT>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Input samplesheet file.")
    parser.add_argument("FILE_OUT", help="Output file.")
    return parser.parse_args(args)


def make_dir(path):
    if len(path) > 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise exception


def print_error(error, context="Line", context_str=""):
    error_str = f"ERROR: Please check samplesheet -> {error}"
    if context and context_str:
        error_str = f"ERROR: Please check samplesheet -> {error}\n{context}: '{context_str.strip()}'"
    print(error_str)
    sys.exit(1)


def check_samplesheet(file_in, file_out):
    """
    This function checks that the samplesheet follows the following structure:
    sample,flowcell_id,barcode,genome
    MCF7,FLO-MIN106,1,ASM584v2
    MCF8,FLO-MIN106,2,ASM584v2
    """

    sample_info_dict = {}
    with open(file_in, "r") as fin:
        # Check header
        HEADER = ["sample", "flowcell_id", "barcode", "genome"]
        header = fin.readline().strip().split(",")
        if header != HEADER:
            print_error(f"Header mismatch: {','.join(header)} != {','.join(HEADER)}")

        # Check sample entries
        for line in fin:
            if not line.strip():
                continue  # Skip empty lines
            lspl = [x.strip() for x in line.strip().split(",")]

            # Check valid number of columns per row
            if len(lspl) < len(HEADER):
                print_error(f"Invalid number of columns (minimum = {len(HEADER)})!", "Line", line)

            # Extract values
            sample, flowcell, barcode, genome = lspl[: len(HEADER)]

            # Check sample name entries
            if not sample or " " in sample:
                print_error("Sample entry is either missing or contains spaces!", "Line", line)

            # Check flowcell entry
            if not flowcell or " " in flowcell:
                print_error("Flowcell ID is either missing or contains spaces!", "Line", line)

            # Check barcode entry
            if not barcode.isdigit():
                print_error("Barcode entry is not an integer!", "Line", line)

            # Check genome entry
            if not genome or " " in genome:
                print_error("Genome entry is either missing or contains spaces!", "Line", line)

            # Create sample mapping dictionary = {sample: {flowcell : [ barcode, genome ]}}
            sample_info = [barcode, genome]
            if sample not in sample_info_dict:
                sample_info_dict[sample] = {}
            if flowcell not in sample_info_dict[sample]:
                sample_info_dict[sample][flowcell] = sample_info
            else:
                print_error("Same flowcell ID provided multiple times for the same sample!", "Line", line)

    # Write validated samplesheet with appropriate columns
    if sample_info_dict:
        out_dir = os.path.dirname(file_out)
        make_dir(out_dir)
        with open(file_out, "w") as fout:
            fout.write(",".join(["sample", "flowcell_id", "barcode", "genome"]) + "\n")
            for sample in sorted(sample_info_dict.keys()):
                for flowcell in sorted(sample_info_dict[sample].keys()):
                    formatted_barcode = f"barcode{sample_info_dict[sample][flowcell][0].zfill(2)}"
                    fout.write(",".join([sample, flowcell, formatted_barcode, sample_info_dict[sample][flowcell][1]]) + "\n")


def main(args=None):
    args = parse_args(args)
    check_samplesheet(args.FILE_IN, args.FILE_OUT)


if __name__ == "__main__":
    sys.exit(main())
