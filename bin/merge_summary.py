# TODO - Create a function for this code
# Inputs: basecall_summary, barcoding_summary, folder with alignment summaries
# Outputs: final_summary

import os
import pandas as pd


def merge_summaries(basecall, barcoding, alignment=None):
    # load the basecall summary
    basecall_summary = pd.read_csv(basecall, sep="\t")
    # load the barcoding summary
    barcoding_summary = pd.read_csv(barcoding, sep="\t")

    merged_summary = pd.merge(
        basecall_summary.drop(columns=["barcode"]),
        barcoding_summary[["read_id", "barcode"]],
        on="read_id",
        how="left",
    )

    if alignment:
        # Directory containing the aligned files
        aligned_dir = alignment

        # List to store each alignment summary DataFrame
        alignment_summaries = []

        # Walk through the directory to find files ending with ".txt"
        for root, dirs, files in os.walk(aligned_dir):
            for file in files:
                if file.endswith("alignment_summary.txt"):
                    # Create the full path to the file
                    full_path = os.path.join(root, file)
                    # Read the file into a DataFrame
                    alignment_summary = pd.read_csv(full_path, sep="\t")
                    # Append the DataFrame to the list
                    alignment_summaries.append(alignment_summary)
        # Concatenate the alignment summaries
        concat_summary = pd.concat(alignment_summaries)
        # Drop the filename column from the concatenated alignment summaries
        concat_summary = concat_summary.drop(columns=["filename"])

        # Join the basecall summary with the alignment summary
        merged_summary = pd.merge(
            merged_summary, concat_summary, on="read_id", how="left"
        )

    return merged_summary


# This should be the input for the code
basecall = "basecall_summary.txt"
barcoding = "barcoding_summary.txt"
alignment = "alignment_summaries"

merged_summary = merge_summaries(basecall, barcoding, alignment)
merged_summary.to_csv("final_summary.txt", sep="\t", index=False)
