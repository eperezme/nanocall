import pandas as pd

# TODO - Create a function for this code
# Inputs: basecall_summary, barcoding_summary, folder with alignment summaries
# Outputs: final_summary


# Load basecalling summary file
basecall_summary = pd.read_csv("../basecall_summary.tsv", sep="\t")

# Load barcoding summary file
barcoding_summary = pd.read_csv("../barcoding_summary.txt", sep="\t")

# Load alignment summary file
alignment_summary = pd.read_csv("../aligned/barcode01/alignment_summary.txt", sep="\t")

# Save the final summary to a new file
# Substitute the barcode column of the basecall_summary.tsv with the barcode column of the barcoding summary for the same read_id
basecall_summary = pd.merge(
    basecall_summary.drop(columns=["barcode"]),
    barcoding_summary[["read_id", "barcode"]],
    on="read_id",
    how="left",
)
basecall_summary.to_csv("edited_summary_with_barcodes.tsv", sep="\t", index=False)

# Now we need to merge all the alignment summaries.
# for each alignment summary file, we need to concatenate the summary from the barcode1 with the barcode2, barcode3, etc.
# we will use the pandas library to do this
# First, let's concatenate the alignment summaries for barcode1, barcode2, barcode3, etc.

# For each alignment_summary.txt in the aligned directory we will concatenate the summary for each barcode
# list all .txt files in the aligned directory
import os

# Directory containing the aligned files
aligned_dir = "../aligned"

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
print(alignment_summaries)
concat_summary = pd.concat(alignment_summaries)
# Drop the filename column from the concatenated alignment summaries
concat_summary = concat_summary.drop(columns=["filename"])

# Join the basecall summary with the alignment summary
final_summary = pd.merge(basecall_summary, concat_summary, on="read_id", how="left")
final_summary.to_csv("final_summary.tsv", sep="\t", index=False)


# Now merge the alignment summary with the edited basecall summary, take
