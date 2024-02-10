import pandas as pd
import argparse
import pysam
import os

# Set up argument parsing
parser = argparse.ArgumentParser(description="Generate and process coverage data from a BAM file.")
parser.add_argument("-i", "--input_bam", required=True, help="Input BAM file")
parser.add_argument("-out", "--output_csv", required=True, help="Output CSV file name")
args = parser.parse_args()

# Function to ensure BAM file is indexed
def ensure_bam_index(bam_file):
    index_file = bam_file + ".bai"
    if not os.path.exists(index_file):
        print("Index file not found. Creating index...")
        pysam.index(bam_file)
    else:
        print("Index file found.")

# Function to calculate coverage from BAM file
def calculate_coverage(bam_file):
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        coverage_data = []
        for pileupcolumn in bam.pileup():
            chrom = pileupcolumn.reference_name
            pos = pileupcolumn.pos + 1  # converting to 1-based
            depth = pileupcolumn.n
            coverage_data.append([chrom, pos, depth])
    return pd.DataFrame(coverage_data, columns=['chr', 'pos', 'depth'])

# Ensure the BAM file is indexed
ensure_bam_index(args.input_bam)

# Calculate coverage from the BAM file
coverage_df = calculate_coverage(args.input_bam)

# List to hold the rows before creating the DataFrame
rows_list = []

# Process each chromosome/scaffold separately
for chr_name, group in coverage_df.groupby('chr'):
    # Initialize start position and previous depth
    start_pos = None
    prev_depth = None
    for row in group.itertuples(index=False):
        # If this is the start of a new range or depth changes, save the previous range
        if start_pos is not None and (prev_depth != row.depth or row.pos != prev_pos + 1):
            rows_list.append({'chr': chr_name, 'start': start_pos, 'end': prev_pos, 'value1': prev_depth})
            start_pos = row.pos
        # If this is the start of a new contiguous range, update start_pos
        if start_pos is None:
            start_pos = row.pos
        prev_pos = row.pos
        prev_depth = row.depth
    # Save the last range for this chromosome/scaffold
    if start_pos is not None:
        rows_list.append({'chr': chr_name, 'start': start_pos, 'end': prev_pos, 'value1': prev_depth})

# Create DataFrame from the accumulated rows
result_df = pd.DataFrame(rows_list)

# Optional: Convert columns to appropriate data types
result_df['start'] = result_df['start'].astype(int)
result_df['end'] = result_df['end'].astype(int)
result_df['value1'] = result_df['value1'].astype(int)

# Check when end > start
result_df = result_df[result_df['end'] > result_df['start']]

# Save the result to a new file
result_df.to_csv(args.output_csv, sep='\t', index=False)
