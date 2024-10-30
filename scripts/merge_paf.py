#!/usr/bin/env python3

import pandas as pd
import argparse
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def merge_rows(df, distance):
    logging.info('Merging rows...')
    # Sort the DataFrame by 'Path', 'Feature', 'Path_start' and 'Path_end'
    df = df.sort_values(['Path', 'Feature', 'Path_strand', 'Path_start', 'Path_end'])

    # Initialize an empty list to store the rows
    rows = []

    # Initialize variables to store the current row and its details
    current_row = df.iloc[0].copy()
    current_path = current_row['Path']
    current_feature = current_row['Feature']
    current_strand = current_row['Path_strand']
    current_start = current_row['Path_start']
    current_end = current_row['Path_end']
    current_alignment_length = current_row['Alignment_length']
    current_estimated_identity = current_row['Estimated_identity'] * current_alignment_length
    current_jaccard_index = current_row['Jaccard_index'] * current_alignment_length
    current_self_coverage = current_row['Self_coverage'] * current_alignment_length

    # Iterate over the rows of the DataFrame
    for _, row in df.iloc[1:].iterrows():
        # If the current row can be merged with the previous row
        if row['Path'] == current_path and row['Feature'] == current_feature and row['Path_strand'] == current_strand and row['Path_start'] - current_end <= distance:
            # Update the end of the current row
            current_end = max(current_end, row['Path_end'])
            current_row['Path_end'] = current_end

            # Update the alignment_length and the weighted averages
            current_alignment_length += row['Alignment_length']
            current_estimated_identity += row['Estimated_identity'] * row['Alignment_length']
            current_jaccard_index += row['Jaccard_index'] * row['Alignment_length']
            current_self_coverage += row['Self_coverage'] * row['Alignment_length']
        else:
            # Calculate the final weighted averages
            current_row['Estimated_identity'] = current_estimated_identity / current_alignment_length
            current_row['Jaccard_index'] = current_jaccard_index / current_alignment_length
            current_row['Self_coverage'] = current_self_coverage / current_alignment_length

            # Append the current row to the list
            rows.append(current_row)

            # Update the current row and its details
            current_row = row.copy()
            current_path = row['Path']
            current_feature = row['Feature']
            current_strand = row['Path_strand']
            current_start = row['Path_start']
            current_end = row['Path_end']
            current_alignment_length = row['Alignment_length']
            current_estimated_identity = row['Estimated_identity'] * current_alignment_length
            current_jaccard_index = row['Jaccard_index'] * current_alignment_length
            current_self_coverage = row['Self_coverage'] * current_alignment_length

    # Calculate the final weighted averages for the last row
    current_row['Estimated_identity'] = current_estimated_identity / current_alignment_length
    current_row['Jaccard_index'] = current_jaccard_index / current_alignment_length
    current_row['Self_coverage'] = current_self_coverage / current_alignment_length

    # Append the last row to the list
    rows.append(current_row)

    # Convert the list of rows into a DataFrame
    results = pd.DataFrame(rows)

    return results

if __name__ == "__main__":
    logging.info('Starting merge_paf.py...')
    # Set up command-line argument parsing
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_file', required=True, help='The input PAF file')
    parser.add_argument('-o', '--output_file', required=True, help='The output file')
    parser.add_argument('-d', '--distance', type=int, required=True, help='The maximum distance between rows to be merged')
    args = parser.parse_args()

    # Load the input file into a DataFrame
    df = pd.read_csv(args.input_file, sep='\t')

    # Merge the rows of the DataFrame
    df = merge_rows(df, args.distance)

    # Save the DataFrame to the output file
    df.to_csv(args.output_file, sep='\t', index=False)
    logging.info('Results saved to %s', args.output_file)