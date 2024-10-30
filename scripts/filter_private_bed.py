#!/usr/bin/env python3

import pandas as pd
import numpy as np
import argparse
import logging
import os

# Set up logging
logging.basicConfig(filename='filter_private_bed.log', filemode='w', level=logging.INFO)


def check_file_exists(filename):
    """Check if a file exists."""
    if not os.path.exists(filename):
        logging.error(f"File {filename} does not exist.")
        return False
    return True


def parse_args():
    """Parse the command line arguments."""
    # Create the parser
    parser = argparse.ArgumentParser(description='Inspect and filter the private genome in a pangenome graph.') 

    # Add the arguments
    parser.add_argument('-i', '--input_file', type=str, required=True, help='The input BED file containing non reference ranges, obtained with the command `odgi paths -i graph.og --non-reference-ranges reference.txt > non_reference.bed`')
    parser.add_argument('-thr', '--threshold', type=float, required=True, help='The length threshold (bp) for filtering the BED file, i.e. exclude sequences shorter than threshold.')
    parser.add_argument('-o', '--output_file', type=str, required=True, help='The output file name for the filtered BED file.')

    # Parse the arguments
    return parser.parse_args()

def process_bed(input_file, threshold, output_file):
    """ Process the BED file."""
    # Read the input file into a DataFrame
    df = pd.read_csv(input_file, sep='\t', header=0)
    # Calculate the length of each path
    df['path_length'] = df['end'] - df['start']

    # Define the boundaries and labels for the categories of path length - this is useful to understand how to set the threshold
    bins = [0, 1, 50, 99, 499, 999, 4999, 9999, 19999, 99999, 499999, 999999, np.inf]
    labels = ['1', '2-50', '51-99', '100-499', '500-999', '1000-4999', '5000-9999', '10000-19999', '20000-99999', '100000-499999', '500000-999999', 'from_1Mbp']

    # Create a new column 'size_category' that categorizes 'path_length'
    df['size_category'] = pd.cut(df['path_length'], bins=bins, labels=labels)

    # Count the number of rows in each category
    count_per_category = df.groupby('size_category').size()

    # Convert the Series to a DataFrame and reset the index
    count_df = count_per_category.reset_index(name='count')

    # Print the DataFrame to stdout
    print(count_df.to_string(index=False))

    # Log the DataFrame
    logging.info("\n" + count_df.to_string(index=False))

    # Filter the DataFrame
    df_filtered = df[df['path_length'] >= threshold]
    # Remove the 'path_length' and 'size_category' columns
    df_filtered = df_filtered.drop(['path_length', 'size_category'], axis=1)
    # Now df_filtered only contains rows where 'path_length' is equal to or greater than 100

    #save
    df_filtered.to_csv(output_file, index=False, sep='\t', header=False)

    # Filter out the lines that do not meet the threshold
    df_filtered_out = df[df['path_length'] < threshold]
    # Remove the 'path_length' and 'size_category' columns
    df_filtered_out = df_filtered_out.drop(['path_length', 'size_category'], axis=1)

    # Save the filtered out lines to a CSV file
    df_filtered_out.to_csv('filtered_out.bed', index=False, sep='\t', header=False)

def main():
    """Main function."""
    # Parse the command line arguments
    args = parse_args()

    # Check if the input file exists
    if not check_file_exists(args.input_file):
        return

    # Process the BED file
    logging.info("Processing BED file...")
    process_bed(args.input_file, args.threshold, args.output_file)
    logging.info("Done. Enjoy your pangenome! 🎉🎉🎉")

if __name__ == "__main__":
    main()