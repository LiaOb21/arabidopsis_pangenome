#!/usr/bin/env python3

import pandas as pd
import numpy as np
import argparse
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def load_and_preprocess(input_file):
    logging.info('Loading and preprocessing data')
    # Load the PAF file into a DataFrame
    df = pd.read_csv(input_file, sep='\t', header=None)

    # Drop unnecessary columns
    df_drp = df[[0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 13, 14]]

    # name the columns
    df_drp.columns = ['Path', 'Path_length', 'Path_start', 'Path_end', 'Path_strand', 'Feature', 'Feature_length', 'Feature_start', 'Feature_end', 'Alignment_length', 'Estimated_identity', 'Jaccard_index', 'Self_coverage']

    # Remove the prefixes
    df_drp.loc[:, 'Estimated_identity'] = df_drp['Estimated_identity'].str.replace('id:f:', '')
    df_drp.loc[:, 'Jaccard_index'] = df_drp['Jaccard_index'].str.replace('jc:f:', '')
    df_drp.loc[:, 'Self_coverage'] = df_drp['Self_coverage'].str.replace('sc:f:', '')

    # Convert the columns to numeric types
    df_drp.loc[:, 'Estimated_identity'] = pd.to_numeric(df_drp['Estimated_identity'])
    df_drp.loc[:, 'Alignment_length'] = pd.to_numeric(df_drp['Alignment_length'])
    df_drp.loc[:, 'Jaccard_index'] = pd.to_numeric(df_drp['Jaccard_index'])
    df_drp.loc[:, 'Self_coverage'] = pd.to_numeric(df_drp['Self_coverage'])

    return df_drp

def assign_block_ids(block):
    logging.info('Assigning block IDs')
    block_id = 0
    block['block_id'] = 0
    for i in range(1, len(block)):
        if block.loc[i, 'Path_start'] != block.loc[i-1, 'Path_end']:
            block_id += 1
        block.loc[i, 'block_id'] = block_id
    return block

def collapse_paf(df_drp, output_file):
    logging.info('Collapsing PAF')
    # Group the DataFrame by 'Path' and 'Feature'
    groups = df_drp.groupby(['Path','Feature'])

    # Create an empty DataFrame to store the results
    results = pd.DataFrame(columns=df_drp.columns.tolist())

    # Now, each group is a block
    for (path, feature), block in groups:
        # Sort the block by 'Feature_start' and 'Path_start'
        block = block.sort_values(by=['Path_start', 'Feature_start'])
        
        # Reset the index of the block
        block.reset_index(drop=True, inplace=True)

        # Assign block IDs
        block = assign_block_ids(block)
        
        # Process each sub-block separately
        for block_id, sub_block in block.groupby('block_id'):
            # Create a new row with the desired values
            new_row = {
                'Path': path,
                'Path_length': sub_block['Path_length'].iloc[0],
                'Path_start': sub_block['Path_start'].iloc[0],
                'Path_end': sub_block['Path_end'].iloc[-1],
                'Path_strand': sub_block['Path_strand'].iloc[0],
                'Feature': feature,
                'Feature_length': sub_block['Feature_length'].iloc[0],
                'Feature_start': sub_block['Feature_start'].iloc[0],
                'Feature_end': sub_block['Feature_end'].iloc[-1],
                'Alignment_length': sub_block['Alignment_length'].sum(),
                'Estimated_identity': np.average(sub_block['Estimated_identity'], weights=sub_block['Alignment_length']),
                'Jaccard_index': np.average(sub_block['Jaccard_index'], weights=sub_block['Alignment_length']),
                'Self_coverage': np.average(sub_block['Self_coverage'], weights=sub_block['Alignment_length'])
            }
            
            # Append the new row to the results DataFrame
            results = pd.concat([results, pd.DataFrame([new_row])], ignore_index=True)

    # Save the results to a file
    results.to_csv(output_file, sep='\t', index=False)

if __name__ == "__main__":
    logging.info('Starting the process...')
    parser = argparse.ArgumentParser(description='Process PAF file.')
    parser.add_argument('-i', '--input', type=str, required=True, help='Path to the input PAF file')
    parser.add_argument('-o', '--output', type=str, default='collapsed.paf', help='Path to the output file')

    args = parser.parse_args()

    df_drp = load_and_preprocess(args.input)
    collapse_paf(df_drp, args.output)
    logging.info('Done!:-)')