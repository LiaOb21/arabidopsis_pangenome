#!/usr/bin/env python3

import pandas as pd
import numpy as np
import re
import argparse

def parse_ryo(exonerate_output):
    '''
    Parse exonerate output file (ryo) and return a DataFrame with the following columns:
    - Query_ID
    - Target_ID
    - Percent_Identity
    - Query_Length
    - Equivalenced_Mismatches
    - Percent_Similarity
    - Query_Alignment_Begin
    - Query_Alignment_End
    - Target_Alignment_Begin
    - Target_Alignment_End
    - Rank
    - Raw_Score
    - Alignment_Length (calculated)
    - Coverage (calculated)
    '''
    # open exonerate output file
    with open(exonerate_output, 'r') as file:
        lines = file.readlines()
    # Filter lines that start with 'UniRef100_'
    ryo_lines = [line for line in lines if line.startswith('UniRef100_')]
    # Create a DataFrame from the lines
    ryo = pd.DataFrame([line.split('\t') for line in ryo_lines])

    # Name columns
    ryo.columns = ['Query_ID', 'Target_ID', 'Percent_Identity', 'Query_Length', 'Equivalenced_Mismatches', 'Percent_Similarity', 'Query_Alignment_Begin', 'Query_Alignment_End', 'Target_Alignment_Begin', 'Target_Alignment_End', 'Rank', 'Raw_Score']
    # Remove the '\n' character from the last column
    ryo['Raw_Score'] = ryo['Raw_Score'].str.rstrip('\n')
    # Ensure that the columns are of the correct type
    ryo['Query_Alignment_End'] = ryo['Query_Alignment_End'].astype(int)
    ryo['Query_Alignment_Begin'] = ryo['Query_Alignment_Begin'].astype(int)
    ryo['Query_Length'] = ryo['Query_Length'].astype(int)
    # Calculate the length of the alignment (N.B. this length is calculated in number of amino acids, not nucleotides)
    ryo['Alignment_Length'] = ryo['Query_Alignment_End'] - ryo['Query_Alignment_Begin']
    # Calculate the coverage (N.B. this is calculated based on number of amino acids, not nucleotides)
    ryo['Coverage'] = (ryo['Alignment_Length'] / ryo['Query_Length']) * 100
    # Ensure that all the lines used for the merging are of the same type
    ryo['Raw_Score'] = ryo['Raw_Score'].astype(int)
    ryo['Target_Alignment_Begin'] = ryo['Target_Alignment_Begin'].astype(int)
    ryo['Target_Alignment_End'] = ryo['Target_Alignment_End'].astype(int)
    return ryo


def parse_c4_alignments_from_file(C4_alignments):
    with open(C4_alignments, 'r') as file:
        text = file.read()
    return parse_c4_alignments(text)


def parse_c4_alignments(text):
    '''
    Parse the text output of C4 alignments and return a list of dictionaries with the following keys:
    - Query
    - Target
    - Raw_score
    - Query_range
    - Target_range
    - Stop_codons
    '''
    # Split the text into blocks
    blocks = text.split("C4 Alignment:")
    codons = []

    for block in blocks[1:]:  # Skip the first block because it's empty
        # Extract the Query, Target, Raw score, Query range, and Target range
        query = re.search(r'Query: (.*)', block).group(1)
        target = re.search(r'Target: (.*)', block).group(1)
        raw_score = int(re.search(r'Raw score: (\d+)', block).group(1))
        query_range = re.search(r'Query range: (.*)', block).group(1)
        target_range = re.search(r'Target range: (.*)', block).group(1)
        # Count the number of stop codons
        # Count the total number of asterisks and divide by 3
        stop_codons = block.count('*') // 3
        codons.append({
            'Query': query,
            'Target': target,
            'Raw_score': raw_score,
            'Query_range': query_range,
            'Target_range': target_range,
            'Stop_codons': stop_codons
        })

    return pd.DataFrame(codons)


def modify_codons_df(codons):
    # Split the 'Query range' column
    codons[['Query_start', 'Query_end']] = codons['Query_range'].str.split(' -> ', expand=True)
    # Split the 'Target range' column
    codons[['Target_start', 'Target_end']] = codons['Target_range'].str.split(' -> ', expand=True)
    # Convert the new columns to int
    codons[['Query_start', 'Query_end', 'Target_start', 'Target_end']] = codons[['Query_start', 'Query_end', 'Target_start', 'Target_end']].astype(int)
    # Drop the 'Query range' and 'Target range' columns (they are duplicated)
    codons = codons.drop(['Query_range', 'Target_range'], axis=1)
    # Split the Target column on spaces to keep the revcomp info
    codons[['Target', 'Target_Orientation']] = codons['Target'].str.split(' ', expand=True)
    # Replace the [revcomp] with 'revcomp'
    codons['Target_Orientation'] = codons['Target_Orientation'].str.replace(r'[revcomp]', 'revcomp')
    # Fill NaN values in the Orientation column with 'forward'
    codons['Target_Orientation'] = codons['Target_Orientation'].fillna('forward')
    # Ensure that the lines used for the merging are of the same type
    codons['Raw_score'] = codons['Raw_score'].astype(int)
    codons['Query_start'] = codons['Query_start'].astype(int)
    codons['Query_end'] = codons['Query_end'].astype(int)
    codons['Target_start'] = codons['Target_start'].astype(int)
    codons['Target_end'] = codons['Target_end'].astype(int)
    return codons


def merge_dfs(ryo, codons):
    # Merge the two DataFrames
    merged = ryo.merge(codons, left_on=['Query_ID', 'Target_ID', 'Raw_Score', 'Query_Alignment_Begin', 'Query_Alignment_End', 'Target_Alignment_Begin', 'Target_Alignment_End'], right_on=['Query', 'Target', 'Raw_score', 'Query_start', 'Query_end', 'Target_start', 'Target_end'], how='inner')
    # Drop the redundant columns from the merged DataFrame
    merged = merged.drop(columns=['Query', 'Target', 'Raw_score', 'Query_start', 'Query_end', 'Target_start', 'Target_end'])
    # Split the 'target_id' column on the ':' character and expand the result into a DataFrame
    coordinates_df = merged['Target_ID'].str.split(':', expand=True)
    # The start and end coordinates are in the second column of coordinates_df
    # Split this column on the '-' character and expand the result into a DataFrame
    start_end_df = coordinates_df[1].str.split('-', expand=True)
    # Add the start and end coordinates as new columns in the merged DataFrame (these will be needed later)
    merged['Target_Start'] = start_end_df[0]
    merged['Target_End'] = start_end_df[1]
    # Convert 'Target_Start' and 'Target_End' to integers
    merged['Target_Start'] = merged['Target_Start'].astype(int)
    merged['Target_End'] = merged['Target_End'].astype(int)
    # Calculate target length
    merged['Target_Length'] = merged['Target_End'] - merged['Target_Start']
    # For forward orientation, check if Target_Alignment_Begin is 0 and Target_Alignment_End is Target_Length
    # For reverse orientation, check if Target_Alignment_Begin is Target_Length and Target_Alignment_End is 0
    merged['Is_Edge_Alignment'] = np.where(
        ((merged['Target_Orientation'] == 'forward') & (merged['Target_Alignment_Begin'] == 0) | (merged['Target_Alignment_End'] == merged['Target_Length'])) |
        ((merged['Target_Orientation'] == 'reverse') & (merged['Target_Alignment_Begin'] == merged['Target_Length']) | (merged['Target_Alignment_End'] == 0)),
        True, False)
    # Keep only rows where Coverage is equal or greater than 20 and where Percent_Identity is equal or greater than 50. We will consider these as valid alignments
    merged['Coverage'] = merged['Coverage'].astype(float)
    merged['Percent_Identity'] = merged['Percent_Identity'].astype(float)
    merged = merged[(merged['Coverage'] >= 20) & (merged['Percent_Identity'] >= 50)] # 20% coverage and 50% identity
    # Create a new column 'Classification' based on the coverage and stop codons to distinguish genes from pseudogenes
    merged['Classification'] = ['gene' if coverage == 100 and stop_codons == 0 else 'pseudogene' for coverage, stop_codons in zip(merged['Coverage'], merged['Stop_codons'])]
    # Create a new column 'UniRef100_Query_ID' by removing 'UniRef100_' prefix from 'Query_ID' column (needed for the conversion table later on)
    merged['UniRef100_Query_ID'] = merged['Query_ID'].str.replace('UniRef100_', '')
    return merged

def reclassify_ambiguous_results(merged):
    # Reclassify certain pseudogenes as 'unclassified'. This is necessary because for edge alignments We don't know if expanding further the region the gene would be complete
    for index, row in merged.iterrows():
        # Check if the row meets the conditions for reclassification
        if row['Classification'] == 'pseudogene' and row['Stop_codons'] == 0 and row['Is_Edge_Alignment']:
            # Update the classification to 'unclassified'
            merged.at[index, 'Classification'] = 'unclassified'
    return merged

def add_headers_info(merged, info_file):
    # Read the headers file into a DataFrame
    headers = pd.read_csv(info_file, sep='\t', header=None, names=['UniRef100_ID', 'Info'])
    
    # Merge the headers DataFrame with the existing DataFrame on 'UniRef100_Query_ID'
    merged = pd.merge(merged, headers, left_on='Query_ID', right_on='UniRef100_ID', how='left')
    # Drop the 'UniRef100_ID' column (it's redundant)
    merged = merged.drop(columns='UniRef100_ID')
    
    return merged

def load_conversion_table(conversion_table_file):
    # Load the conversion table into a DataFrame
    conversion_table = pd.read_csv(conversion_table_file, sep='\t', header=None, names=['UniProt_code', 'TAIR_locus', 'TAIR_gene'])
    # Split the UniProt_code column by "-", expanding into two columns
    conversion_table[['UniProt_code', 'UniProt_extension']] = conversion_table['UniProt_code'].str.split('-', expand=True) 
    # Group by UniProt_code and combine the unique values in the other columns
    conversion_table = conversion_table.groupby('UniProt_code').agg({
        'TAIR_locus': lambda x: ', '.join(set(x)),
        'TAIR_gene': lambda x: ', '.join(set(x)),
        'UniProt_extension': lambda x: ', '.join(set(x.dropna()))  # Exclude NaN values and join unique extensions
    }).reset_index()
    # save conversion table
    conversion_table.to_csv('conversion_table.tsv', index=False, sep='\t')
    return conversion_table

def add_TAIR_genes_to_merged(merged, conversion_table):
    # Merge the original DataFrame with the conversion table DataFrame
    uniref_merged = pd.merge(merged, conversion_table, left_on='UniRef100_Query_ID', right_on='UniProt_code', how='left')
    return uniref_merged

def save_sure_tair_rows(uniref_merged):
    # Extract rows where TAIR_gene is not null
    sure_tair_rows = uniref_merged[uniref_merged['TAIR_gene'].notnull()]
    # Save the rows where TAIR_gene is not null
    sure_tair_rows.to_csv('corresponding_to_tair.tsv', index=False, sep='\t')
    return sure_tair_rows

def save_not_tair_rows(uniref_merged):
    # Extract rows where TAIR_gene is null
    not_tair_rows = uniref_merged[uniref_merged['TAIR_gene'].isnull()]
    # Save the rows where TAIR_gene is null
    not_tair_rows.to_csv('new_results_not_tair.tsv', index=False, sep='\t')
    return not_tair_rows

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Parse exonerate and C4 alignment output for the non reference sequences extracted from the pangenome.')
    parser.add_argument('-e', '--exonerate_file', required=True, help='Exonerate output file to parse.')
    parser.add_argument('-c', '--c4_alignments_file', required=True, help='C4 alignments file to parse.')
    parser.add_argument('-t', '--conversion_table_file', required=True, help='Conversion table file to load.')
    parser.add_argument('-i', '--info_file', required=True, help='Information file obtained from uniref100 headers.')
    args = parser.parse_args()

    # call the functions
    df_ryo = parse_ryo(args.exonerate_file)
    df_c4 = parse_c4_alignments_from_file(args.c4_alignments_file)
    df_c4 = modify_codons_df(df_c4)
    conversion_table = load_conversion_table(args.conversion_table_file)
    merged_df = merge_dfs(df_ryo, df_c4)
    merged_df = reclassify_ambiguous_results(merged_df)
    merged_df = add_headers_info(merged_df, args.info_file)
    uniref_merged = add_TAIR_genes_to_merged(merged_df, conversion_table)
    sure_tair_rows = save_sure_tair_rows(uniref_merged)
    not_tair_rows = save_not_tair_rows(uniref_merged)