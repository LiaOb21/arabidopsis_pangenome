#!/usr/bin/env python3

import pandas as pd
import argparse

def matrix_prep(input_file, prefix, output_file):
    # upload the data
    node_mat = pd.read_csv(input_file, sep='\t')
    # drop unnecessary columns: path.name  path.length  path.step.count
    node_mat = node_mat.drop(columns=['path.name', 'path.length', 'path.step.count'])
    # sort the data by group.name
    node_mat = node_mat.sort_values(by='group.name')
    # Step 1: Group by 'group.name'
    # Step 2: Sum the node columns for each group
    # Step 3: Reset the index to get a clean DataFrame
    node_mat_aggregated = node_mat.groupby('group.name').sum().reset_index()
    # Renaming columns
    node_mat_aggregated.rename(columns={col: prefix + col if col.startswith('node') else col for col in node_mat_aggregated.columns}, inplace=True)
    #save the data
    node_mat_aggregated.to_csv(output_file, sep='\t', index=False)

    return node_mat_aggregated

def count_and_label(node_mat_aggregated, number_of_assemblies, softcore_limit):
    # For each node, count how many rows have a count greater than 0 - this will say how many assemblies have that node
    rows_greater_than_zero = (node_mat_aggregated.loc[:, node_mat_aggregated.columns.str.startswith('chr')] > 0).sum()
    # For each node, calculate the total sum of the column - this will say how many times that node appears in all assemblies
    total_sum = node_mat_aggregated.loc[:, node_mat_aggregated.columns.str.startswith('chr')].sum()
    # Create a new DataFrame with the information
    new_df = pd.DataFrame({
        'PAV': rows_greater_than_zero,
        'CNV': total_sum
    })
    # Optionally, reset the index to make the node names a column if they're not already
    new_df.reset_index(inplace=True)
    new_df.rename(columns={'index': 'Node_Name'}, inplace=True)

    # Create a new column 'Label' based on the count
    labels = []
    for i, row in new_df.iterrows():
        count = row['PAV']
        feature = row['Node_Name']
        if count > number_of_assemblies:
            raise ValueError(f"Count greater than {number_of_assemblies} for feature '{feature}'")
        elif count == number_of_assemblies:
            labels.append('core')
        elif softcore_limit <= count < number_of_assemblies:
            labels.append('softcore')
        elif count == 1:
            labels.append('private')
        else:
            labels.append('dispensable')
    new_df['Label'] = labels

    #save the data
    new_df.to_csv('PAV_CNV.csv', sep='\t', index=False)
    
    return new_df

def get_statistics(new_df):
    """Calculate statistics and return them."""
    # Calculate statistics
    # Number of features and percentage of annotated features in the pangenome compared to the total number of features in the GFF file
    total_features = len(new_df)

    # Number and percentage of core, softcore, dispensable and reference_private features calculated on annotated features
    core_features = len(new_df[new_df['Label'] == 'core'])
    core_percentage = (core_features / total_features) * 100

    softcore_features = len(new_df[new_df['Label'] == 'softcore'])
    softcore_percentage = (softcore_features / total_features) * 100

    dispensable_features = len(new_df[new_df['Label'] == 'dispensable'])
    dispensable_percentage = (dispensable_features / total_features) * 100

    private_features = len(new_df[new_df['Label'] == 'private'])
    private_percentage = (private_features / total_features) * 100

    # Create a dictionary to store the statistics
    stats = {
        'total_features': total_features,
        'core_features': core_features,
        'core_percentage': core_percentage,
        'softcore_features': softcore_features,
        'softcore_percentage': softcore_percentage,
        'dispensable_features': dispensable_features,
        'dispensable_percentage': dispensable_percentage,
        'private_features': private_features,
        'private_percentage': private_percentage,
    }

    # Convert the statistics to a DataFrame
    stats_df = pd.DataFrame([stats])

    # Save the statistics to a CSV file
    stats_df.to_csv('statistics.csv', sep='\t', index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='This script aggregates the node matrix data to have one row per assembly')
    parser.add_argument('-i', '--input', type=str, required=True, help='The node matrix file')
    parser.add_argument('-p', '--prefix', type=str, help='The custom prefix to add to node columns', required=True)
    parser.add_argument('-o', '--output', type=str, required=True, help='The output file name')
    parser.add_argument('-na', '--number_of_assemblies', type=int, required=True, help='The number of assemblies in the pangenome')
    parser.add_argument('-sl', '--softcore_limit', type=int, required=True, help='Number of assemblies above which a gene is considered softcore (sl< softcore < na). For example, if you have 10 assemblies in your pangenome and you want to consider genes that are present in 8-9 assemblies as softcore, you should use 8 as the value for this argument.')
    args = parser.parse_args()
    
    node_mat_aggregated = matrix_prep(args.input, args.prefix, args.output)
    new_df = count_and_label(node_mat_aggregated, args.number_of_assemblies, args.softcore_limit)
    get_statistics(new_df)