#!/usr/bin/env python3

import pandas as pd
import argparse
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def untangle_paf(untangle_paf_file):
    logging.info('Processing file: %s', untangle_paf_file)
    # read in the paf file -this file must be "filter_passed_features.csv" coming from core_dispensable_genes.py
    untangle_paf = pd.read_csv(untangle_paf_file, sep='\t')
    # remove all what comes after ':' in feature column
    untangle_paf['Feature'] = untangle_paf['Feature'].str.split(':').str[0]
    # add column called 'Source' to the dataframe and set all values to 'odgi_untangle' 
    untangle_paf['Source'] = 'odgi_untangle'
    # drop columns 'Jaccard_index', 'Self_coverage'
    untangle_paf = untangle_paf.drop(columns=['Jaccard_index', 'Self_coverage'])
    return untangle_paf

def untangle_paf_attributes(pangenome_screening_file):
    logging.info('Processing file: %s', pangenome_screening_file)
    # read in the pangenome screening file - this file must be "pangenome_screening.csv" coming from core_dispensable_genes.py
    pangenome_screening = pd.read_csv(pangenome_screening_file, sep='\t')
    # select the columns 'Feature' and 'Attributes
    untangle_attributes = pangenome_screening[['Feature', 'Attributes']]
    return untangle_attributes


def exonerate_results_tair(exonerate_results_tair_file):
    logging.info('Processing file: %s', exonerate_results_tair_file)
    # read in the exonerate results file - this file must be "reviewed_corresponding_to_tair.tsv" coming from parse_exonerate.py and reviewed with review_exonerate_results.py
    exonerate_results_tair = pd.read_csv(exonerate_results_tair_file, sep='\t')
    # create a assembly column by splitting Target_ID at '#'
    exonerate_results_tair.loc[:, 'Assembly'] = exonerate_results_tair['Target_ID'].str.split('#').str[0]
    # create a 'Path' column by splitting Target_ID at ':'
    exonerate_results_tair.loc[:, 'Path'] = exonerate_results_tair['Target_ID'].str.split(':').str[0]
    #reorder columns
    exonerate_results_tair = exonerate_results_tair[['Assembly', 'Path', 'Target_Length', 'Target_Start', 'Target_End', 'Target_Orientation', 'TAIR_gene', 'Query_Length', 'Query_Alignment_Begin', 'Query_Alignment_End', 'Alignment_Length', 'Percent_Identity', 'Coverage', 'Classification', 'Info']]
    #rename columns
    exonerate_results_tair.columns = ['Assembly', 'Path', 'Path_length', 'Path_start', 'Path_end', 'Path_strand', 'Feature', 'Feature_length', 'Feature_start', 'Feature_end', 'Alignment_length', 'Estimated_identity', 'Coverage', 'Classification', 'Attributes']
    # keep only rows where Classification is 'gene'
    exonerate_tair_genes = exonerate_results_tair[exonerate_results_tair['Classification'] == 'gene']
    # drop column 'Classification'
    exonerate_tair_genes = exonerate_tair_genes.drop(columns=['Classification'])
    # add column called 'Source' to the dataframe and set all values to 'exonerate' 
    exonerate_tair_genes['Source'] = 'exonerate'
    # create a new df 'exonerate_tair_genes_attributes' by selecting the columns 'Feature' and 'Attributes'
    exonerate_tair_genes_attributes = exonerate_tair_genes[['Feature', 'Attributes']]
    # drop column 'Attributes'
    exonerate_tair_genes = exonerate_tair_genes.drop(columns=['Attributes'])
    return exonerate_tair_genes, exonerate_tair_genes_attributes

def exonerate_results_new(exonerate_results_new_file):
    logging.info('Processing file: %s', exonerate_results_new_file)
    # read in the exonerate results file - this file must be "reviewed_not_tair_results.tsv" coming from parse_exonerate.py and reviewed with review_exonerate_results.py
    exonerate_results_new = pd.read_csv(exonerate_results_new_file, sep='\t')
    # create a assembly column by splitting Target_ID at '#'
    exonerate_results_new.loc[:, 'Assembly'] = exonerate_results_new['Target_ID'].str.split('#').str[0]
    # create a gene column by splitting Target_ID at ':'
    exonerate_results_new.loc[:, 'Path'] = exonerate_results_new['Target_ID'].str.split(':').str[0]
    #reorder columns
    exonerate_results_new = exonerate_results_new[['Assembly', 'Path', 'Target_Length', 'Target_Start', 'Target_End', 'Target_Orientation', 'Query_ID', 'Query_Length', 'Query_Alignment_Begin', 'Query_Alignment_End', 'Alignment_Length', 'Percent_Identity', 'Coverage', 'Classification', 'Info']]
    #rename columns
    exonerate_results_new.columns = ['Assembly', 'Path', 'Path_length', 'Path_start', 'Path_end', 'Path_strand', 'Feature', 'Feature_length', 'Feature_start', 'Feature_end', 'Alignment_length', 'Estimated_identity', 'Coverage', 'Classification', 'Attributes']
    # keep only rows where Classification is 'gene'
    exonerate_new_genes = exonerate_results_new[exonerate_results_new['Classification'] == 'gene']
    # drop column 'Classification'
    exonerate_new_genes = exonerate_new_genes.drop(columns=['Classification'])
    # add column called 'Source' to the dataframe and set all values to 'exonerate_results_new' 
    exonerate_new_genes['Source'] = 'exonerate'
    # create a new df 'exonerate_new_genes_attributes' by selecting the columns 'Feature' and 'Attributes'
    exonerate_new_genes_attributes = exonerate_new_genes[['Feature', 'Attributes']]
    # drop column 'Attributes'
    exonerate_new_genes = exonerate_new_genes.drop(columns=['Attributes'])
    return exonerate_new_genes, exonerate_new_genes_attributes

def merge_pafs_and_count(untangle_paf, exonerate_tair_genes, exonerate_new_genes, number_of_assemblies, softcore_limit):
    logging.info('Merging PAF files and counting...')
    # merge the untangle_paf with the exonerate_tair_genes and exonerate_new_genes
    all_data = pd.concat([untangle_paf, exonerate_tair_genes, exonerate_new_genes])
    # Group by 'Feature' and 'Assembly', then count unique 'Assembly' for each 'Feature'
    counts_total = all_data.groupby(['Feature', 'Assembly']).size().reset_index().groupby('Feature').size()
    # Convert the Series to a DataFrame
    counts_total = counts_total.to_frame(name='Total_count').reset_index()
    # Count how many times a Feature is present in each Assembly, and store ID, Estimated_identity, Coverage and source
    counts_x_assembly = (all_data.groupby(['Feature', 'Assembly', 'Estimated_identity', 'Coverage', 'Source'])
            .size()
            .reset_index(name='count')
            .groupby('Feature')
            .apply(lambda x: list(zip(x['Assembly'], x['Estimated_identity'], x['Coverage'], x['Source'], x['count'])), include_groups=False)
            .reset_index(name='Assembly_count'))
    # Merge the two dataframes
    counts = pd.merge(counts_total, counts_x_assembly, on='Feature')
    # Create a new column 'Label' based on the count
    labels = []
    for i, row in counts.iterrows():
        count = row['Total_count']
        feature = row['Feature']
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
    counts['Label'] = labels
    # Add a new column 'Assembly_names' to store only the names of the assemblies
    counts['Assembly_names'] = counts['Assembly_count'].apply(lambda x: ', '.join(sorted(set([assembly[0] for assembly in x]))))
    # save the dataframe to a csv file
    counts.to_csv('counts_with_labels.csv', sep='\t', index=False)
    return counts

def merge_attributes(untangle_attributes, exonerate_tair_genes_attributes, exonerate_new_genes_attributes):
    logging.info('Merging attributes...')
    # Merge the untangle_attributes with the exonerate_tair_genes_attributes and exonerate_new_genes_attributes
    all_attributes = pd.concat([untangle_attributes, exonerate_tair_genes_attributes, exonerate_new_genes_attributes])
    
    # Group by 'Feature' and merge 'Attributes'
    all_attributes['Attributes'] = all_attributes['Attributes'].fillna('')
    all_attributes = all_attributes.groupby('Feature')['Attributes'].apply(lambda x: ', '.join(set(', '.join(x).split(', ')))).reset_index()

    # Save the DataFrame to a CSV file
    all_attributes.to_csv('all_genes_attributes.csv', sep='\t', index=False)
    
    return all_attributes

def merge_paf_and_attributes(counts, all_attributes):
    logging.info('Merging PAF and attributes...')
    # Merge the all_data with the all_attributes
    all_data_final = pd.merge(counts, all_attributes, on='Feature')
    # Save the DataFrame to a CSV file
    all_data_final.to_csv('final_genes_screening.csv', sep='\t', index=False)

def get_statistics(counts):
    logging.info('Calculating statistics...')
    # Calculate statistics
    total_features = len(counts)
        # Number and percentage of core, softcore, dispensable and reference_private features
    core_features = len(counts[counts['Label'] == 'core'])
    core_percentage = (core_features / total_features) * 100

    softcore_features = len(counts[counts['Label'] == 'softcore'])
    softcore_percentage = (softcore_features / total_features) * 100

    dispensable_features = len(counts[counts['Label'] == 'dispensable'])
    dispensable_percentage = (dispensable_features / total_features) * 100

    private_features = len(counts[counts['Label'] == 'private'])
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
        'private_percentage': private_percentage
    }

    # Convert the statistics to a DataFrame
    stats_df = pd.DataFrame([stats])
    # Save the DataFrame to a CSV file
    stats_df.to_csv('statistics.csv', sep='\t', index=False)


def get_pav_matrix(counts):
    logging.info('Creating presence-absence matrix...')
    # Step 1: Extract unique assembly names
    all_assemblies = set()
    counts['Assembly_names'].str.split(', ').apply(all_assemblies.update)
    all_assemblies = sorted(all_assemblies)  # Sort the assembly names for consistent column ordering

    # Step 2: Initialize the matrix
    presence_absence_matrix = pd.DataFrame(0, index=counts['Feature'], columns=all_assemblies)

    # Step 3: Populate the matrix
    for index, row in counts.iterrows():
        assemblies = row['Assembly_names'].split(', ')
        presence_absence_matrix.loc[row['Feature'], assemblies] = 1

    # save the matrix to a CSV file
    presence_absence_matrix.to_csv('presence_absence_matrix.csv', sep='\t')

def main():
    parser = argparse.ArgumentParser(description='Process and merge data from multiple CSV files.')
    parser.add_argument('-u', '--untangle_paf_file', required=True, help='Path to the untangle_paf file. This must be "filter_passed_features.csv" coming from core_dispensable_genes.py.')
    parser.add_argument('-p', '--pangenome_screening_file', required=True, help='Path to the pangenome_screening file. This file must be "pangenome_screening.csv" coming from core_dispensable_genes.py')
    parser.add_argument('-t', '--exonerate_results_tair_file', required=True, help='Path to the exonerate_results_tair file. This file must be "corresponding_to_tair.csv" coming from parse_exonerate.py')
    parser.add_argument('-n', '--exonerate_results_new_file', required=True, help='Path to the exonerate_results_new file. This file must be "reviewed_not_tair_results.tsv" coming from parse_exonerate.py')
    parser.add_argument('-a', '--number_of_assemblies', type=int, required=True, help='Number of assemblies.')
    parser.add_argument('-s', '--softcore_limit', type=int, required=True, help='Softcore limit.')
    
    args = parser.parse_args()

    untangle_paf_data = untangle_paf(args.untangle_paf_file)
    untangle_attributes = untangle_paf_attributes(args.pangenome_screening_file)
    exonerate_tair_genes, exonerate_tair_genes_attributes = exonerate_results_tair(args.exonerate_results_tair_file)
    exonerate_new_genes, exonerate_new_genes_attributes = exonerate_results_new(args.exonerate_results_new_file)
    
    counts = merge_pafs_and_count(untangle_paf_data, exonerate_tair_genes, exonerate_new_genes, args.number_of_assemblies, args.softcore_limit)
    all_attributes = merge_attributes(untangle_attributes, exonerate_tair_genes_attributes, exonerate_new_genes_attributes)
    all_data_final = merge_paf_and_attributes(counts, all_attributes)
    get_statistics(counts)
    get_pav_matrix(counts)

    logging.info('All done!')

if __name__ == "__main__":
    main()
