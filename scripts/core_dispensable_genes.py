#!/usr/bin/env python3

import os
import pandas as pd
import argparse
import logging

# Set up logging
logging.basicConfig(level=logging.INFO)


def check_file_exists(filename):
    """Check if a file exists."""
    if not os.path.exists(filename):
        logging.error(f"File {filename} does not exist.")
        return False
    return True


def parse_args():
    """Parse the command line arguments."""
    # Create the parser
    parser = argparse.ArgumentParser(description='Find core, softcore, and dispensable genes in a pangenome.')

    # Add the arguments
    parser.add_argument('-i', '--input_file', type=str, required=True, help='The input PAF file (output of odgi untangle)')
    parser.add_argument('-id', '--id_threshold', type=float, required=True, help='The ID (estimated identity) threshold for filtering the PAF file (0 <= idf <= 100). For example, if you want to keep only alignments with IDF >= 95, you should use 0.95 as the value for this argument.')
    parser.add_argument('-jc', '--jc_threshold', type=float, required=True, help='The JC (Jaccard index) threshold for filtering the PAF file (0 <= jcf <= 1). For example, if you want to keep only alignments with JCF >= 0.95, you should use 0.95 as the value for this argument.')
    parser.add_argument('-cov', '--coverage', type=int, required=True, help='The coverage for filtering the PAF file (coverage = Alignment_length / Feature_length * 100). For example, if you want to keep only alignments with coverage >= 95%%, you should use 95 as value for this argument.')
    parser.add_argument('-na', '--number_of_assemblies', type=int, required=True, help='The number of assemblies in the pangenome')
    parser.add_argument('-sl', '--softcore_limit', type=int, required=True, help='Number of assemblies above which a gene is considered softcore (sl< softcore < na). For example, if you have 10 assemblies in your pangenome and you want to consider genes that are present in 8-9 assemblies as softcore, you should use 8 as the value for this argument.')
    parser.add_argument('-gff', '--gff_file', type=str, required=True, help='The GFF used to untangle the features across the assemblies of the pangenome')
    parser.add_argument('-ft', '--feature_type', type=str, required=True, help='The feature type to keep in the GFF file (e.g. gene, pseudogene, etc.)')

    # Parse the arguments
    return parser.parse_args()


def save_df_to_csv(df, filename):
    """Save a DataFrame to a CSV file."""
    df.to_csv(filename, sep='\t', index=False)


def process_paf(input_file, id_threshold, jc_threshold, coverage):
    """Process the PAF (collapsed_and_merged.paf) file and filter the features based on the IDF, JCF, and coverage. Save the filtered features to a CSV file and return the filtered DataFrame. Save the filtered out features to a CSV file."""
    # Read pre-processed paf file with collapse_paf.py and merge_paf.py
    df = pd.read_csv(input_file, sep='\t')

    # Add assembly_name column: split the path_name at the '#' and keep the first part
    df["Assembly"] = df["Path"].str.split("#").str[0]

    # Move the Assembly column to the first position
    cols = df.columns.tolist()
    cols.insert(0, cols.pop(cols.index('Assembly')))
    df = df[cols]

    # Calculate gene coverage
    df["Coverage"] = ((df["Feature_end"] - df["Feature_start"]) / df["Feature_length"]) * 100

    # Remove rows with coverage > 100 and save them to a new file. These rows are ambiguous
    high_coverage = df[df['Coverage'] > 100]
    high_coverage.to_csv('excluded_for_high_coverage.paf', index=False, sep='\t')

    # Exclude rows where Coverage is higher than 100 from the DataFrame
    filtered_df = df[~(df['Coverage'] > 100)]

    # Filter DataFrame based on the conditions: this allows to remove alignments with low estimated identity and low coverage
    # First save features that don't pass the filter to a new file
    filtered_out_df = filtered_df[(filtered_df['Estimated_identity'] < id_threshold) | (filtered_df['Jaccard_index'] < jc_threshold) | (filtered_df['Coverage'] < coverage)]
    # Save the filtered out features to a CSV file
    save_df_to_csv(filtered_out_df, 'filter_excluded_features.csv')
    # Keep only the features that pass the filter
    filtered_df = filtered_df[(filtered_df['Estimated_identity'] >= id_threshold) & (filtered_df['Jaccard_index'] >= jc_threshold) & (filtered_df['Coverage'] >= coverage)]

    # Save the DataFrame to a CSV file
    save_df_to_csv(filtered_df, 'filter_passed_features.csv')

    return filtered_df


def label_genes(filtered_df, number_of_assemblies, softcore_limit):
    """Label the features as core, softcore, or dispensable based on the number of assemblies they are present in. Save the result to a CSV file and return the DataFrame."""
    # Group by 'Feature' and 'Assembly', then count unique 'Assembly' for each 'Feature'
    counts_total = filtered_df.groupby(['Feature', 'Assembly']).size().reset_index().groupby('Feature').size()
    # Convert the Series to a DataFrame
    counts_total = counts_total.to_frame(name='Total_count').reset_index()

    # Count how many times a Feature is present in each Assembly
    counts_x_assembly = (filtered_df.groupby(['Feature', 'Assembly'])
          .size()
          .reset_index(name='count')
          .groupby('Feature')
          .apply(lambda x: list(zip(x['Assembly'], x['count'])), include_groups=False)
          .reset_index(name='Assembly_count'))

    # Merge the two DataFrames
    counts = pd.merge(counts_total, counts_x_assembly, on='Feature')

    # Save the DataFrame immediately after it's created to check eventual errors
    save_df_to_csv(counts, 'counts_before_labeling.csv')

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

    # Save the result to a CSV file
    save_df_to_csv(counts, 'counts_with_labels.csv')

    return counts


def merge_gff(counts, gff_file, feature_type):
    """Merge the counts_df with the selected columns from the GFF file and save the result to a CSV file."""

    # Extract 'Feature_start' and 'Feature_end' using a regular expression from the 'Feature' column of the counts DataFrame
    counts[['Feature', 'Feature_start', 'Feature_end']] = counts.iloc[:, 0].str.extract(r'(.*):(.*)-(.*)')

    # change columns order
    counts = counts[['Feature', 'Feature_start', 'Feature_end', 'Total_count', 'Label', 'Assembly_count']]

    # Load the GFF file into a DataFrame
    gff_df = pd.read_csv(gff_file, sep='\t', header=None, comment='#')

    # Rename the columns you're interested in
    gff_df.columns = ['Chromosome', 'Source', 'Type', 'Start', 'End', 'Score', 'Strand', 'Phase', 'Attributes']

    # Keep only rows where 'type' is feature_type
    gff_df = gff_df[gff_df['Type'] == feature_type]

    # Extract 'ID' and 'Name' from the 'attributes' column
    gff_df['Feature_name'] = gff_df['Attributes'].str.extract(r'ID=gene-([^;]*)')

    # Make sure that the fields used for merging are of the same type
    counts.loc[:, 'Feature_start'] = counts['Feature_start'].astype(int)
    counts.loc[:, 'Feature_end'] = counts['Feature_end'].astype(int)

    gff_df.loc[:, 'Start'] = gff_df['Start'].astype(int)
    gff_df.loc[:, 'End'] = gff_df['End'].astype(int)

    # Merge counts_df with the selected columns from gff_df
    merged_counts_gff = pd.merge(counts, gff_df, 
                        left_on=['Feature', 'Feature_start', 'Feature_end'], 
                        right_on=['Feature_name', 'Start', 'End'],
                        how='left')

    # Drop unnecessary columns and the redundant 'start' and 'end' columns
    merged_counts_gff = merged_counts_gff.drop(columns=['Feature_name', 'Source', 'Start', 'End', 'Score', 'Phase'])

    # sort merged_df_gff by 'Feature_start'
    merged_counts_gff = merged_counts_gff.sort_values(by=['Feature_start'])

    # order columns
    merged_counts_gff = merged_counts_gff[['Chromosome', 'Feature', 'Type', 'Feature_start', 'Feature_end', 'Strand', 'Total_count', 'Label', 'Attributes', 'Assembly_count']]
    
    # Save the DataFrame to a CSV file
    save_df_to_csv(merged_counts_gff, 'pangenome_screening.csv') 

    # Extract the 'chromosome' value from the first row
    chromosome_name = merged_counts_gff['Chromosome'].iloc[0]

    # Check if all 'chromosome' values are the same
    if not merged_counts_gff['Chromosome'].nunique() == 1:
        # Find the index of the first row where 'chromosome' is different
        diff_chromosome_index = merged_counts_gff[merged_counts_gff['Chromosome'] != chromosome_name].index[0]
        raise ValueError(f"Error: 'Chromosome' value is different at line {diff_chromosome_index + 1} in 'pangenome_screening.csv'")

    return counts, gff_df, merged_counts_gff, chromosome_name

def unmerged_features(counts, gff_df, chromosome_name):
    """Find the features in the GFF file that were not merged with counts_df."""
    # Filter gff_df to keep only the 'chromosome' values that are present in merged_df_gff
    filtered_gff_df = gff_df[gff_df['Chromosome'] == chromosome_name]

    # Merge counts_df with the selected columns from gff_df
    unmerged_features_df = pd.merge(counts, filtered_gff_df, 
                        left_on=['Feature', 'Feature_start', 'Feature_end'], 
                        right_on=['Feature_name', 'Start', 'End'],
                        how='outer', indicator=True)

    # Create a DataFrame of rows from gff_df that didn't merge with counts_df
    unmerged_gff = unmerged_features_df[unmerged_features_df['_merge'] == 'right_only']

    # Drop the columns that came from filtered_gff_df
    unmerged_gff = unmerged_gff.dropna(axis=1)

    # Drop the merge indicator column
    unmerged_gff = unmerged_gff.drop(columns=['_merge'])

    # Save unmerged features to a CSV file
    save_df_to_csv(unmerged_gff, 'not_annotated_features.csv')

    return filtered_gff_df, unmerged_gff


def get_statistics(merged_counts_gff, counts, filtered_gff_df, unmerged_gff):
    """Calculate statistics and return them."""
    # Calculate statistics
    # Number of features and percentage of annotated features in the pangenome compared to the total number of features in the GFF file
    total_features = len(filtered_gff_df)
    annotated_features = len(merged_counts_gff)
    annotated_percentage = (annotated_features / total_features) * 100

    # Number and percentage of non-annotated features in the pangenome compared to the total number of features in the GFF file
    non_annotated_features = len(unmerged_gff)
    non_annotated_percentage = (non_annotated_features / total_features) * 100

    # Check if non_annotated_features is different from total_features - annotated_features
    if non_annotated_features != total_features - annotated_features:
        raise ValueError("Non-annotated features is not equal to total features - annotated features")

    # Number and percentage of core, softcore, dispensable and reference_private features calculated on annotated features
    core_features = len(counts[counts['Label'] == 'core'])
    core_percentage = (core_features / annotated_features) * 100

    softcore_features = len(counts[counts['Label'] == 'softcore'])
    softcore_percentage = (softcore_features / annotated_features) * 100

    dispensable_features = len(counts[counts['Label'] == 'dispensable'])
    dispensable_percentage = (dispensable_features / annotated_features) * 100

    private_features = len(counts[counts['Label'] == 'private'])
    private_percentage = (private_features / annotated_features) * 100

    # Create a dictionary to store the statistics
    stats = {
        'total_features': total_features,
        'annotated_features': annotated_features,
        'annotated_percentage': annotated_percentage,
        'core_features': core_features,
        'core_percentage': core_percentage,
        'softcore_features': softcore_features,
        'softcore_percentage': softcore_percentage,
        'dispensable_features': dispensable_features,
        'dispensable_percentage': dispensable_percentage,
        'private_features': private_features,
        'private_percentage': private_percentage,
        'non_annotated_features': non_annotated_features,
        'non_annotated_percentage': non_annotated_percentage
    }

    # Convert the statistics to a DataFrame
    stats_df = pd.DataFrame([stats])

    # Save the statistics to a CSV file
    save_df_to_csv(stats_df, 'statistics.csv')

def get_matrix(merged_counts_gff):
    """Generate a matrix from the merged_counts_gff (i.e. pangenome_screening.csv) DataFrame and save it to a CSV file."""
    # Step 1: Extract 'Feature' and 'Assembly_count' columns
    matrix_df = merged_counts_gff[['Feature', 'Assembly_count']]
    # Step 2: Explode the 'Assembly_Count' column
    data_exploded = matrix_df.explode('Assembly_count')
    # Step 3: Split the tuples into two columns: 'Assembly' and 'Count'
    data_exploded[['Assembly', 'Count']] = pd.DataFrame(data_exploded['Assembly_count'].tolist(), index=data_exploded.index)
    # Step 4: Pivot the table
    matrix = data_exploded.pivot_table(index='Feature', columns='Assembly', values='Count', fill_value=0)
    # save
    matrix.to_csv('matrix.csv', index=True, sep='\t')
    # The following is for debugging purposes
    # Step 5: Find the maximum value
    max_value = matrix.max().max()
    # Step 6: Locate the position of the maximum value
    max_position = matrix.stack().idxmax()
    # Step 7: Print the details
    print(f"Maximum value: {max_value} found for Feature: {max_position[0]} and Assembly: {max_position[1]}")

def main():
    """Main function to run the script."""
    # Parse the arguments
    args = parse_args()

    # Check if the input file exists
    if not check_file_exists(args.input_file) or not check_file_exists(args.gff_file):
        return
    
    # Call the functions
    logging.info("Processing PAF file...")
    filtered_df = process_paf(args.input_file, args.id_threshold, args.jc_threshold, args.coverage)
    
    logging.info("Labeling genes...")
    counts = label_genes(filtered_df, args.number_of_assemblies, args.softcore_limit)

    logging.info("Merging with GFF file...")
    counts, gff_df, merged_counts_gff, chromosome_name = merge_gff(counts, args.gff_file, args.feature_type)

    logging.info("Finding not annotated features...")
    filtered_gff_df, unmerged_gff = unmerged_features(counts, gff_df, chromosome_name)

    logging.info("Calculating statistics...")
    get_statistics(merged_counts_gff, counts, filtered_gff_df, unmerged_gff)

    logging.info("Generating matrix...")
    get_matrix(merged_counts_gff)

    logging.info("Script finished successfully. Enjoy your pangenome! 🎉🎉🎉")

if __name__ == "__main__":
    main()