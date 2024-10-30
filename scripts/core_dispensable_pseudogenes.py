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
    parser.add_argument('-chr', '--reference_chromosome', type=str, required=True, help='The reference chromosome name corresponding to the analysed community in the GFF file. E.g. CP002684.1')
    parser.add_argument('-ref', '--reference', type=str, required=True, help='The reference genome name in the GFF file (must be the same as in the graph, i.e what is coming before the first # in the path name, e.g. GCA_000001735)')

    # Parse the arguments
    return parser.parse_args()


def save_df_to_csv(df, filename):
    """Save a DataFrame to a CSV file."""
    df.to_csv(filename, sep='\t', index=False)


def process_paf(input_file, id_threshold, jc_threshold, coverage):
    """Process the PAF (collapsed_and_merged.paf)file and filter the features based on the IDF, JCF, and alignment length ratio. Save the filtered features to a CSV file and return the filtered DataFrame. Save the filtered out features to a CSV file."""
    # Read pre-processed paf file with collapse_paf.py and merge_paf.py
    df = pd.read_csv(input_file, sep='\t')

    # Add assembly_name column: split the path_name at the '#' and keep the first part
    df["Assembly"] = df["Path"].str.split("#").str[0]

    # Extract coordinates from the 'Feature' column and add them as separate columns
    # Feature looks like: GCA_000001735#1#CP002684.1:10003207-10003479
    df[['path_dup', 'Start', 'End']] = df['Feature'].str.extract(r'(.*):(.*)-(.*)')

    #Drop the 'path_dup' column (we already have this information in the 'Path' column)
    df.drop(columns=['path_dup'], inplace=True)

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
    save_df_to_csv(filtered_out_df, 'filter_id_excluded_features.csv')
    # Keep only the features that pass the filter
    filtered_df = filtered_df[(filtered_df['Estimated_identity'] >= id_threshold) & (filtered_df['Jaccard_index'] >= jc_threshold) & (filtered_df['Coverage'] >= coverage)]

    # Save the DataFrame to a CSV file
    save_df_to_csv(filtered_df, 'filter_id_passed_features.csv')

    return filtered_df


def load_and_process_gff(gff_file, feature_type, reference_chromosome):
    """Load the GFF file into a DataFrame and filter it based on the feature type. Save the filtered DataFrame to a CSV file and return it."""
    # Load the GFF file into a DataFrame
    gff_df = pd.read_csv(gff_file, sep='\t', header=None, comment='#')
    # Rename the columns you're interested in
    gff_df.columns = ['Chromosome', 'Source', 'Type', 'Start', 'End', 'Score', 'Strand', 'Phase', 'Attributes']
    # Keep only rows where 'type' is feature_type
    gff_df = gff_df[gff_df['Type'] == feature_type]
    # Extract Pater and Stop_codon columns from the Attributes column
    # attributes look like:
    # coverage=0.09;Genewise_score=56.39;position=UTR;size=129;Type=FRAG;Pater=AT1G01010.1;Species=A.thaliana;Stop_codon=0;Frameshift=0
    gff_df['Pater'] = gff_df['Attributes'].str.extract(r'Pater=(.*?);')
    gff_df['Stop_codon'] = gff_df['Attributes'].str.extract(r'Stop_codon=(.*?);')
    # remove what comes after the first dot in the Pater column
    gff_df['Pater'] = gff_df['Pater'].str.split('.').str[0]
    # From gff_df, keep only the rows where chromosome = reference_chromosome. This is a further filter to make sure that the pseudogenes 
    # are on the same chromosome as the alignments, because we are merging only based on coordinates as pseudogenes don't have a real name 
    # but only a pater name
    gff_df_chr = gff_df[gff_df['Chromosome'] == reference_chromosome]
    #save the filtered gff_df_chr to a new file
    gff_df_chr.to_csv('chr.gff', index=False, sep='\t')
    return gff_df_chr


def merge_gff(filtered_df, gff_df_chr):
    # Merge the two DataFrames: filtered_df using 'Start' and 'End', and gff_df using 'Start' and 'End'
    # make sure that they are the same type first
    filtered_df['Start'] = filtered_df['Start'].astype(int)
    filtered_df['End'] = filtered_df['End'].astype(int)
    # Correctly setting values using .loc to avoid SettingWithCopyWarning
    gff_df_chr.loc[:, 'Start'] = gff_df_chr['Start'].astype(int)
    gff_df_chr.loc[:, 'End'] = gff_df_chr['End'].astype(int)
    # Merge the two DataFrames: filtered_df using 'Start' and 'End', and gff_df using 'Start' and 'End'
    merged_df = pd.merge(filtered_df, gff_df_chr, on=['Start', 'End'], how='left')
    # Save the merged DataFrame to a new file
    merged_df.to_csv('filtered_id_with_ambiguous.csv', index=False, sep='\t')
    return merged_df

def remove_ambiguous_features(merged_df, reference):
    # exclude rows where Assembly is not the reference, Stop_codons = 0, and coverage = 100. 
    # These three conditions must be met at the same time
    merged_df['Stop_codon'] = merged_df['Stop_codon'].astype(int)
    # Define the conditions
    condition1 = merged_df['Assembly'].str.contains(reference) # Keep all the rows where Assembly is the reference
    condition2 = merged_df['Stop_codon'] > 0  # Keep always if there's a stop codon, regardless of other conditions
    condition3 = (~merged_df['Assembly'].str.contains(reference)) & (merged_df['Coverage'] < 100) & (merged_df['Stop_codon'] == 0) # When Assembly is not the reference, rows with Stop_codon = 0 are kept only if Coverage < 100
    # Combine conditions
    final_condition = condition1 | condition2 | condition3
    # Filter the DataFrame
    merged_df_no_ambiguous = merged_df[final_condition]
    # Save the filtered DataFrame to a new file
    merged_df_no_ambiguous.to_csv('filtered_id_no_ambiguous.csv', index=False, sep='\t')
    return merged_df_no_ambiguous

def label_genes(merged_df_no_ambiguous, number_of_assemblies, softcore_limit):
    """Label the features as core, softcore, or dispensable based on the number of assemblies they are present in. Save the result to a CSV file and return the DataFrame."""
    # Group by 'Pater' and 'Assembly', then count unique 'Assembly' for each 'Pater'
    counts_total = merged_df_no_ambiguous.groupby(['Pater', 'Assembly']).size().reset_index().groupby('Pater').size()
    # Convert the Series to a DataFrame
    counts_total = counts_total.to_frame(name='Total_count').reset_index()

    # Count how many times a Pater is present in each Assembly
    counts_x_assembly = (merged_df_no_ambiguous.groupby(['Pater', 'Assembly'])
          .size()
          .reset_index(name='count')
          .groupby('Pater')
          .apply(lambda x: list(zip(x['Assembly'], x['count'])), include_groups=False)
          .reset_index(name='Assembly_count'))

    # Merge the two DataFrames
    counts = pd.merge(counts_total, counts_x_assembly, on='Pater')

    # Save the DataFrame immediately after it's created to check eventual errors
    save_df_to_csv(counts, 'counts_before_labeling.csv')

    # Create a new column 'Label' based on the count
    labels = []
    for i, row in counts.iterrows():
        count = row['Total_count']
        pater = row['Pater']
        if count > number_of_assemblies:
            raise ValueError(f"Count greater than {number_of_assemblies} for pater '{pater}'")
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

def get_final_files(merged_df_no_ambiguous, gff_df_chr, counts):
    """This is to keep compatibility with the downstream scripts after changing the classification strategy"""
    # From merged_df_no_ambiguous, keep columns for filter_passed_features.csv, keeping Pater rather than Feature
    final_paf = merged_df_no_ambiguous[['Assembly', 'Path', 'Path_length', 'Path_start', 'Path_end', 'Path_strand', 'Pater', 'Feature_length', 'Feature_start', 'Feature_end', 'Alignment_length', 'Estimated_identity', 'Jaccard_index', 'Self_coverage', 'Coverage']]
    # Save the final paf to a CSV file
    save_df_to_csv(final_paf, 'filter_passed_features.csv')
    # From gff_df_chr keep columns to be merged to counts_with_labels.csv to obtain pangenome_screening.csv
    columns_to_merge = gff_df_chr[['Chromosome', 'Pater', 'Start', 'End', 'Type', 'Strand', 'Attributes']]
    # Before merging these columns with counts, we have to group for unique pater names
    columns_to_merge_grouped = columns_to_merge.groupby('Pater').apply(
        lambda x: pd.Series({
            'Chromosome': ', '.join(sorted(set(x['Chromosome']))),
            'Start': ', '.join(sorted(set(x['Start'].astype(str)))),
            'End': ', '.join(sorted(set(x['End'].astype(str)))),
            'Type': ', '.join(sorted(set(x['Type']))),
            'Strand': ', '.join(sorted(set(x['Strand'].astype(str)))),
            'Attributes': ', '.join(sorted(set(x['Attributes'])))
        })
    ).reset_index()
    # Merge the two DataFrames
    merged_counts_gff = pd.merge(counts, columns_to_merge_grouped, on=['Pater'], how='left')
    # order columns
    merged_counts_gff = merged_counts_gff[['Chromosome', 'Pater', 'Type', 'Start', 'End', 'Strand', 'Total_count', 'Label', 'Attributes', 'Assembly_count']]
    # Save the DataFrame to a CSV file
    save_df_to_csv(merged_counts_gff, 'pangenome_screening.csv')
    return columns_to_merge_grouped, merged_counts_gff

def unmerged_features(counts, columns_to_merge_grouped):
    """Find the features in the GFF file that were not merged with counts_df."""
    # Merge counts with the selected columns from gff_df grouped by 'Pater'
    unmerged_features_df = pd.merge(counts, columns_to_merge_grouped, 
                        on=['Pater'], 
                        how='outer')

    # Filter rows where 'Total_count' is NaN 
    unmerged = unmerged_features_df[unmerged_features_df['Total_count'].isna()]

    # Save unmerged features to a CSV file
    save_df_to_csv(unmerged, 'not_annotated_features.csv')
    return unmerged


def get_statistics(merged_counts_gff, counts, gff_df_chr, unmerged):
    """Calculate statistics and return them."""
    # Calculate statistics
    # Number of features and percentage of annotated features in the pangenome compared to the total number of features in the GFF file (based on Pater)
    #extract unique Pater values from gff_df_chr
    filtered_gff_df = gff_df_chr.drop_duplicates(subset=['Pater'])
    total_features = len(filtered_gff_df)
    annotated_features = len(merged_counts_gff)
    annotated_percentage = (annotated_features / total_features) * 100

    # Number and percentage of non-annotated features in the pangenome compared to the total number of features in the GFF file
    non_annotated_features = len(unmerged)
    non_annotated_percentage = (non_annotated_features / total_features) * 100

    # Check if non_annotated_features is different from total_features - annotated_features
    if non_annotated_features != total_features - annotated_features:
        raise ValueError("Non-annotated features is not equal to total features - annotated features")

    # Number and percentage of core, softcore, dispensable and private features calculated on the annotated features
    core_features = len(counts[counts['Label'] == 'core'])
    core_percentage = (core_features / annotated_features) * 100

    softcore_features = len(counts[counts['Label'] == 'softcore'])
    softcore_percentage = (softcore_features /annotated_features) * 100

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
    stats_df.to_csv('statistics.csv', index=False, sep='\t')

def get_matrix(merged_counts_gff):
    """Generate a matrix from the merged_counts_gff (i.e. pangenome_screening.csv) DataFrame and save it to a CSV file."""
    # Step 1: Extract 'Pater' and 'Assembly_count' columns
    matrix_df = merged_counts_gff[['Pater', 'Assembly_count']]
    # Step 2: Explode the 'Assembly_Count' column
    data_exploded = matrix_df.explode('Assembly_count')
    # Step 3: Split the tuples into two columns: 'Assembly' and 'Count'
    data_exploded[['Assembly', 'Count']] = pd.DataFrame(data_exploded['Assembly_count'].tolist(), index=data_exploded.index)
    # Step 4: Pivot the table
    matrix = data_exploded.pivot_table(index='Pater', columns='Assembly', values='Count', fill_value=0)
    # save
    matrix.to_csv('matrix.csv', index=True, sep='\t')
    # The following is for debugging purposes
    # Step 5: Find the maximum value
    max_value = matrix.max().max()
    # Step 6: Locate the position of the maximum value
    max_position = matrix.stack().idxmax()
    # Step 7: Print the details
    print(f"Maximum value: {max_value} found for Pater: {max_position[0]} and Assembly: {max_position[1]}")

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

    logging.info("Loading and processing GFF file...")
    gff_df_chr = load_and_process_gff(args.gff_file, args.feature_type, args.reference_chromosome)

    logging.info("Merging PAF and GFF files...")
    merged_df = merge_gff(filtered_df, gff_df_chr)

    logging.info("Removing ambiguous features...")
    merged_df_no_ambiguous = remove_ambiguous_features(merged_df, args.reference)
    
    logging.info("Labeling genes...")
    counts = label_genes(merged_df_no_ambiguous, args.number_of_assemblies, args.softcore_limit)

    logging.info("Getting final files...")
    columns_to_merge_grouped, merged_counts_gff = get_final_files(merged_df_no_ambiguous, gff_df_chr, counts)

    logging.info("Finding not annotated features...")
    unmerged = unmerged_features(counts, columns_to_merge_grouped)

    logging.info("Calculating statistics...")
    get_statistics(merged_counts_gff, counts, gff_df_chr, unmerged)

    logging.info("Generating matrix...")
    get_matrix(merged_counts_gff)

    logging.info("Script finished successfully. Enjoy your pangenome! 🎉🎉🎉")

if __name__ == "__main__":
    main()