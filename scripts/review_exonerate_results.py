#!/usr/bin/env python3

import pandas as pd
import argparse
import logging

# Set up the logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def get_uniprot_results(id_mapping_file):
    logging.info('Getting UniProt proteins IDs to keep...')
    # upload file of UniRef100 to UniProtKB/Swiss-Prot mapping containing reviewed proteins
    uniprot = pd.read_csv(id_mapping_file, sep='\t', header=0)
    # Keep unique values in 'From' column, i.e. UniRef100 IDs which correspond to reviewed proteins
    uniprot_reviewed_ids = uniprot['From'].unique()
    return uniprot_reviewed_ids


def review_results(not_ref_results, uniprot_reviewed_ids, output_file):
    logging.info('Filtering results...')
    # upload new_results_not_tair.tsv (this file is coming from parse_exonerate.py)
    not_ref_results = pd.read_csv(not_ref_results, sep='\t')
    # exclude rows from not_ref_results if 'Query_ID' doesn't have a correspondent in 'From' column from uniprot_reviewed_ids
    not_ref_results_filtered = not_ref_results[not_ref_results['Query_ID'].isin(uniprot_reviewed_ids)]
    # save the DataFrame in a file
    not_ref_results_filtered.to_csv(output_file, sep='\t', index=False)


def main():
    # Set up the argument parser
    parser = argparse.ArgumentParser(description='Review results excluding proteins that were not reviewed in UniProt')
    parser.add_argument('-i', '--id_mapping_file', type=str, required=True, help='UniProtKB/Swiss-Prot mapping file')
    parser.add_argument('-a', '--not_ref_results', type=str, required=True, help='This must be the file new_results_not_tair.tsv (this file is coming from parse_exonerate.py)')
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Output file name for the reviewed results')
    args = parser.parse_args()

    uniprot_reviewed_ids = get_uniprot_results(args.id_mapping_file)

    # Review the results
    review_results(args.not_ref_results, uniprot_reviewed_ids, args.output_file)

if __name__ == '__main__':
    main()