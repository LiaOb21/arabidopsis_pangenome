#!/bin/bash

# Function to display help message
show_help() {
    echo "Usage: $0 -i INPUT_FILE -l LINES_PER_CHUNK -s SEQ_FASTA -p PROTEIN_FASTA"
    echo "  -i INPUT_FILE         Path to the diamond results file"
    echo "  -l LINES_PER_CHUNK    Number of lines per chunk"
    echo "  -s SEQ_FASTA          Path to the sequence FASTA file"
    echo "  -p PROTEIN_FASTA      Path to the protein FASTA file"
}

# Parse command-line options
while getopts "hi:l:s:p:" opt; do
    case ${opt} in
        h )
            show_help
            exit 0
            ;;
        i )
            input_file=$OPTARG
            ;;
        l )
            lines_per_chunk=$OPTARG
            ;;
        s )
            seq_fasta=$OPTARG
            ;;
        p )
            protein_fasta=$OPTARG
            ;;
        \? )
            show_help
            exit 1
            ;;
    esac
done

# Check if all required options are provided
if [ -z "$input_file" ] || [ -z "$lines_per_chunk" ] || [ -z "$seq_fasta" ] || [ -z "$protein_fasta" ]; then
    show_help
    exit 1
fi

# Split the diamond results file into chunks
split -l "$lines_per_chunk" "$input_file" diamond_chunk_

# Loop over each chunk to extract IDs and corresponding FASTA sequences
for chunk in diamond_chunk_*
do
    # Extract unique sequence and protein IDs
    awk '{print $1}' $chunk | sort | uniq > ${chunk}_seq_ids.txt
    awk '{print $2}' $chunk | sort | uniq > ${chunk}_protein_ids.txt

    # Extract corresponding FASTA sequences
    seqtk subseq "$seq_fasta" ${chunk}_seq_ids.txt > ${chunk}_seqs.fasta
    seqtk subseq "$protein_fasta" ${chunk}_protein_ids.txt > ${chunk}_proteins.fasta

    # Optional: Remove intermediate ID files to save space
    rm ${chunk}_seq_ids.txt
    rm ${chunk}_protein_ids.txt
done

