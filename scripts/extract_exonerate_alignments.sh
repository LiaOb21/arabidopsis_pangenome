#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 3 ]; then
  echo "Usage: $0 <exonerate_output_file> <identifiers_file> <output_file>"
  exit 1
fi

EXONERATE_FILE=$1
IDENTIFIERS_FILE=$2
OUTPUT_FILE=$3

# Read identifiers into an array
readarray -t IDENTIFIERS < "$IDENTIFIERS_FILE"

# Create an empty output file
> "$OUTPUT_FILE"

# Extract alignments for each identifier
for ID in "${IDENTIFIERS[@]}"; do
  awk -v id="$ID" '
  BEGIN { RS="-- completed exonerate analysis" }
  $0 ~ id { print $0 "\n-- completed exonerate analysis" }
  ' "$EXONERATE_FILE" >> "$OUTPUT_FILE"
done

