#!/bin/bash

# Input and Output file names
input_file="pav_data_set_R.txt"
output_file="pav_collapsed_output.txt"

# First, extract all unique GCA codes
awk '{ print $1 }' "$input_file" | sort -u | while read gca; do
  # Initialize empty variables for each category
  core="NA"
  softcore="NA"
  dispensable="NA"
  private="NA"
  
  # Now, find the values for each category
  while read line; do
    gca_code=$(echo "$line" | awk '{print $1}')
    category=$(echo "$line" | awk '{print $2}')
    value=$(echo "$line" | awk '{print $3}')
    
    # Assign the value to the corresponding variable based on the category
    if [ "$gca_code" == "$gca" ]; then
      case $category in
        core)
          core=$value
          ;;
        softcore)
          softcore=$value
          ;;
        dispensable)
          dispensable=$value
          ;;
        private)
          private=$value
          ;;
      esac
    fi
  done < "$input_file"

  # Print the row for this GCA code
  echo -e "$gca\t$core\t$softcore\t$dispensable\t$private"
done > "$output_file"

# Add the header to the output file
sed -i '1i GCA\tcore\tsoftcore\tdispensable\tprivate' "$output_file"
