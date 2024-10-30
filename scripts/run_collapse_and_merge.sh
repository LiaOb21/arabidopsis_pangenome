#!/bin/bash

# Initialize variables
input_file=""
collapse_script_path=""
merge_script_path=""
distance=""
job_limit=""

# Parse command-line options
while getopts i:c:m:d:j: flag; do
    case "${flag}" in
        i) input_file="${OPTARG}";;
        c) collapse_script_path="${OPTARG}";;
        m) merge_script_path="${OPTARG}";;
        d) distance="${OPTARG}";;
        j) job_limit="${OPTARG}";;
        ?) echo "Invalid option: -$OPTARG"
           exit 1;;
        *) echo "Usage: $0 -i <input_paf_file> -c <collapse_paf.py_script_path> -m <merge_paf.py_script_path> -d <max_distance_for_merging> -j <job_limit>"
           exit 1;;
    esac
done

# Check if all options were provided
if [ -z "$input_file" ] || [ -z "$collapse_script_path" ] || [ -z "$merge_script_path" ] || [ -z "$distance" ] || [ -z "$job_limit" ]; then
    echo "Usage: $0 -i <input_paf_file> -c <collapse_paf.py_script_path> -m <merge_paf.py_script_path> -d <max_distance_for_merging> -j <job_limit>"
    exit 1
fi

echo "Input file: $input_file"
echo "Collapse script path: $collapse_script_path"
echo "Merge script path: $merge_script_path"
echo "Max distance for merging: $distance"
echo "Job limit: $job_limit"

echo "Extracting assembly names from the input file..."

awk '$1 != $6' $input_file | cut -f1 | awk -F'[#]' '{print $1}' | sort | uniq > assembly_names.txt

echo "Splitting the input file into chunks based on the assembly names..."

while read assembly; do
    grep "^$assembly" $input_file > "${assembly}_chunk.paf" &
    if [ $? -ne 0 ]; then
        echo "Error: grep failed for assembly $assembly"
        exit 1
    fi
    # Limit the number of background jobs to avoid overloading the system
    while [ $(jobs | wc -l) -ge 4 ]; do
        sleep 1
    done
done < assembly_names.txt
wait

echo "Processing the chunks..."

# Loop over each chunk
for chunk in *_chunk.paf; do
    # Generate the output file names
    collapse_output="${chunk%.paf}_collapsed.paf"
    merge_output="${chunk%.paf}_merged.paf"

    # Run the collapse script and then the merge script in the background
    python3 "$collapse_script_path" -i "$chunk" -o "$collapse_output" && python3 "$merge_script_path" -i "$collapse_output" -o "$merge_output" -d "$distance" &
    if [ $? -ne 0 ]; then
        echo "Error: Python script failed for chunk $chunk"
        exit 1
    fi

    # Limit the number of background jobs to avoid overloading the system
    while [ $(jobs | wc -l) -ge "$job_limit" ]; do
        sleep 1
    done
done
wait

echo "All chunks have been processed successfully."
echo "Joining the merged files..."

# concatenate all the merged files
awk 'FNR > 1 || NR == 1' *_merged.paf > collapsed_and_merged.paf

echo "Cleaning up the intermediate files..."

# Cleanup
mkdir -p temp
mv *_chunk.paf *_collapsed.paf *_merged.paf temp/

echo "All done! The final merged output file is 'collapsed_and_merged.paf'."