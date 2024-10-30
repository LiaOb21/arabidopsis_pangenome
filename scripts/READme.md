Just a quick explanation of what each script in this directory does. 

# bed_to_gff3.ipynb

In this Jupyter notebook there is the explanation of how the original file of the pseudogenes annotations was transformed in GFF3 to be compatible with downstream analysis.

# collapse_paf.py

This Python script processes a PAF (Pairwise mApping Format) file, which is a format often used in bioinformatics to store sequence alignment data. The PAF file produced by `odgi untangle` is slightly different from the standard. This is the input file for this script. Here's a step-by-step breakdown:

1. import necessary libraries: pandas for data manipulation, numpy for numerical operations, argparse for handling command-line arguments, logging for log info.

2. The `load_and_preprocess` function reads a PAF file into a pandas DataFrame, drops unnecessary columns, renames the remaining columns, removes certain prefixes from some columns, and converts these columns to numeric types.

3. The `assign_block_ids` function assigns a unique ID to each block of data in the DataFrame. A new block ID is assigned whenever the `Path_start` value of the current row is not equal to the `Path_end` value of the previous row.

4. The collapse_paf function groups the DataFrame by `Path` and `Feature`, and processes each group (or block) separately. For each block, it sorts the data by `Path_start` and `Feature_start`, assigns block IDs, and then processes each sub-block separately. For each sub-block, it creates a new row with certain values (like the sum of `Alignment_length` and the weighted averages of `Estimated_identity`, `Jaccard_index`, and `Self_coverage`) and appends this row to a results DataFrame. Finally, it saves the results to the output file.

5. The script's main function sets up command-line argument parsing, calls `load_and_preprocess` to load and preprocess the input PAF file, and then calls `collapse_paf` to process the data and save the results to an output file.

# merge_paf.py

This Python script is designed to merge rows of a DataFrame based on certain conditions, i.e. if the gaps between the alignments are smaller than a distance given from command line. Here's a step-by-step breakdown of what the script does:

1. **Imports necessary modules**: The script imports pandas for data manipulation, argparse for command-line argument parsing, and logging for logging information.

2. **Sets up logging**: The logging level is set to INFO, and the format for the logging messages is defined.

3. **Defines the `merge_rows` function**: This function takes a DataFrame and a distance as input and returns a DataFrame where rows have been merged based on certain conditions.

    - The function first sorts the DataFrame by 'Path', 'Feature', 'Path_strand', 'Path_start', and 'Path_end'.
    - It then initializes an empty list to store the rows and variables to store the details of the current row.
    - The function then iterates over the rows of the DataFrame. If the current row can be merged with the previous row (based on the 'Path', 'Feature', 'Path_strand', and the difference between 'Path_start' and the current end), it updates the end of the current row and the alignment_length and the weighted averages. If the current row cannot be merged with the previous row, it calculates the final weighted averages, appends the current row to the list, and updates the current row and its details.
    - After iterating over all the rows, the function calculates the final weighted averages for the last row, appends the last row to the list, and converts the list of rows into a DataFrame.

4. **Defines the main part of the script**: it sets up command-line argument parsing, loads the input file into a DataFrame, merges the rows of the DataFrame using the `merge_rows` function, and saves the DataFrame to the output file.

The script is designed to be run from the command line with the input file, output file, and distance as command-line arguments. The input file should be a tab-separated file that can be read into a DataFrame. The output file will be a tab-separated file created from the DataFrame after the rows have been merged. The distance is used to determine whether rows can be merged.


# run_collapse_and_merge.sh

This bash script is designed to automatise the consequent running of `collapse_paf.py` and `merge_paf.py`. It allows the parallelisation of the jobs, which are executed in chunks, and rejoins together the results at the end of the process.

1. The script starts by initializing some variables and then parsing command-line options. The options specify the input file, the paths to the collapse and merge scripts, the maximum distance for merging, and the job limit.

2. It checks if all options were provided. If not, it prints a usage message and exits.

3. It prints the values of the options to the console.

4. It extracts the assembly names from the input file. It does this by selecting lines where the first and sixth fields are not equal (this is to avoid the inclusion of genes mapping against each other; for pseudogenes in our case a previous step is required to achieve the same thing - see Notebook.md), cutting out the first field, splitting it on the '#' character, sorting the results, and removing duplicates. The assembly names are saved to a file called `assembly_names.txt`.

5. It splits the input file into chunks based on the assembly names. For each assembly name, it selects lines from the input file that start with that assembly name and writes them to a new file. The new file's name is the assembly name followed by `_chunk.paf`. This is done in parallel, with a limit on the number of parallel jobs to avoid overloading the system.

6. It processes the chunks. For each chunk, it generates the output file names, runs the collapse script on the chunk, and then runs the merge script on the result. This is also done in parallel, with a limit on the number of parallel jobs.

7. It waits for all background jobs to finish.

8. It prints a message saying that all chunks have been processed successfully.

9. It concatenates all the merged files into a single output file, skipping the header line in all but the first file.

10. It cleans up the directory moving intermediate files to `temp` directory. All these files are kept for double checking of the results, but the entire directory can be removed if no anomalies are observed.

11. It prints a message saying that the final merged output file is `collapsed_and_merged.paf`.

In summary, this script splits an input file into chunks based on assembly names, processes each chunk in parallel using the collapse and merge scripts, and then concatenates the results into a single output file.


# core_dispensable_genes.py

This scripts analyses the output PAF from `odgi untangle` that was obtained for the gene annotations, after the processing with `run_collapse_and_merge.sh`.

The script apply some filters to the PAF file:
- coverage (which is calculated first) (`-cov`)
- estimated identity (`-id`)
- jaccard index (`-jc`)

The thresholds for these filters are decided by the user and given as input from command line. 
The alignments that do not meet the criteria decided by the user are saved in `filter_excluded_features.csv`. The alignments which are considered valid are saved in `filter_passed_features.csv`. The script also outputs a file called `excluded_for_high_coverage.paf` which contains annotations where the coverage is higher than 100%, which should be empty if everything worked fine. 

In order to identify core, softcore, dispensable and private genes the user should also input the total number of assemblies (`-na`) and the softcore limit (`-sl`). To classify genes into these categories they are grouped by Feature and Assembly. Depending on how many assembly contains a given gene (info saved in `Total_count` column), they are labelled with the category name. The scripts also calculates how many times that gene is present in a given assembly, giving the opportunity to explore copy number variations (info saved in `Assembly_count` column). The attributes for each gene are then extacted from the original GFF3 file and merged with the counts based on the gene name. The merged dataset is saved as `pangenome_screening.csv`. Genes from the original GFF3 file which were not merged with the counts dataframe (and thus not annotated in the pangenome) are saved in `not_annotated_features.csv`.  

The script then calculates some statistics (`statistics.csv`) and generate a matrix (`matrix.csv`) where the columns are the assemblies and the rows are the genes. The matrix is populated with the counts for each gene in each assembly. The scripts outputs to stdout the maximum value found in the matrix that can be used as sanity check. 


# core_dispensable_pseuogenes.py

The idea behind this script is exactly the same as `core_dispensable_genes.py` but there are substantial differences between these two scripts. 

In the first part of the script the PAF file is processed in the same way to filter the results based on coverage, estimated identity and jaccard index. 

Afterwards the GFF3 file is loaded to extract some information from the attributes, particularly the name of the pater gene and the number of stop codons. There are two more inputs to give to this script compared to the one designed for the genes: the name of the reference chromosome and the name of the reference assembly (e.g. GCA_000001735#1#CP002688.1 - chromosome name: CP002688.1 ; reference name: GCA_000001735). The chromosome name is used to keep only the pseudogenes annotated in the processed chromosome, i.e if we are processing the PAF obtained from the community corresponding to chr1, we will keep only the pseudogenes annotated on chr1. This is a further filter to make sure that the pseudogenes are on the same chromosome as the alignments, because we are merging only based on coordinates as pseudogenes don't have a real name but only a pater name.

Subsequently ambiguous pseudogenes are removed. The pseudogenes which are considered ambiguous are those that are not annotated on the reference (this is why we give the reference name as input), do not contain stop codons (info extracted from attributes) and have coverage greater or equal to 100. These are not considered valid as we don't really know if in the other assemblies these pseudogenes are really pseudogenes (they may contain the full functional gene).

The script continues performing the counting and the classification of genes as explained above. The only difference here is that we are not grouping by Feature and then for Assembly, but we group for Pater and then for Assembly. This approach allowed to map pseudogenes coming from the same pater but in different position as "different" pseudogenes in the `untangle` step, but to count at this stage how many times a "pseudogenised" gene was found in each assembly. We have in this way a copy number variation based on the pater genes for the pseudogenes. 

The scripts then produces files in the same format as the outputs of `core_dispensable_genes.py` to keep compatibility with downstream scripts and analysis. The statistics and the matrix are then generated as explained above. 


# filter_private_bed.py 

This Python script is designed to process a BED (Browser Extensible Data) file, which typically contains genomic intervals. The script filters the BED file based on a length threshold, categorises the genomic intervals into size categories, and outputs the filtered data to a new file. Here's a step-by-step explanation:

- The script imports necessary libraries: `pandas`, `numpy`, `argparse`, `logging`
- **Function `check_file_exists(filename)`**: Checks if the specified file exists. If not, logs an error message and returns `False`; otherwise, returns `True`.
-Function `parse_args()`: The script requires three arguments:
     - `-i`/`--input_file`: The input BED file path.
     - `-thr`/`--threshold`: A length threshold for filtering the BED file. Sequences shorter than this threshold will be separated and saved in `filtered_out.bed`.
     - `-o`/`--output_file`: The output file path for the sequences longer than threshold.
- **Function `process_bed(input_file, threshold, output_file)`**:
   - Reads the input BED file into a pandas DataFrame.
   - Calculates the length of each genomic interval (path) and categorizes these lengths into predefined size categories.
   - Counts the number of intervals in each size category and logs this information.
   - Filters the DataFrame to include only the intervals with lengths equal to or greater than the specified threshold. It then saves this filtered data to the output file.
   - Additionally, intervals that do not meet the threshold are saved to a separate file named `filtered_out.bed`.
- **Function `main()`**: calls the other functions


# split_and_extract_for_exonerate.sh

This Bash script is designed to process a "diamond results file" by splitting it into smaller chunks, extracting unique sequence and protein IDs from each chunk, and then using those IDs to extract corresponding sequences from FASTA files. Here's a step-by-step explanation:

- **Parse Command-Line Options**: The script accepts four options:
   - `-i` for specifying the path to the input diamond results file.
   - `-l` for specifying the number of lines each chunk should contain.
   - `-s` for specifying the path to the sequence FASTA file.
   - `-p` for specifying the path to the protein FASTA file.
   The options are parsed using a `while` loop and `getopts`.

- **Check Required Options**: The script checks if all required options (`-i`, `-l`, `-s`, `-p`) are provided. If any are missing, it displays the help message and exits.

- **Split the Input File**: The `split` command is used to divide the input diamond results file into smaller files (chunks), each containing a specified number of lines (`$lines_per_chunk`). These chunks are prefixed with `diamond_chunk_`.

- **Process Each Chunk**:
   - For each chunk, unique sequence IDs (first column of the chunk) and protein IDs (second column of the chunk) are extracted, sorted, and stored in separate files (`${chunk}_seq_ids.txt` and `${chunk}_protein_ids.txt`).
   - The `seqtk subseq` command is then used to extract sequences from the sequence FASTA file and protein FASTA file based on the IDs in `${chunk}_seq_ids.txt` and `${chunk}_protein_ids.txt`, respectively. The extracted sequences are saved in FASTA format files (`${chunk}_seqs.fasta` and `${chunk}_proteins.fasta`).
   - Optionally, the intermediate ID files are removed to save space.

This script is useful for handling large diamond results files by breaking them down into manageable pieces, extracting relevant data, and organizing that data into separate files for further analysis.


# PW_exonerate.pl

This Perl script is designed to facilitate the analysis of protein and genome sequences by using the `exonerate` tool. It processes input FASTA files and diamond results, then runs `exonerate` for each entry in the diamond results. Here's a breakdown of its functionality:

1. **Import Required Modules**: The script starts by importing the `strict` module for good coding practices and `Getopt::Long` for processing command-line options.

2. **Display Help Message**: The `show_help` function prints usage information for the script, explaining the required command-line arguments.

3. **Parse Command-Line Options**: It retrieves command-line options for protein FASTA file, genome FASTA file, diamond results file, and a suffix for naming output files. These options are essential for the script's operation.

4. **Check Required Options**: The script checks if all necessary command-line arguments are provided. If any are missing, it displays the help message and exits.

5. **Load FASTA Files**: The `load_fasta_as_hash` function reads sequences from FASTA files (protein and genome) and stores them in hashes. This allows for efficient access to sequences by their IDs.

6. **Process Diamond Results**: The `load_file_as_AoA` function loads the diamond results into an array of arrays (AoA), where each sub-array represents a line from the diamond results file, split by a separator (default is tab).

7. **Main Processing Loop**:
   - For each entry in the diamond results, it generates unique filenames for protein, genome, and output files using the provided suffix.
   - It writes the relevant protein and genome sequences (retrieved from the hashes loaded earlier) to their respective files.
   - It constructs and executes an `exonerate` command using the generated files. The `exonerate` tool is run with specific parameters to align the protein sequences to the genome sequences. The results are appended to a GFF file named according to the provided suffix.

8. **Utility Functions**:
   - `load_fasta_as_hash`: Reads a FASTA file and stores sequences in a hash, allowing for retrieval by sequence ID. Sequences can be stored as strings or arrays of characters.
   - `load_file_as_AoA`: Loads a file into an array of arrays, splitting each line by a specified separator.
   - `load_file_as_array`: Reads a file into an array, with each element representing a line from the file.
   - `make_array_from_scalar`: Splits a scalar (string) into an array using a specified separator.

This script is useful for batch processing of protein and genome sequences with `exonerate`, automating the generation of input files and collection of results.

# parse_exonerate.py 

This script was designed to parse exonerate results with the following goals:
- filter out results that do not meet minimal criteria (20% coverage and 50% identity)
- distinguish genes from pseudogenes (it considers pseudogenes those alignments that contain stop codons or do not reach 100% coverage)
- classify pseudogenes which align close to the edge of the alignments as "unclassified" (we do not know if these are fully contained in the assembly or not)
- rename the UniRef100 genes which have a correspondent in Araport annotation with the Araport ID

To achieve so, the scripts extract the `ryo` lines from exonerate results and creates a dataframe, that allows to calculate the coverage. It then parses separately the C4 alignments previously extracted (from where we get the information about the presence of stop codons) and creates a dataframe. These two dataframes are merged and: 
- the information about edge alignments is extracted (it's considered an edge alignment when 1. if it is in forward orientation, the target alignment begins at 0 or the target alignment end is equal to the target length; 2. if it is in reverse orientation, the target alignment begin is equal to the target length or the target alignment end is equal to 0)
- the filter of coverage at least 20% and identityy at least 50% is applied
- the results are classified as genes if the coverage is equal to 100% and stop codons are 0
- the results are classified as pseudogenes if the coverage is lower than 100% or there are stop codons.

Subsequently, ambiguous pseudogenes are removed. Pseudogenes are considered ambiguous if they do not contain stop codons and they are edge alignments (in this situation we are not able to know if they are really pseudogenes or full genes). These are classified as 'unclassified'.

The information previously extracted from the UniRef100 headers are added back to the dataframe. 

The UniProt to Araport conversion table is loaded. For the conversion, only Arapost loci names are taken into account (splicing variants are ignored). The data are then grouped by UniProt_code, in order to have single UniProt IDs, and eventually join Araport loci if more than one locus correspond to the same UniProt ID. The processed table is saved in `conversion_table.tsv`.

The conversion table is merged with the dataframe. From here, those features that have a correspondent in TAIR10 annotations are subset and saved in `corresponding_to_tair.tsv`. Those features that do not have a TAIR10 correspondent annotation are subset and saved to `new_results_not_tair`.


# review_exonerate_results.py


This Python script is designed to filter and review exonerate results based on UniProt protein IDs. It uses the pandas library for data manipulation and argparse for command-line argument parsing. Here's a detailed breakdown of its functionality:

### 1. **Imports and Logging Setup**
- Imports necessary libraries: pandas for data handling, argparse for parsing command-line arguments, and logging for logging messages.
- Configures logging to display timestamps, log level, and messages.

### 2. `get_uniprot_results`
- **Purpose**: To read a UniProt ID mapping file and extract unique UniRef100 IDs that correspond to reviewed proteins.
- **Parameters**: `id_mapping_file`- The path to the UniProtKB/Swiss-Prot mapping file.
- **Process**:
  - Reads the mapping file into a pandas DataFrame.
  - Extracts unique values from the 'From' column, which contains UniRef100 IDs.
- **Returns**: A list of unique UniRef100 IDs.

### 3. `review_results`
- **Purpose**: To filter a results file, excluding entries that do not match the reviewed UniProt IDs.
- **Parameters**:
  - `not_ref_results` - Path to the results file generated by `parse_exonerate.py`.
  - `uniprot_reviewed_ids` - List of reviewed UniProt IDs to keep.
  - `output_file` - Path for the output file to save filtered results.
- **Process**:
  - Reads the results file into a pandas DataFrame.
  - Filters the DataFrame, keeping only rows where the 'Query_ID' matches any of the reviewed UniProt IDs.
  - Saves the filtered DataFrame to the specified output file.

### 4. **Main Function**
- Sets up command-line arguments for the script:
  - `-i`/`--id_mapping_file`: The path to the UniProtKB/Swiss-Prot mapping file.
  - `-a`/`--not_ref_results`: The path to the results file to be reviewed.
  - `-o`/`--output_file`: The path for the output file where reviewed results will be saved.
- Parses command-line arguments.
- Calls  the other functions.

The script is used to refine the output of the previous exonerate analysis by excluding proteins that have not been reviewed in UniProt. This is achieved by filtering the results against a list of reviewed UniProt IDs, ensuring that the final output contains only entries of interest. 


# final_genes_screening.py

This script joins the results coming from `odgi untangle` analysis with those obtained from `exonerate`. It first formats all the results to obtain a PAF-like file for all of them, and then concatenates the dataframes to count the occurrences of the genes across the assemblies. It also keeps track of the estimated identity and coverage for each annotated gene, plus the way it was annotated (untangle and/or exonerate). The genes are then labelled according to their presence across the assemblies as core (present in all assemblies `-na`), softcore (not found in all the assemblies but in a number of assemblies above or equal to the softcore limit `-sl`), dispensable (found in a nuumber of assemblies below the softcore limit but at least in 2), private (when they are found in only one assembly). The attributes for each gene are then added back from the input files. 

The scripts then calculates statistics about how many genes were found in each compartment and their proportions. 

The counts per assembly (i.e. how many times that gene is found in a specific assembly), in this case is not precise, because in some cases the same gene is annotated by both untangle and exonerate. The analysis with exonerate was done only for the sequences not covered by the reference and some of these regions were enlarged to include all the possible variations. This implied a change in coordinates and the impossibility of determine if a gene annotated with exonerate was the same already annotated by untangle. For all these reasons, the matrix is produced by this script as a presence-absence matrix, and does not take into account CNVs.

# final_pseudogenes_screening.py

The logic behind this script is exactly the same used for `final_genes_screening.py`, with minor differences to adapt the script to the format of the previous analysis on pseudogenes. 

# nodes_matrix_processing.py

This Python script is designed for processing and analysing a node coverage matrix. 

1. **Imports and Setup**: The script imports the necessary libraries (`pandas` for data manipulation and `argparse` for command-line argument parsing).

2. `matrix_prep`:
   - **Purpose**: Prepares the node matrix by cleaning and aggregating the data.
   - **Steps**:
     - Reads a tab-separated values (TSV) file into a DataFrame.
     - Drops unnecessary columns specified by name.
     - Sorts the DataFrame based on the 'group.name' column.
     - Aggregates the data by summing up node columns for each group, identified by 'group.name'.
     - Renames node columns by adding a specified prefix.
     - Saves the aggregated data to a new TSV file.

3. `count_and_label`
   - **Purpose**: Counts occurrences of nodes across assemblies and labels them based on their prevalence.
   - **Steps**:
     - Counts how many rows have a count greater than 0 for each node.
     - Calculates the total sum for each node across all assemblies.
     - Creates a new DataFrame with 'PAV' (Presence/Absence Variation) and 'CNV' (Copy Number Variation) columns.
     - Labels each node as 'core', 'softcore', 'private', or 'dispensable' based on its count relative to the total number of assemblies and a specified softcore limit.
     - Saves this information to a TSV file.

4. `get_statistics`:
   - **Purpose**: Calculates and saves statistics about the node classifications.
   - **Steps**:
     - Calculates the total number of features and the number and percentage of features classified as core, softcore, dispensable, and private.
     - Stores these statistics in a dictionary, converts it to a DataFrame, and saves it to a TSV file.


# transpose.pl

This Perl script is designed to transpose a matrix stored in a file.

1. **Use Strict Mode**: The script starts with `use strict;` to enforce some good programming practices, such as requiring all variables to be declared.

2. **Load the Matrix**: It loads the matrix from a file specified by the first command-line argument (`$ARGV[0]`) into an array of arrays (`@AoA`) using the `load_file_as_AoA` subroutine. This subroutine reads the file into an array of arrays, where each sub-array represents a row of the matrix.

3. **Transpose the Matrix**:
   - The script iterates over the columns of the matrix (outer loop with index `$j`) and then over the rows (inner loop with index `$i`).
   - For each element, it prints the value located at the `[row][column]` position in the original matrix, effectively accessing the `[column][row]` position in the transposed matrix.
   - It prints a tab character (`\t`) between elements of the same row and a newline character (`\n`) at the end of each row to format the output as a transposed matrix.

4. **Subroutines**:
   - `load_file_as_array`: Reads the entire file into an array, where each element of the array corresponds to a line in the file.
   - `make_array_from_scalar`: Splits a scalar (string) into an array using a specified separator (default is a tab character). This is used to split each line of the file into individual elements of a row.
   - `load_file_as_AoA`: Combines `load_file_as_array` and `make_array_from_scalar` to read the file into an array of arrays. Each line is split into its constituent elements based on the separator, forming the rows of the matrix.

The script assumes that the input file is formatted as a matrix with rows separated by newline characters and columns separated by tabs (or another specified separator).

# comp_node_path_CNV.pl

This Perl script reads input from two files specified by the user:
- the list of features (genes, pseudogenes or node) with their classification into genomic classes (core, softcore, dispensable, private)
- CNV matrix.

It outputs the results in a format ready to be visualised in R. 

1. **Load a Table as a Hash**: The script starts by loading a table from the first input file (`$ARGV[0]`) into a hash (`%hash`) using the `load_table_as_hash` subroutine. This table likely contains genomic data mapped to specific categories (e.g., core, softcore, dispensable, private, missing).

2. **Process the Second Input File**: It opens the second input file (`$ARGV[1]`) and reads it line by line. The first line is treated as a header, from which it extracts column names into an array (`@array`), excluding the first column name.

3. **Data Processing Loop**:
   - For each subsequent line, it extracts the path name (i.e. assembly name) using `take_name` and the nodes using `take_nodes`.
   - It then counts the types of nodes (`core`, `softcore`, `dispensable`, `private`) based on the previously loaded hash and the nodes extracted from the current line, using the `count_node_types` subroutine.
   - The counts are joined into a string and printed alongside the path name.
   - Additionally, it calls `print_R_ready` to append the counts in a format suitable for R analysis to a file named "data_set_R.txt".

4. **Subroutines**:
   - `take_name` and `take_nodes` are utility functions to extract the path name and nodes from a line of text, respectively.
   - `count_node_types` computes the frequency of each node type based on the input nodes and the hash loaded from the first file.
   - `load_table_as_hash` loads a table from a file into a hash, where each key-value pair corresponds to a row in the table, with the key being the first column's value.
   - `load_file_as_array` and `load_file_as_AoA` (Array of Arrays) are helper functions to read a file into an array or an array of arrays, respectively, for easier data manipulation.
   - `print_R_ready` formats the counts of node types for each path into a string and appends it to a file, presumably for analysis in R.

# comp_node_path_PAV.pl

This Perl script reads input from two files specified by the user:
- the list of features (genes, pseudogenes or node) with their classification into genomic classes (core, softcore, dispensable, private)
- PAV matrix.

It outputs the results in a format ready to be visualised in R. 

1. **Initialization and Input Handling**:
   - It loads a table from the first input file specified by the user, i.e. the list file, (`$ARGV[0]`) into a hash (`%hash`) using the `load_table_as_hash` subroutine. 

2. **Reading and Processing the Second Input File**:
   - Opens the second input file (`$ARGV[1]`), i.e. the matrix, and reads its first line to extract column names (excluding the first column) into an array (`@array`).
   - Iterates over each subsequent line of the file, processing it to extract the path name (`$path_name`) and nodes (`@nodes`) using `take_name` and `take_nodes` subroutines, respectively.

3. **Node Type Counting and Output**:
   - For each line (representing a path, i.e. an assembly), it counts the types of nodes (e.g., core, softcore, dispensable, private) based on the previously loaded hash and the nodes extracted from the current line. This is done using the `count_node_types` subroutine.
   - The counts are joined into a string (`$stampo`) and printed alongside the path name.
   - Additionally, it calls `print_R_ready` to append the counts in a format suitable for R analysis to a file named "data_set_R.txt".

4. **Subroutines**:
   - `take_name` and `take_nodes` extract the path name and nodes from a line of text, respectively.
   - `count_node_types` computes the frequency of each node type based on the input nodes and the hash loaded from the first file.
   - `load_table_as_hash` loads a table from a file into a hash, where each key-value pair corresponds to a row in the table, with the key being the first column's value.
   - `load_file_as_array` and `load_file_as_AoA` (Array of Arrays) are helper functions to read a file into an array or an array of arrays, respectively, for easier data manipulation.
   - `print_R_ready` formats the counts of node types for each path into a string and appends it to a file, presumably for analysis in R.

# combine_transposed_matrices.py

This Python script combines multiple transposed matrices into a single matrix using the pandas library. 

1. **Import pandas**: The script starts by importing the pandas library, which is a powerful tool for data manipulation and analysis.

2. **Load DataFrames**:
   - It loads transposed matrices from five text files (`chr1_transp_matrix.txt`, `chr2_transp_matrix.txt`, etc.) into pandas DataFrames ([`df1`], [`df2`], etc.). Each matrix file is expected to be tab-separated ([`sep='\t'`], and the first column of each file is set as the index of the DataFrame ([`index_col=0`]).

3. **Print Columns Information**:
   - Before concatenating, it prints the number of columns in each DataFrame to give an overview of the data structure. 

4. **Concatenate DataFrames**:
   - The script then concatenates these DataFrames horizontally ([`axis=1`]), meaning it aligns them side by side based on their index. This results is placed in a single DataFrame ([`result`]) that combines all the columns from the individual matrices.

5. **Save the Result**:
   - The concatenated DataFrame is saved to a new file (`concatenated_matrix.csv`), using a tab as the separator ([`sep='\t'`]).

6. **Print Final DataFrame Information**:
   - Finally, it prints the number of columns in the concatenated DataFrame to confirm the successful combination of the matrices.

This script is very specific to *A. thaliana* data and to our analysis.

# combine_notransp_matrices.py

This Python script combines multiple matrices stored in text files into a single matrix and saves the result to a new file. 

1. **Import pandas library**: The script starts by importing the pandas library, which is a powerful tool for data manipulation and analysis.

2. **Load matrices into DataFrames**:
   - It loads five matrices from text files (`matrix_chr1.txt` to `matrix_chr5.txt`). Each matrix is stored in a separate file.
   - Each file is read into a pandas DataFrame using [`pd.read_csv`], with tab (`'\t'`) as the separator and the first column set as the index of the DataFrame.

3. **Print the number of rows for each DataFrame**:
   - It creates a list of the loaded DataFrames.
   - Then, it iterates over this list, printing the number of rows ([`shape[0]`]) for each DataFrame.

4. **Concatenate the DataFrames horizontally**:
   - The script concatenates the five DataFrames along the rows (axis=0) using [`pd.concat`]. This means it stacks the DataFrames on top of each other, combining them into a single DataFrame.
   - Note: The comment suggests horizontal concatenation, but the code actually performs vertical concatenation (stacking them vertically). For horizontal concatenation, [`axis`] should be set to 1.

5. **Save the concatenated DataFrame**:
   - The combined DataFrame is saved to a new file named `concatenated_matrix.csv`, using tab as the separator.

6. **Print the number of rows in the final DataFrame**:
   - Finally, it prints the number of rows in the concatenated DataFrame, indicating the size of the final combined matrix.

This script is very specific to *A. thaliana* data and to our analysis.

# BP_TopGO.R

This is the script used to conduct the gene ontology enrichment analysis for biological processes in R studio.

# shuffle_count.pl

This script was used to obtain the dataset for the simulation of the pangenome growth curve based on gene data (PAV).

It performs the following tasks:
1. Loads a matrix from a file.
2. Iterates over different sample sizes and performs multiple iterations for each sample size.
3. Shuffles the matrix, counts genes, and classifies genes.
4. Calculates averages and standard deviations of the counts and classifications.
5. Writes the results to an output file.


# curve_plot2.py

This script was used to obtain the figure of the simulated pangenome growth curve.

It reads data from a specified file, processes it to extract sample sizes and datasets (average values and standard deviations), and plots the data with error bars using Matplotlib. The script is designed to be run from the command line with the filename as an argument.