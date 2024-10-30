#!/usr/bin/env python3

import pandas as pd

# Step 1: Load each matrix into a DataFrame, setting the first column as the index
df1 = pd.read_csv('matrix_chr1.txt', sep='\t', index_col=0)
df2 = pd.read_csv('matrix_chr2.txt', sep='\t', index_col=0)
df3 = pd.read_csv('matrix_chr3.txt', sep='\t', index_col=0)
df4 = pd.read_csv('matrix_chr4.txt', sep='\t', index_col=0)
df5 = pd.read_csv('matrix_chr5.txt', sep='\t', index_col=0)
# Add more DataFrames as needed

# Print number of rows for each DataFrame
dataframes = [df1, df2, df3, df4, df5]
for i, df in enumerate(dataframes, start=1):
    print(f"DataFrame {i} has {df.shape[0]} rows.")

# Step 2: Concatenate the DataFrames vertically
result = pd.concat([df1, df2, df3, df4, df5], axis=0)  # Add more DataFrames as needed

# Step 3: Save the concatenated DataFrame to a new file
result.to_csv('concatenated_matrix.csv', sep='\t')

# Print number of columns and head for the final DataFrame
print(f"Final concatenated DataFrame has {result.shape[0]} rows.")
