# Script that generates a gene expression matrix from Cufflinks ouputs  
# Reads the FPKM values from multiples samples and megres them into a single matrix 
# Last edited: 04/04/2026. Rohit

# Importing the necessary libraries 
import pandas as pd
import glob
import os

# Defining al the paths
path = "/scratch/alice/r/rt334/steered_sra/cufflinks_output_v2/*_cufflinks/genes.fpkm_tracking"
output = "/scratch/alice/r/rt334/steered_sra/expression_matrix_full.csv"

# Finding the files
files = sorted(glob.glob(path))
total_files = len(files)

print(f"Found {total_files} of the samples. Starting merge...")

# Processing the files
df_list = []
gene_map = None
count = 0

for f in files:
    sample = os.path.basename(os.path.dirname(f)).replace("_cufflinks", "")

    temp_df = pd.read_csv(
        f,
        sep='\t',
        usecols=['gene_id', 'gene_short_name', 'FPKM']
    )

    # Remove duplicate gene IDs
    temp_df = temp_df.drop_duplicates(subset='gene_id')

    # Save gene_id → gene_name mapping (once)
    if gene_map is None:
        gene_map = temp_df[['gene_id', 'gene_short_name']].drop_duplicates()
        gene_map.set_index('gene_id', inplace=True)

    # Keep expression only
    temp_df = temp_df[['gene_id', 'FPKM']]
    temp_df.columns = ['gene_id', sample]
    temp_df.set_index('gene_id', inplace=True)

    df_list.append(temp_df)

    count += 1
    if count % 50 == 0:
        print(f"Processed {count}/{total_files} samples...")

# Building the matrix
print("Merging all samples now")

expression_matrix = pd.concat(df_list, axis=1, sort=False)
expression_matrix = expression_matrix.fillna(0)

# Adds the gene names in the right position of the expression matrix
expression_matrix = expression_matrix.merge(
    gene_map,
    left_index=True,
    right_index=True,
    how='left'
)

# Reordering the columns
cols = ['gene_short_name'] + [c for c in expression_matrix.columns if c != 'gene_short_name']
expression_matrix = expression_matrix[cols]

# Resetting the index so gene_id becomes column
expression_matrix.reset_index(inplace=True)

# Rename columns for clarity
expression_matrix.rename(columns={
    'gene_id': 'gene_id',
    'gene_short_name': 'gene_name'
}, inplace=True)

# Saving the expression matrix as a csv file
expression_matrix.to_csv(output, index=False)

print(f"SUCCESS! Matrix saved to: {output}")
