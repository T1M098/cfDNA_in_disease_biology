import pandas as pd
import gzip
import numpy as np
import glob
import os
import argparse

# Set up argument parsing
parser = argparse.ArgumentParser(description="Process expression data for a specific sample")
parser.add_argument('sample', type=str, help='Sample name, e.g., EE87888')
args = parser.parse_args()

# Use the sample from the argument
sample = args.sample

# Define paths
dir_path = f"/main/body/{sample}/smoothed/" 
file_pattern = "*.tsv.gz"
outdir_path = "/main/body/corr_tissues/" 

# Load expression data
expression = pd.read_csv("/avg_expression_Tabula_Sapiens.tsv", sep='\t')
# Reset the index and rename it to 'GeneID'
expression.reset_index(inplace=True)
expression.rename(columns={'index': 'GeneID'}, inplace=True)

# Load gene annotation file, containing ENSG, start and end position
gene_mapping = pd.read_csv("/final_annotation_file.tsv",
                           names=['GeneID','chromosome', 'start', 'end', 'strand'],
                           header=None, delimiter='\t')
gene_mapping = gene_mapping.iloc[1:].reset_index(drop=True)

# Add gene start and end positions to gene expression matrix
merged_data = pd.merge(gene_mapping, expression, on='GeneID')
filtered_merged_data = merged_data[~merged_data['chromosome'].str.contains('random|Un')]
expression_data_final = filtered_merged_data[~filtered_merged_data['chromosome'].isin(['M', 'X', 'Y'])]
expression_data_final = expression_data_final.drop_duplicates(subset='GeneID', keep='first')
expression_data_final['start'] = pd.to_numeric(expression_data_final['start'], errors='coerce')
expression_data_final['end'] = pd.to_numeric(expression_data_final['end'], errors='coerce')

# Add new column for mean_center_score
final_results = expression_data_final.copy()
final_results.insert(5, 'mean_center_score', np.nan)
final_results.reset_index(inplace=True)

left_out = []
# Go through files, calculate mean FCC for each gene and fill in values mean_center_score
for ENSG_file in os.listdir(dir_path):
    ENSG = ENSG_file.split('_')[1].split('.')[0]
    header = ["chromosome", "position", "covCount", "startCount", "CenterCount", "smoothed_CenterCount"]
    fragment_data = pd.read_csv(dir_path + ENSG_file, compression="gzip", sep="\t", header=None, names=header)
    fragment_data['chromosome'] = "chr" + fragment_data['chromosome'].astype(str)
    if ENSG in final_results["GeneID"].values:
        for index, row in final_results[final_results['GeneID'] == ENSG].iterrows():
            scores_within_range = fragment_data[
                (fragment_data['position'] >= row['start']) &
                (fragment_data['position'] <= row['end'])
            ]
            mean_score = scores_within_range['smoothed_CenterCount'].mean() if not scores_within_range.empty else np.nan
            final_results.at[index, 'mean_center_score'] = mean_score                
    else:
        left_out.append(ENSG)

final_results.dropna(inplace=True)

# Calculate correlation of 'mean_center_score' with each column from the 7th column onward
correlations_mean = final_results.iloc[:, 7:].apply(lambda x: final_results['mean_center_score'].corr(x))

# Convert to DataFrame and reset index
correlations_mean = correlations_mean.reset_index()

# Rename columns for clarity
correlations_mean.columns = ['cell_type', 'correlation_mean']
correlations_mean['rank_mean'] = correlations_mean['correlation_mean'].rank(method='first', ascending=True).astype(int)

# Sort by correlation values
correlations_mean = correlations_mean.sort_values(by='correlation_mean', ascending=True)

correlations_mean.to_csv(outdir_path +  f'{sample}_correlations_center.csv', index=False)