#!/usr/bin/env python
# coding: utf-8
# Linear regression on ranks for top up regulated cell types with sampled tumor fraction

# Load packages
import re
import pandas as pd
import numpy as np
import math
import os

input_dir = "/preeclampsia/main_centers/body/corr_tissues_new/"
output_dir = "/preeclampsia/main_centers/body/outputs/"

# in meta data preeclampsia 2 with ATL concentrations samp 57 is healthy (in previous meta data it was preeclampsia)
meta_dat = pd.read_csv("/mnt/DATA3/timo/data/preeclampsia_metadata_2.tsv", sep = "\t")
meta_dat = meta_dat[["sample name", "ALT_(U/l)_high_above_33", "group"]]
meta_dat = meta_dat.rename(columns={"group": "status", "sample name" : "sample", "ALT_(U/l)_high_above_33" : "ALT"})


meta_dat[meta_dat['sample'] == 'samp57']

# Create the empty DataFrame
columns = ["sample", "status", "ALT", "cell_type", "correlation", "rank"]
empty_df = pd.DataFrame(columns=columns)

for file in os.listdir(input_dir):
    sample = file.split('_')[0]
    filename = input_dir + file
    data = pd.read_csv(filename)
    data["sample"] = sample
    data = data.rename(columns={"correlation_mean": "correlation", "rank_mean" : "rank"})
    data_merged = pd.merge(data, meta_dat, on = "sample")
    data_merged = data_merged[columns]
    empty_df = pd.concat([empty_df, data_merged], ignore_index = True)
empty_df = empty_df.rename(columns={"cell_type": "cell_type_tissue"})
empty_df

# Remove rows where the cell_type_tissue contains 'RNA.NA.' 
empty_df = empty_df[~empty_df['cell_type_tissue'].str.startswith('RNA.NA')]

empty_df[empty_df['sample'] == 'samp57']

# Function to replace dots before specific words and the last dot if no underscores are present
def replace_dot_before_words(s):
    # Replace dot before specific words (e.g., Bone, Small, etc.)
    s = re.sub(r'\.(?=Bone|Small|Large|Salivary|Lymph)', '_', s)
    
    # If no underscore is present, replace the last dot in the string with an underscore
    if '_' not in s:
        # Replace the last dot before the last word with an underscore
        s = re.sub(r'\.(?=[^\._]*$)', '_', s)
    
    return s

empty_df.loc[:, 'cell_type_tissue'] = empty_df['cell_type_tissue'].apply(replace_dot_before_words)


# Function to clean the cell type strings
def clean_cell_type(cell_type):
    # Replace all "." and ".." with spaces
    cell_type = re.sub(r'\.\.+', ' ', cell_type)
    cell_type = cell_type.replace('.', ' ')
    
    # Remove "RNA" from the string
    cell_type = cell_type.replace("RNA", "")
     
    return cell_type

# Apply the function to the 'cell_type' column
empty_df.loc[:, 'cell_type_tissue'] = empty_df['cell_type_tissue'].apply(clean_cell_type)

len(empty_df["cell_type_tissue"].unique())

empty_df.to_csv(output_dir + "Ranks_FCC_RNA_cell_type_tissue_preeclampsia.csv", index=False)