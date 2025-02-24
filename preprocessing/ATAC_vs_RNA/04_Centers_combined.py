#!/usr/bin/env python
# coding: utf-8

import re
import pandas as pd
import numpy as np
import math
import os

rdataset_type = "ATAC" # RNA
input_dir = "ATAC/correlation_files/" + rdataset_type + "/"
output_dir = "ATAC/outputs/"
samples_to_exclude =["EE87920", "EE87921",
    "EE87947", "EE87959", "EE87963", "EE87998", 
    "EE88009", "EE88015", "EE88020", "EE88025", 
    "EE88027", "EE88082", "EE88093", "EE88127", 
    "EE88136", "EE88141", "EE88152", "EE88164", "EE87892", "EE88269"]

meta_dat = pd.read_csv("Cristiano_samplemap.tsv", sep = "\t")
meta_dat = meta_dat[["sample", "tfx", "Patient Type", "Gender"]]
meta_dat = meta_dat.rename(columns={"Patient Type": "status", "Gender" : "gender"})

# Create the empty DataFrame
columns = ["sample", "status", "tfx", "gender", "cell_type", "correlation", "rank"]
empty_df = pd.DataFrame(columns=columns)

for file in os.listdir(input_dir):
    sample = file[0:7]
    if sample not in samples_to_exclude:
        filename = input_dir + file
        data = pd.read_csv(filename)
        data["sample"] = sample
        data = data.rename(columns={"correlation_mean": "correlation", "rank_mean" : "rank"})
        data_merged = pd.merge(data, meta_dat, on = "sample")
        data_merged = data_merged[columns]
        empty_df = pd.concat([empty_df, data_merged], ignore_index = True)

empty_df["sample"].nunique()

empty_df.to_csv(output_dir + "Ranks_FCC_" + rdataset_type + "_cell_type.csv", index=False)