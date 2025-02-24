#!/usr/bin/env python
# coding: utf-8
import csv
import os
import re

# output file name and path
corr_files = '/preeclampsia/main_centers/body/outputs/corr_files_preeclampsia.tsv'

directory = '/preeclampsia/main_centers/body/corr_tissues'
data=[]

files = os.listdir(directory)
samp_numbers = []
for file in files:
    match = re.search(r'samp\d+', file)  # Regex to find patterns like EE followed by numbers
    if match:
        samp_numbers.append(match.group())
# print(samp_numbers)

for i in range(len(samp_numbers)):
    path = "/preeclampsia/main_centers/body/outputs/corr_tissues_new/" + files[i]
    #print(path)
    samp_number = samp_numbers[i]
    #print(samp_number)
    data.append((samp_number,path))
# print(data)

with open(corr_files, mode = 'w', newline='') as file:
    writer = csv.writer(file, delimiter='\t')  # Set delimiter to tab
    writer.writerow(['ID', 'path'])  # Write header row
    writer.writerows(data)  # Write data rows
    print("done")