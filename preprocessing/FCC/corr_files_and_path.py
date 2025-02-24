#!/usr/bin/env python
# coding: utf-8

# In[1]:


import csv
import os
import re


# In[2]:


# output file name and path
corr_files = '/mnt/DATA3/timo/main/body/outputs/corr_files.tsv'


# In[3]:


directory = 'main/body/outputs/corr_tissues/'
data=[]

files = os.listdir(directory)

ee_numbers = []
for file in files:
    match = re.search(r'EE\d+', file)
    if match:
        ee_numbers.append(match.group())


# In[4]:


for i in range(len(ee_numbers)):
    path = "main/body/outputs/corr_tissues/" + files[i]
    ee_number = ee_numbers[i]
    data.append((ee_number,path))


# In[5]:


with open(corr_files, mode = 'w', newline='') as file:
    writer = csv.writer(file, delimiter='\t')  # Set delimiter to tab
    writer.writerow(['ID', 'path'])  # Write header row
    writer.writerows(data)  # Write data rows
    print("done")

