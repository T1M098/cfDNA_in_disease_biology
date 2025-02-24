#!/usr/bin/env python
# -*- coding: ASCII -*-

"""

:Author: Martin Kircher
:Contact: mkircher@uw.edu
:Date: *05.07.2014

"""

"""
:Modified Version

:Author: Timo Rieger
:Contact: trieger@student.ethz.ch
:Date: *24.09.2024
"""


import sys, os
import pysam
import gzip
from optparse import OptionParser
from collections import defaultdict
from pathlib import Path

parser = OptionParser()
parser.add_option("-i","--individual", dest="individual", help="Individual to use (comma-separate, def 'A11')",default="A11")
parser.add_option("-p","--project", dest="project", help="Project name to use (def 'TSS')",default="TSS")
parser.add_option("-t","--temp", dest="temp", help="Temp directory with blocks (def '/net/shendure/vol1/home/mkircher/nucleosomes/expression')",default="/net/shendure/vol1/home/mkircher/nucleosomes/expression")
parser.add_option("-r","--root", dest="root", help="Root directory for outputs (def '/net/shendure/vol1/home/mkircher/nucleosomes/expression')",default="/net/shendure/vol1/home/mkircher/nucleosomes/expression")
parser.add_option("-a","--annotation", dest="annotation", help="Use regions transcript file (def transcriptAnno.tsv)",default="transcriptAnno.tsv")
parser.add_option("-s","--suffix", dest="suffix", help="Suffix for fft_summaries folder (def '')",default="")
(options, args) = parser.parse_args()

rootDir = options.root
tempDir = options.temp

project = options.project
ind = options.individual

outfileCov = gzip.open("%s/%s/fft_summaries%s/fft_%s_cov.tsv.gz"%(rootDir,project,options.suffix,ind),'w')
outfileStarts = gzip.open("%s/%s/fft_summaries%s/fft_%s_starts.tsv.gz"%(rootDir,project,options.suffix,ind),'w')
outfileWPS = gzip.open("%s/%s/fft_summaries%s/fft_%s_WPS.tsv.gz"%(rootDir,project,options.suffix,ind),'w')

geneIDs = set()
if os.path.exists(options.annotation):
  infile = open(options.annotation)
  for line in infile:
    cid,chrom,start,end,strand = line.split() # positions are 1-based and inclusive
    geneIDs.add(cid)
  infile.close()

wrote_header = False
# goes through all cid of the annotation file / through all block files
cid_count = 0
cid_counter = 0
sample_dir = os.path.join(tempDir, project, ind, "fft")
sample_dir=Path(sample_dir)
# Count the number of files in the folder
file_count = len(list(sample_dir.glob('*.tsv.gz')))

print("In progress...")
for cid in sorted(geneIDs):
  # checking if the block file exists 
  if os.path.exists("%s/%s/%s/fft%s/block_%s.tsv.gz"%(tempDir,project,ind,options.suffix,cid)):
    # writing the header if not already done
    if not wrote_header:
      # changed it so "#Region" is in bytes from: header =  ['#Region']
      header = [b'#Region']
      # read the block file containing these columns: frequency, coverage, starts and WPS
      # added mode = 'rt' so the content are strings instead of bytes
      infile = gzip.open("%s/%s/%s/fft%s/block_%s.tsv.gz"%(tempDir,project,ind,options.suffix,cid))#, mode = 'rt')
      infile.readline()
      # for each line in the block file take the frequency and append to header list
      for line in infile:
        fields = line.split()
        header.append(fields[0])
      infile.close()
      #print(fields)
      #print(header)
      #print(type(header))
      #print(header[1])
      #print(type(header[0]))
      # changed these lines due to error when writing in bytes, need to define it also for "\t" and "\n" -> b"\t"
      outfileCov.write(b"\t".join(header)+b"\n")
      outfileStarts.write(b"\t".join(header)+b"\n")
      outfileWPS.write(b"\t".join(header)+b"\n")
      wrote_header=True
    infile = gzip.open("%s/%s/%s/fft%s/block_%s.tsv.gz"%(tempDir,project,ind,options.suffix,cid))  
    # changed cid to bytes (from strings)
    #print(type(cid), cid)
    # Convert cid to bytes
    cid_bytes = cid.encode('utf-8')
    #print(type(cid_bytes), cid_bytes)
    covs = [cid_bytes]
    starts = [cid_bytes]
    wps = [cid_bytes]
    infile.readline()
    for line in infile:
      fields = line.split()
      covs.append(fields[1])
      starts.append(fields[2])
      wps.append(fields[3])
    infile.close()
    #print(covs)
    # changed these lines due to error when writing from bytes
    outfileCov.write(b"\t".join(covs)+b"\n")
    outfileStarts.write(b"\t".join(starts)+b"\n")
    outfileWPS.write(b"\t".join(wps)+b"\n")
    
    # display progress
    cid_count += 1
    cid_counter += 1
    if cid_counter == 1000:
        print(f"{cid_count}/{file_count}")
        cid_counter = 0

print("Done.")    
outfileCov.close()
outfileStarts.close()
outfileWPS.close()
