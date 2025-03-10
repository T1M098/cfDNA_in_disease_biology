#!/usr/bin/env python

"""
:Original Version

:Author: Martin Kircher
:Contact: mkircher@uw.edu
:Date: *03.06.2014
"""

# modified
"""
:Modified Version

:Author: Timo Rieger
:Contact: trieger@student.ethz.ch
:Date: *24.09.2024
"""

import sys, os
from optparse import OptionParser
import gzip
import pysam
import random

from collections import defaultdict
from bx.intervals.intersection import Intersecter, Interval

def isSoftClipped(cigar):
    # Op BAM Description
    # M 0 alignment match (can be a sequence match or mismatch)
    # I 1 insertion to the reference
    # D 2 deletion from the reference
    # N 3 skipped region from the reference
    # S 4 soft clipping (clipped sequences present in SEQ)
    # H 5 hard clipping (clipped sequences NOT present in SEQ)
    # P 6 padding (silent deletion from padded reference)
    # = 7 sequence match
    # X 8 sequence mismatch
    for (op, count) in cigar:
        if op in [4, 5, 6]:
            return True
    return False

def aln_length(cigarlist):
    tlength = 0
    for operation, length in cigarlist:
        if operation == 0 or operation == 2 or operation == 3 or operation >= 6:
            tlength += length
    return tlength

parser = OptionParser()
parser.add_option("-i", "--input", dest="input", help="Use regions transcript file (def transcriptAnno.tsv)", default="transcriptAnno.tsv")
parser.add_option("-l", "--lengthSR", dest="lengthSR", help="Length of full reads (default 76)", default=76, type="int")
parser.add_option("-m", "--merged", dest="merged", help="Assume reads are merged (default Off)", default=False, action="store_true")
parser.add_option("-t", "--trimmed", dest="trimmed", help="Assume reads are trimmed (default Off)", default=False, action="store_true")
parser.add_option("-w", "--protection", dest="protection", help="Base pair protection assumed for elements (default 120)", default=120, type="int")
parser.add_option("-o", "--outfile", dest="outfile", help="Outfile prefix (def 'block_%s.tsv.gz')", default='block_%s.tsv.gz')  # reserve at least 6 digits
parser.add_option("-e", "--empty", dest="empty", help="Keep files of empty blocks (def Off)", default=False, action="store_true")
parser.add_option("--minInsert", dest="minInsSize", help="Minimum read length threshold to consider (def None)", default=-1, type="int")
parser.add_option("--maxInsert", dest="maxInsSize", help="Minimum read length threshold to consider (def None)", default=-1, type="int")
parser.add_option("--max_length", dest="max_length", help="Assumed maximum insert size (default 1000)", default=1000, type="int")
parser.add_option("--downsample", dest="downsample", help="Ratio to down sample reads (default OFF)", default=None, type="float")
parser.add_option("-v", "--verbose", dest="verbose", help="Turn debug output on", default=False, action="store_true")
(options, args) = parser.parse_args()

minInsSize, maxInsSize = None, None
if options.minInsSize > 0 and options.maxInsSize > 0 and options.minInsSize < options.maxInsSize:
    minInsSize = options.minInsSize
    maxInsSize = options.maxInsSize
    sys.stderr.write("Using min/max length cutoffs: %d/%d\n" % (minInsSize, maxInsSize))

options.outfile = options.outfile.strip("""\'""")

protection = options.protection // 2

validChroms = set(map(str,list(range(1,23))+["X","Y"]))

step_process = 0
total_files_removed = 0

if os.path.exists(options.input):
    infile = open(options.input)
    total_regions = sum(1 for _ in infile)  # Count total regions for progress indication
    infile.seek(0)  # Reset file pointer to the beginning of the file
    
    print(f"Total regions to process: {total_regions}")

 
#else:
    #print("Path to annotation file doesn't exist.")

    
    for idx, line in enumerate(infile):
        cid, chrom, start, end, strand = line.split()  # positions should be 0-based and end non-inclusive
        if chrom not in validChroms:
            print(f"Skipping invalid chromosome: {chrom}")
            continue
    
        regionStart, regionEnd = int(start), int(end)
    
        if regionStart < 1:
            print(f"Skipping region with invalid start: {regionStart}")
            continue
    
        #print(f"Processing region {idx + 1}/{total_regions}: {cid} ({chrom}:{start}-{end})")
        
        posRange = defaultdict(lambda: [0, 0])
        filteredReads = Intersecter()
        
        step_process += 1
        
        for bamfile in args:
            if options.verbose:
                sys.stderr.write(f"Reading {bamfile}\n")
            bamfile = bamfile.strip("""\'""")
            # Check if BAM and BAI files exist
            if os.path.exists(bamfile) and (os.path.exists(bamfile.replace(".bam", ".bai")) or os.path.exists(bamfile + ".bai")):
                #print(f"Processing BAM file: {bamfile}")
                input_file = pysam.Samfile(bamfile, "rb")
                prefix = ""
                # Determine if chromosome names have "chr" prefix
                for tchrom in input_file.references:
                    if tchrom.startswith("chr"):
                        prefix = "chr"
                        break
    
                for read in input_file.fetch(prefix + chrom, regionStart - protection - 1, regionEnd + protection + 1):
                    if read.is_duplicate or read.is_qcfail or read.is_unmapped:
                        continue
                    if isSoftClipped(read.cigar):
                        continue
    
                    if read.is_paired:
                        if read.mate_is_unmapped:
                            continue
                        if read.rnext != read.tid:
                            continue
                        if read.is_read1 or (read.is_read2 and read.pnext + read.qlen < regionStart - protection - 1):
                            if read.isize == 0:
                                continue
                            if options.downsample is not None and random.random() >= options.downsample:
                                continue
                            rstart = min(read.pos, read.pnext) + 1  # 1-based
                            lseq = abs(read.isize)
                            rend = rstart + lseq - 1  # end included
                            if minInsSize is not None and ((lseq < minInsSize) or (lseq > maxInsSize)):
                                continue
    
                            filteredReads.add_interval(Interval(rstart, rend))
                            for i in range(rstart, rend + 1):
                                if regionStart <= i <= regionEnd:
                                    posRange[i][0] += 1
                            if regionStart <= rstart <= regionEnd:
                                posRange[rstart][1] += 1
                            if regionStart <= rend <= regionEnd:
                                posRange[rend][1] += 1
                    else:
                        if options.downsample is not None and random.random() >= options.downsample:
                            continue
                        rstart = read.pos + 1  # 1-based
                        lseq = aln_length(read.cigar)
                        rend = rstart + lseq - 1  # end included
                        if minInsSize is not None and ((lseq < minInsSize) or (lseq > maxInsSize)):
                            continue
    
                        filteredReads.add_interval(Interval(rstart, rend))
                        for i in range(rstart, rend + 1):
                            if regionStart <= i <= regionEnd:
                                posRange[i][0] += 1
                        if ((options.merged or read.qname.startswith('M_')) or ((options.trimmed or read.qname.startswith('T_')) and read.qlen <= options.lengthSR - 10)):
                            if regionStart <= rstart <= regionEnd:
                                posRange[rstart][1] += 1
                            if regionStart <= rend <= regionEnd:
                                posRange[rend][1] += 1
                        elif read.is_reverse:
                            if regionStart <= rend <= regionEnd:
                                posRange[rend][1] += 1
                        else:
                            if regionStart <= rstart <= regionEnd:
                                posRange[rstart][1] += 1
    
            #print(f"Gone through Bam file: {bamfile}.")
        
        
        filename = options.outfile % cid
        # changed mode from "w" to "wt"
        outfile = gzip.open(filename, 'wt')
        cov_sites = 0
        outLines = []
        for pos in range(regionStart, regionEnd + 1):
            rstart, rend = pos - protection, pos + protection
            centerCount = 0
            for read in filteredReads.find(rstart, rend):
                center = (read.start + read.end) / 2
                if rstart <= center <= rend:
                    centerCount += 1
            covCount, startCount = posRange[pos]
            cov_sites += covCount
            outLines.append("%s\t%d\t%d\t%d\t%d\n" % (chrom, pos, covCount, startCount, centerCount))
    
        if strand == "-":
            outLines = outLines[::-1]
        for line in outLines:
            outfile.write(line)
        outfile.close()
    
        if cov_sites == 0 and not options.empty:
            os.remove(filename)
            # Print progress for each processed region
            #print(f"removed outfile: {cid}, {idx + 1}/{total_regions}")
            total_files_removed +=1
            
        if step_process == 1000:
            print(f"{idx + 1}/{total_regions} files done.")
            print(f"{total_files_removed} files have been removed.")
            step_process = 0
    
print("All regions processed.")