# Calculate FFT WPS for WPS
# for all samples
# find body -mindepth 1 -maxdepth 1 -type d | sed 's|^body/||' parallel --progress --eta --tmpdir /mnt/DATA3/timo/main/tmp/ -j 30 -I{} "(cd /mnt/DATA3/timo/main/body/{}/counts && find . -name 'block_*.tsv.gz' -printf '%f\0' | xargs -0 -n 1000 Rscript /mnt/DATA3/timo/scripts/data_pre-processing_Snyder/2_fft_path_modified.R /mnt/DATA3/timo/main/body/{}/counts /mnt/DATA3/timo/main/body/{}/fft)"

"""
:Original Version

:Author: -
:Source: https://github.com/shendurelab/cfDNA/blob/master/expression/fft_path.R
:Date: *03.06.2014
"""

"""
:Modified Version

:Author: Timo Rieger
:Contact: trieger@student.ethz.ch
:Date: *24.09.2024
"""


library(data.table)  # Make sure to load the data.table package for fread() function

args <- commandArgs(trailingOnly = TRUE)                                                                                                                                                                                                     
ipath <- args[1]                                                                                                                                                                                                                             
opath <- args[2]                                                                                                                                                                                                                             
args <- args[c(-1,-2)]                                                                                                                                                                                                                       

#print(paste("Input path: ", ipath))

#print(paste("Output path: ", opath))

# for testing take smallest file: chr21, only 127M (chr1 has 670M)
#args<- "genome_wide_chr21.tsv.gz"
#print(length(args))

#file_counter <- 0
#total_files <- length(args)
#print(paste("Total files to process: ", total_files))

for (filename in args) {                                                                                                                                                                                                             
  #print(paste("Reading file:", paste(ipath, filename, sep="/")))
  
  # Read the gzipped TSV file
  dataA1 <- fread(paste(ipath, filename, sep="/"), header = FALSE)
  #dataA1 <- read.table(paste(ipath, filename, sep="/"), as.is = TRUE, header = FALSE)
  #dataA1 <- read.table(gzfile(paste(ipath, filename, sep="/")), as.is = TRUE, header = FALSE)
  #print(paste("Finished reading in file:", paste(ipath, filename, sep="/"))) 

  # Convert the third column to a time series
  tdataA1 <- as.ts(dataA1$V3)                                                                                                                                                                                                                
  #print("Time series conversion completed for column V3.")

  # Filtering operations
  tdataA1 <- filter(c(window(tdataA1, 1, 300), tdataA1), filter = 1/seq(5, 100, 4), method = "recursive")[-(1:300)]
  #print("Applied filtering to the time series for column V3.")
  
  tdataA1 <- tdataA1 - mean(tdataA1, trim = 0.1)                                                                                                                                                                                                  
  #print("Mean centering applied to time series for column V3.")
  
  # Spectral analysis
  resA1 <- spec.pgram(tdataA1, pad = 0.3, tap = 0.3, span = 2, plot = FALSE, detrend = TRUE, demean = TRUE)                                                                                                                                                        
  resA1$freq <- 1/resA1$freq
  #print("Spectral analysis completed for column V3.")
  
  # Create results table
  restbl <- data.frame(Freq = round(resA1$freq[resA1$freq >= 120 & resA1$freq <= 280]), 
                       Cov = resA1$spec[resA1$freq >= 120 & resA1$freq <= 280])
  
  # Process V4
  valsA1 <- dataA1$V4/dataA1$V3
  valsA1[is.na(valsA1)] <- 0                                                                                                                                                                                                                  
  tdataA1 <- as.ts(valsA1)                                                                                                                                                                                                                    
  tdataA1 <- filter(c(window(tdataA1, 1, 300), tdataA1), filter = 1/seq(5, 100, 4), method = "recursive")[-(1:300)]                                                                                                                                     
  tdataA1 <- tdataA1 - mean(tdataA1, trim = 0.1)                                                                                                                                                                                                  
  #print("Processed column V4: filtering and mean centering completed.")
  
  # Spectral analysis for V4
  resA1 <- spec.pgram(tdataA1, pad = 0.3, tap = 0.3, span = 2, plot = FALSE, detrend = TRUE, demean = TRUE)                                                                                                                                                        
  resA1$freq <- 1/resA1$freq
  restbl$Starts <- resA1$spec[resA1$freq >= 120 & resA1$freq <= 280]
  #print("Spectral analysis completed for column V4.")
  
  # Process V5
  tdataA1 <- as.ts(dataA1$V5)                                                                                                                                                                                                                
  tdataA1 <- filter(c(window(tdataA1, 1, 300), tdataA1), filter = 1/seq(5, 100, 4), method = "recursive")[-(1:300)]                                                                                                                                     
  tdataA1 <- tdataA1 - mean(tdataA1, trim = 0.1)                                                                                                                                                                                                  
  #print("Processed column V5: filtering and mean centering completed.")
  
  # Spectral analysis for V5
  resA1 <- spec.pgram(tdataA1, pad = 0.3, tap = 0.3, span = 2, plot = FALSE, detrend = TRUE, demean = TRUE)                                                                                                                                                        
  resA1$freq <- 1/resA1$freq
  restbl$WPS <- resA1$spec[resA1$freq >= 120 & resA1$freq <= 280]
  #print("Spectral analysis completed for column V5.")
  
  # Writing output
  gz1 <- gzfile(paste(opath, filename, sep="/"), "w")                                                                                                                                                                                          
  write.table(restbl, gz1, quote = FALSE, sep = "\t", row.names = FALSE)                                                                                                                                                                             
  close(gz1) 
  #print(paste("Finished writing output for:", filename))
    
  #file_counter <- file_counter + 1
   #   if (file_counter %% 50 == 0) {
    #  print(paste(file_counter, "of", total_files, "files done."))
      #}
}
#print("1000 files done.")