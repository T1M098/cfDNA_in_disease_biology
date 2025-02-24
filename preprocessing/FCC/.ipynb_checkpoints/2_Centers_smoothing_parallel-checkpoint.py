import os
import argparse
import pandas as pd
import numpy as np
import scipy.ndimage
from whittaker_eilers import WhittakerSmoother
import gzip
import sys
import matplotlib.pyplot as plt
from multiprocessing import Pool

# Define smoothing function
def smooth_fragment_centers(fragment_center_array):
    whittaker_smoother = WhittakerSmoother(
        lmbda=1000, order=2, data_length=len(fragment_center_array))
    smoothed_fragment_centers = np.array(whittaker_smoother.smooth(fragment_center_array))
    sigma = 30  # Adjust sigma to control the width of the smoothing
    smoothed_fragment_centers = scipy.ndimage.gaussian_filter1d(smoothed_fragment_centers, sigma)
    return smoothed_fragment_centers


# Smooth center counts in each file
def process_files(sample):
    dir_path = f'/main/body/{sample}/counts/'
    outdir_path = f'/main/body/{sample}/smoothed/'   
    os.makedirs(outdir_path, exist_ok=True)

    for filename in os.listdir(dir_path):
        input_path = os.path.join(dir_path, filename)
        output_path = os.path.join(outdir_path, filename)
        # Define header names
        header = ["chr", "pos", "covCount", "startCount", "CenterCount"]
        # Read in File
        centers = pd.read_csv(input_path, compression="gzip", sep="\t", header=None, names=header)
        # Extract relevant column, containing values for smoothing
        fragment_centers = centers["CenterCount"].values
        smoothed_centers = smooth_fragment_centers(fragment_centers)
        # To avoid negative values after smoothing
        smoothed_centers = [max(0, x) for x in smoothed_centers]
        centers["smoothedCenterCount"] = smoothed_centers
        # Write output file
        centers.to_csv(output_path, sep='\t', index=False, header=False)
        
def main():
    # Create an argument parser
    parser = argparse.ArgumentParser(description="Process cfDNA data for a given sample.")

    # Add an argument for the sample
    parser.add_argument("sample", type=str, help="The sample name to process")

    # Parse the arguments
    args = parser.parse_args()

    # Process the given sample
    process_files(args.sample)
    print(args)

if __name__ == "__main__":
    main()

