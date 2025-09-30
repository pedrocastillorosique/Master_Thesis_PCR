"""
@author: Pedro Castillo
Automates the comparison of 12 DLPFC slices using stKeep.
"""

import numpy as np
import os
import scanpy as sc
import sys
import pandas as pd
import importlib
from pathlib import Path 

# Add the parent directory to Python path to allow importing stKeep as a local module
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import stKeep

if __name__ == '__main__':

    # Load default parameter settings from stKeep
    parser  =  stKeep.parameter_setting()
    args    =  parser.parse_args()

    # Define the input folder containing all DLPFC .h5ad slices
    args.inputPath = ".../Dataset/ARTICULOS/DLPFC/"
    
    # Choose histology feature extraction model (alternative to SimCLR)
    Histmodel = "ResNet50"
    
    # Loop over all files in the dataset folder
    for i in os.listdir(args.inputPath):
        if 'h5ad' in i:   # Only process AnnData objects (.h5ad files)
            
            # Create an output folder for each slice inside stKeep/{slice_name}/stKeep/
            name = i.split(".")[0]  # Slice identifier (e.g., 151507)
            outPath = args.inputPath + f'stKeep/{name}/stKeep/'
            Path(outPath).mkdir(parents=True, exist_ok=True)

            # Load the AnnData object
            adata = sc.read_h5ad(args.inputPath + i)

            # Extract the "Region" annotations if available, otherwise assign "Cluster1"
            anno = pd.DataFrame(adata.obs['Region'])
            anno['Region'] = 'Cluster1'   # Default label if no annotation
            anno.to_csv(outPath + 'Image_cell_segmentation.txt', sep='\t', index=True) 

            # Determine number of clusters (ignore NaN values if present)
            if adata.obs['Region'].isna().sum() > 0:
                ncluster = len(adata.obs['Region'].unique()) - 1
            else: 
                ncluster = len(adata.obs['Region'].unique())

            # Save the full input path for current slice
            inputPath = args.inputPath + i

            # Import and reload the analysis module for each slice
            # (ensures fresh execution without caching)
            import analysis
            importlib.reload(analysis)
