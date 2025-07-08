"""
@author: Pedro Castillo
"""

import numpy as np
import os
import scanpy as sc
import sys
import pandas as pd
import importlib

from pathlib import Path 

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import stKeep

if __name__ == '__main__':

    parser  =  stKeep.parameter_setting()
    args    =  parser.parse_args()

    args.inputPath = "D:/Pitagoras/Spatia_seq/Dataset/ARTICULOS/DLPFC/"
    Histmodel = "ResNet50"
    
    for i in os.listdir(args.inputPath):
        if 'h5ad' in i: 
            # Create saving folders
            name = i.split(".")[0]
            outPath        = args.inputPath + f'stKeep/{name}/stKeep/'
            Path(outPath).mkdir(parents=True, exist_ok=True)

            adata = sc.read_h5ad(args.inputPath + i)
            anno = pd.DataFrame(adata.obs['Region'])
            anno['Region'] = 'Cluster1'
            anno.to_csv(outPath + 'Image_cell_segmentation.txt', sep='\t', index=True) 

            if adata.obs['Region'].isna().sum() > 0:
                ncluster = len(adata.obs['Region'].unique())-1
            else: 
                ncluster = len(adata.obs['Region'].unique())

            inputPath = args.inputPath + i
            import analysis
            importlib.reload(analysis)