# Modified GraphST:

This repository is based on the original [GraphST](https://github.com/JinmiaoChenLab/GraphST) implementation by **Jinmiao Chen**, with modifications to extend its functionality.  

The modifications mainly adapt the workflow for the **DLPFC 12 slices** algorithm comparison.  

---

## Modifications in this repository  

To simplify the comparison across the **12 DLPFC slices**, this repository includes an automation script: [`run.py`](run.py).  
This script executes the workflow step by step of **GraphST**:

- Added **`run.py`** to automate the analysis of multiple `.h5ad` slices.  
- The script iterates through all input files in a dataset folder and runs the GraphST pipeline without manual per-slice calls.  
- Output folders are automatically created for each slice.  
- Logging and clustering results are stored in slice-specific directories.  

### Usage
```bash
# Simply run and add the path to the .h5ad files for each sample:
python run.py
```
# Code Pull from Spatially informed clustering, integration, and deconvolution of spatial transcriptomics with GraphST

[![DOI](https://zenodo.org/badge/494373596.svg)](https://zenodo.org/badge/latestdoi/494373596)

![](https://github.com/JinmiaoChenLab/GraphST/blob/main/GraphST.jpg)

## Overview
GraphST is a versatile graph self-supervised contrastive learning model that incorporates spatial location information and gene expression profiles to accomplish three key tasks, spatial clustering, spatial transcriptomics (ST) data integration, and single-cell RNA-seq (scRNA-seq) transfer onto ST. GraphST combines graph neural networks (GNNs) with self-supervised contrastive learning to learn spot representations in the ST data by modeling gene expressions and spatial locaiton information. After the representation learning, the non-spatial alignment algorithm is used to cluster the spots into different spatial domains. Each cluster is regarded as a spatial domain, containing spots with similar gene expression profiles and spatially proximate. GraphST can jointly analyze multiple ST samples while correcting batch effects, which is achieved by smoothing features between spatially adjacent spots across samples. For the scRNA-seq transfer onto ST data, a mapping matrix is trained via an augmentation-free contrastive learning mechanism, where the similarity of spatially adjacent spots are maximized while those of spatially non-adjacent spots are minimized. With the learned mapping matrix, arbitrary cell attributes (e.g., cell type and sample type) can be flexibly projected onto spatial space.   

## Requirements
You'll need to install the following packages in order to run the codes.
* python==3.8
* torch>=1.8.0
* cudnn>=10.2
* numpy==1.22.3
* scanpy==1.9.1
* anndata==0.8.0
* rpy2==3.4.1
* pandas==1.4.2
* scipy==1.8.1
* scikit-learn==1.1.1
* tqdm==4.64.0
* matplotlib==3.4.2
* R==4.0.3

## Tutorial
For the step-by-step tutorial, please refer to:
[https://deepst-tutorials.readthedocs.io/en/latest/](https://deepst-tutorials.readthedocs.io/en/latest/)

## Citation
Long et al. Spatially informed clustering, integration, and deconvolution of spatial transcriptomics with GraphST. _**Nature Communications**_. 14(1), 1155 (2023). 
