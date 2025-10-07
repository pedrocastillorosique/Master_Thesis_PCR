# Implemented stLearn:

This script runs a **minimal stLearn pipeline** to (i) tile the H&E image, (ii) extract **morphology features** with a CNN, (iii) perform **stSME** (spatial–morphology–expression) normalization, and (iv) **cluster** Visium spots.

- **Load AnnData & Spatial Metadata**
  - Reads each `.h5ad`, rescales `obsm["spatial"]` to **hi-res pixel space** using `tissue_hires_scalef`
  - Exposes `imagecol` / `imagerow` in `.obs` for stLear
- **Expression Pre-processing**
  - Library-size normalization → log transform (with a permissive gene filter)
- **Image Tiling (hi-res)**
  - Crops H&E into tiles (`crop_size=50` by default) for CNN feature extraction
- **Morphology Feature Extraction**
  - Runs stLearn’s pretrained CNN to populate `obsm["X_morphology"]` (also saved to `features.npy`)
- **PCA (Expression)**
  - Optional dimensionality reduction on expression for diagnostics
- **stSME Normalization + Clustering**
  - Applies **stSME** on raw counts to obtain `obsm["raw_SME_normalized"]`
  - Replaces `adata.X` with the stSME matrix, scales, runs PCA, builds neighbors, and clusters (Leiden)
- **Save Outputs**
  - Writes `results.h5ad` per sample with stSME matrix and cluster labels

### Usage
```bash
# Simply run and add the path to the .h5ad files for each sample:
python run.py #predefine the datapath
```
# Code Pull from stLearn - A Spatial Transcriptomics Toolkit

<p align="center">
  <img src="https://i.imgur.com/yfXlCYO.png"
    alt="deepreg_logo" title="DeepReg" width="300"/>
</p>

<table align="center">
  <tr>
    <td>
      <b>Package</b>
    </td>
    <td>
      <a href="https://pypi.python.org/pypi/stlearn/">
      <img src="https://img.shields.io/pypi/v/stlearn.svg" alt="PyPI Version">
      </a>
      <a href="https://pepy.tech/project/stlearn">
      <img src="https://static.pepy.tech/personalized-badge/stlearn?period=total&units=international_system&left_color=grey&right_color=orange&left_text=Downloads"
        alt="PyPI downloads">
      </a>
    </td>
  </tr>
  <tr>
    <td>
      <b>Documentation</b>
    </td>
    <td>
      <a href="https://stlearn.readthedocs.io/en/latest/">
      <img src="https://readthedocs.org/projects/stlearn/badge/?version=latest" alt="Documentation Status">
      </a>
    </td>
  </tr>
  <tr>
    <td>
     <b>Paper</b>
    </td>
    <td>
      <a href="https://doi.org/10.1038/s41467-023-43120-6"><img src="https://zenodo.org/badge/DOI/10.1038/s41467-023-43120-6.svg"
        alt="DOI"></a>
    </td>
  </tr>
  <tr>
    <td>
      <b>License</b>
    </td>
    <td>
      <a href="https://github.com/BiomedicalMachineLearning/stLearn/blob/master/LICENSE"><img src="https://img.shields.io/badge/License-BSD-blue.svg"
        alt="LICENSE"></a>
    </td>
  </tr>
</table>


# stLearn - A downstream analysis toolkit for Spatial Transcriptomic data

**stLearn** is designed to comprehensively analyse Spatial Transcriptomics (ST) data to investigate complex biological processes within an undissociated tissue. ST is emerging as the “next generation” of single-cell RNA sequencing because it adds spatial and morphological context to the transcriptional profile of cells in an intact tissue section. However, existing ST analysis methods typically use the captured spatial and/or morphological data as a visualisation tool rather than as informative features for model development. We have developed an analysis method that exploits all three data types: Spatial distance, tissue Morphology, and gene Expression measurements (SME) from ST data. This combinatorial approach allows us to more accurately model underlying tissue biology, and allows researchers to address key questions in three major research areas: cell type identification, spatial trajectory reconstruction, and the study of cell-cell interactions within an undissociated tissue sample.

---

## Getting Started

- [Documentation and Tutorials](https://stlearn.readthedocs.io/en/latest/)

## Citing stLearn

If you have used stLearn in your research, please consider citing us:

> Pham, Duy, et al. "Robust mapping of spatiotemporal trajectories and cell–cell interactions in healthy and diseased tissues."
> Nature Communications 14.1 (2023): 7739.
> [https://doi.org/10.1101/2020.05.31.125658](https://doi.org/10.1038/s41467-023-43120-6)
