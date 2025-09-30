# Modified conGI:

This repository is based on the original [conGI](https://github.com/biomed-AI/ConGI) implementation by **Yuansong Zeng**, with modifications to extend its functionality.  

The modifications mainly adapt the workflow for the **DLPFC 12 slices** algorithm comparison.  

---

## Modifications in this repository  

To simplify the comparison across the **12 DLPFC slices**, this repository includes `run.py` & `eval.py` to train and infer the emmbeding and apply clustering to the DLPFC. 
This script executes the workflow step by step, calling the core modules of **conGI** in sequence:  

1. **[`run.py`](conGI/run.py)**
   - Automatically generates the embedding for the genes `{sample_name}_xg.npy` and image crops `{sample_name}_xi.npy` used for clustering analysis.
3. **[`eval.py`](conGI/eval.py)**
   - Automatically calculates the clustering using Mclust R function and calculates ARI metric. Finally saves the clustering as `{sample_name}.h5ad`.

### Usage
```bash
# Simply run and add the path to the .h5ad files for each sample:
python run.py --path ../data_DLPFC/
# Now the embeddings should be in the same -path folder: 
python eval.py --path ../data_DLPFC/
```
# Code Pull from original alogorihm conGI:

## Overview
We propose a novel method ConGI to accurately decipher spatial domains by integrating gene expression and histopathologi-cal images, where the gene expression is adapted to image infor-mation through contrastive learning. We introduce three contrastive loss functions within and between modalities to learn the common semantic representations across all modalities while avoiding their meaningless modality-private noise information. The learned rep-resentations are then used for deciphering spatial domains through a clustering method. By comprehensive tests on tumor and normal spatial transcriptomics datasets, ConGI was shown to outperform existing methods in terms of spatial domain identification. More importantly, the learned representations from our model have also been used efficiently for various downstream tasks, including trajectory inference, clustering, and visualization.

![(Variational) gcn](framework.bmp)


## Requirements
Please ensure that all the libraries below are successfully installed:

- **torch 1.7.1**
- CUDA Version 10.2.89
- scanpy 1.8.1
- mclust








## Run ConGI on the example data.

```

python train.py --dataset SpatialLIBD  --name 151509 

```


## output

The clustering labels will be stored in the dir `output` /dataname_pred.csv. 


## Tutorial

We also provide a [Tutorial](https://github.com/biomed-AI/ConGI/blob/main/tutorial.ipynb) script for users. 



## Datasets

The spatial transcriptomics datasets that support the findings of this study are available here:
(1) human HER2-positive breast tumor ST data https://github.com/almaan/HER2st/. 
(2) The LIBD human dorsolateral prefrontal cortex (DLPFC) data was acquired with 10X Visium composed of spatial transcriptomics data acquired from twelve tissue slices (http://research.libd.org/spatialLIBD/).
(3) The mouse brain anterior section from 10X Visium (https://www.10xgenomics.com/resources/datasets). 
(4) the human epidermal growth factor receptor (HER) 2-amplified (HER+) invasive ductal carcinoma (IDC) (https://support.10xgenomics.com/spatial-gene-expression/datasets). 




## Citation

Please cite our paper:

```

@article{zengys,
  title={Deciphering Spatial Domains by Integrating Histopathological Image and Tran-scriptomics via Contrastive Learning },
  author={Yuansong Zeng1,#, Rui Yin3,#, Mai Luo1, Jianing Chen1, Zixiang Pan1, Yutong Lu1, Weijiang Yu1* , Yuedong Yang1,2*},
  journal={biorxiv},
  year={2022}
 publisher={Cold Spring Harbor Laboratory}
}

```
