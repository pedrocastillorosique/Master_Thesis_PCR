# Deep learning models that integrate transcriptomic and spatial information enable efficient reconstruction and clonal TCR analysis of the tumor microenvironment

This repository contains the code used to validate the [Louvain](https://doi.org/10.1038/s41598-019-41695-z), [BayesSpace](https://doi.org/10.1038/s41587-021-00935-2), [conGI](https://doi.org/10.1093/bib/bbad048), [GraphST](https://doi.org/10.1038/s41467-023-36796-3), and [stKeep](https://doi.org/10.1038/s41467-024-49171-7) algorithms, and implements GraphST on four 10x Visium–profiled ovarian cancer samples from the [PITAGORAS project](https://cima.cun.es/investigacion/proyecto-pitagoras).  
Validation was performed using the ARI metric on clustering results from 12 slices of the [DLPFC dataset](https://github.com/LieberInstitute/HumanPilot).

## Objectives of this Master’s Thesis

- **OC1.** Conduct a literature review focused on the tumor microenvironment and algorithms available for spatial RNA-seq (spaRNA-seq) data analysis.  
  <div align="center">
    <img width="964" height="360" alt="Overview" src="https://github.com/user-attachments/assets/987422e2-43b6-4fb1-8e80-acf6caec7b4b" />
  </div>

- **OC2.** Compare spatial domain identification (clustering) results obtained with state-of-the-art algorithms based on classical methods, Machine Learning, and Deep Learning, using annotated spaRNA-seq datasets.  
  For this objective, each algorithm was implemented in its native language and environment. You can reproduce the runs by following the scripts provided in the corresponding folders.  
  <div align="center">
    <img width="964" height="500" alt="Algorithm comparison" src="https://github.com/user-attachments/assets/7e9749d3-51aa-400a-abd6-ba2da4341a2b" />
  </div>

- **OC3.** Apply the best-performing clustering algorithm from OC2 to ovarian tumor spaRNA-seq data, integrating genetic and spatial information on T-cell receptor (TCR) clones in tumor-infiltrating lymphocytes (TILs), to identify relevant patterns in cell distribution and gene expression.  
  <div align="center">
    <img width="963" height="587" alt="Application to ovarian cancer" src="https://github.com/user-attachments/assets/ded256a9-5817-420d-a5fb-8d601ed453a6" />
  </div>

    **Figure.** (a) H&E-stained ovarian tumor tissue. (b) Manual annotation into border, stroma, and tumor regions.  
      (c) Clustering results obtained with the selected algorithm. (d) 2D embedding (UMAP) visualization of spatial clusters.
