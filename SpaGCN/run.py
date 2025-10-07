# =========================
# Batch SpaGCN 
# =========================
# This script loops through a directory with .h5ad files (AnnData/Scanpy),
# computes adjacency matrices, trains SpaGCN, refines clusters, and saves results
# (H5AD + PNGs). 

import matplotlib.pyplot as plt
import os
import numpy as np
import random
import scanpy as sc
from SpaGCN_package import SpaGCN as spg
import torch
import warnings
warnings.filterwarnings("ignore")

# === Path containing input .h5ad files ===
datapath = r''

# Loop over each file in the input directory
for file in os.listdir(datapath):

    # Full path to the current file
    data_path_slice = datapath + f'/{file}'

    # Extract filename without extension
    file = file.split('.')[0]

    # Create an output directory for this sample
    a = datapath.split('01_Original')[0]
    results_path_slice = rf'{a}SpaGCN\{file}'
    os.makedirs(results_path_slice, exist_ok=True)

    # Load AnnData object
    adata = sc.read_h5ad(data_path_slice)

    # Extract histology image (Visium data stores images in adata.uns)
    img = adata.uns['spatial'][file]['images']['hires']
    sf = adata.uns['spatial'][file]['scalefactors']

    # === Spatial coordinates ===
    # Array coordinates (grid)
    x_array = adata.obs["array_col"].tolist()
    y_array = adata.obs["array_row"].tolist()

    # Pixel coordinates from obsm['spatial']
    adata.obs['x_pixel'] = (adata.obsm['spatial'] * sf['tissue_hires_scalef'])[:, 1]
    adata.obs['y_pixel'] = (adata.obsm['spatial'] * sf['tissue_hires_scalef'])[:, 0]
    x_pixel = adata.obs["x_pixel"].tolist()
    y_pixel = adata.obs["y_pixel"].tolist()

    # Overlay coordinates on the histology image for QC
    img_new = img.copy()
    for i in range(len(x_pixel)):
        x = x_pixel[i]
        y = y_pixel[i]
        img_new[int(x-20):int(x+20), int(y-20):int(y+20), :] = 0
    plt.imsave(f'{results_path_slice}/{file}_map.png', img_new)

    # === Adjacency matrix ===
    # Example if histology is available:
    adj = spg.calculate_adj_matrix(x=x_array, y=y_array, x_pixel=x_pixel, y_pixel=y_pixel, image=img, beta=49, alpha=1, histology=True, path=results_path_slice)
    # Here: no histology for adjacency
    # adj = spg.calculate_adj_matrix(x=x_pixel, y=y_pixel, histology=False)
    np.savetxt(f'{results_path_slice}/adj.csv', adj, delimiter=',')

    # Load adjacency back (ensures compatibility downstream)
    adj = np.loadtxt(f'{results_path_slice}/adj.csv', delimiter=',')

    # === Gene filtering and normalization ===
    adata.var_names_make_unique()
    spg.prefilter_genes(adata, min_cells=3)  # remove genes expressed in too few cells
    spg.prefilter_specialgenes(adata)        # remove unwanted genes
    sc.pp.normalize_per_cell(adata)          # per-cell normalization
    sc.pp.log1p(adata)                       # log-transform

    # === Hyperparameter search ===
    p = 0.5
    # Search for "l" given p
    l = spg.search_l(p, adj, start=0.01, end=1000, tol=0.01, max_run=100)

    # Number of clusters: if "Region" column exists, count categories
    labels = set(adata.obs['Region'])
    if np.nan in labels:
        n_clusters = len(labels) - 1
    else: 
        n_clusters = len(labels)

    # === Reproducibility: fix seeds ===
    r_seed = t_seed = n_seed = 2025
    random.seed(r_seed)
    torch.manual_seed(t_seed)
    np.random.seed(n_seed)

    # Search for optimal resolution
    res = spg.search_res(
        adata, adj, l, n_clusters,
        start=0.7, step=0.1, tol=5e-3, lr=0.05,
        max_epochs=20, r_seed=r_seed, t_seed=t_seed, n_seed=n_seed
    )

    # === Train SpaGCN ===
    clf = spg.SpaGCN()
    clf.set_l(l)
    clf.train(
        adata, adj, init_spa=True, init="louvain", res=res,
        tol=5e-3, lr=0.05, max_epochs=200
    )
    y_pred, prob = clf.predict()
    _, adata.obs["pred"] = np.unique(y_pred, return_inverse = True)

    # === Cluster refinement (optional) ===
    # Use "hexagon" for Visium data, "square" for ST data.
    adj_2d = spg.calculate_adj_matrix(x=x_array, y=y_array, histology=False)

    refined_pred = spg.refine(
        sample_id=adata.obs.index.tolist(),
        pred=adata.obs["pred"].tolist(),
        dis=adj_2d, shape="hexagon"
    )
    adata.obs["refined_pred"] = refined_pred
    adata.obs["refined_pred"] = adata.obs["refined_pred"].astype('category')

    # Save AnnData with results
    adata.write_h5ad(f"{results_path_slice}/results.h5ad")

    # Reload results
    adata = sc.read(f"{results_path_slice}/results.h5ad")

    # === Colors for plotting ===
    plot_color = [
        "#F56867","#FEB915","#C798EE","#59BE86","#7495D3","#D1D1D1","#6D1A9C",
        "#15821E","#3A84E6","#997273","#787878","#DB4C6C","#9E7A7A","#554236",
        "#AF5F3C","#93796C","#F9BD3F","#DAB370","#877F6C","#268785"
    ]
    
    # === Plot predicted spatial domains ===
    domains = "pred"
    num_celltype = len(adata.obs[domains].unique())
    adata.uns[domains + "_colors"] = list(plot_color[:num_celltype])
    ax = sc.pl.scatter(
        adata, alpha=1, x="y_pixel", y="x_pixel",
        color=domains, title=domains,
        show=False, size=100000/adata.shape[0]  # scaled dot size
    )
    ax.set_aspect('equal', 'box')
    ax.axes.invert_yaxis()
    plt.savefig(f"{results_path_slice}/pred.png", dpi=600)
    plt.close()

    # === Plot refined spatial domains ===
    domains = "refined_pred"
    num_celltype = len(adata.obs[domains].unique())
    adata.uns[domains + "_colors"] = list(plot_color[:num_celltype])
    ax = sc.pl.scatter(
        adata, alpha=1, x="y_pixel", y="x_pixel",
        color=domains, title=domains,
        show=False, size=100000/adata.shape[0]
    )
    ax.set_aspect('equal', 'box')
    ax.axes.invert_yaxis()
    plt.savefig(f"{results_path_slice}/refined_pred.png", dpi=600)
    plt.close()

    # === UMAP visualization (explicit PCA to avoid Scanpy warning) ===
    # Compute PCA first so neighbors uses X_pca instead of raw high-dimensional X
    sc.pp.pca(adata, n_comps=50, random_state=2025)     # control #PCs + reproducibility
    sc.pp.neighbors(adata, n_neighbors=100, use_rep="X_pca")
    sc.tl.umap(adata, random_state=2025)

    ax = sc.pl.umap(adata, color=domains, title=['SpaGCN'], show=False)
    plt.savefig(f"{results_path_slice}/umap_{file}.png", dpi=600)

    plt.close()
