#! /usr/bin/env python
# -*- coding:utf-8 -*-
import os
import pandas as pd
import numpy as np
import scanpy as sc
import anndata
from sklearn.metrics import adjusted_rand_score


def res_search_fixed_clus(adata, fixed_clus_count, increment=0.02):
    for res in sorted(list(np.arange(0.02, 2, increment)), reverse=True):
        sc.tl.leiden(adata, random_state=0, resolution=res)
        count_unique_leiden = len(pd.DataFrame(adata.obs['leiden']).leiden.unique())
        if count_unique_leiden == fixed_clus_count:
            return res

def mclust_R(x, n_clusters, model='EEE', random_seed=2025):
    """\
    Clustering using the mclust algorithm.
    The parameters are the same as those in the R package mclust.
    """
    import rpy2.robjects as robjects
    from rpy2.robjects.packages import importr
    from rpy2.robjects.conversion import localconverter
    from rpy2.robjects import numpy2ri
    
    # Set environment variables if needed
    # os.environ['R_HOME'] = r'D:\Programs\R-4.4.1'
    os.environ['R_USER'] = r'C:\Users\pcastillor\AppData\Local\anaconda3\envs\conGI\Lib\site-packages\rpy2'
    # Set random seed
    np.random.seed(random_seed)
    # Import R package
    mclust = importr("mclust")
    # Set R seed
    robjects.r['set.seed'](random_seed)
    # Convert NumPy array to R object safely
    with localconverter(robjects.default_converter + numpy2ri.converter):
        r_x = robjects.conversion.py2rpy(x)
    # Call Mclust
    res = mclust.Mclust(r_x, G=n_clusters, modelNames=model)
    # Extract cluster assignments (2nd-to-last element in result, typically)
    with localconverter(robjects.default_converter + numpy2ri.converter):
        mclust_res = robjects.conversion.rpy2py(res[-2]).astype(int) - 1

    return mclust_res


def eval_mclust_ari(labels, z, n_clusters):
    raw_preds = mclust_R(z, n_clusters)
    return  raw_preds
