"""
stlearn_pipeline.py
-------------------
Pipeline to tile Visium images, extract morphology features with stLearn,
apply SME normalization, and cluster spots on DLPFC slices.

Requirements:
    Python 3.9+
    scanpy >= 1.9
    stlearn >= 0.4
    numpy

Input:
    A folder containing .h5ad AnnData files (each one a Visium slice).

Output (per slice):
    <OUTROOT>/<sample_id>/tiles/         # image tiles (hires)
    <OUTROOT>/<sample_id>/features.npy   # morphology features from stLearn
    <OUTROOT>/<sample_id>/results.h5ad   # SME-normalized AnnData with clustering
"""


from __future__ import annotations

import warnings
from pathlib import Path
import numpy as np
import scanpy as sc
import stlearn as st

# (Optional) Silence verbose warnings from backends
warnings.filterwarnings("ignore")

# --------- CONFIGURE PATHS ---------
DATAPATH = Path(r"")
OUTROOT = DATAPATH.parent / "stLearn"
OUTROOT.mkdir(exist_ok=True, parents=True)

for h5_file in DATAPATH.glob("*.h5ad"):
    sample_id = h5_file.stem
    outdir = OUTROOT / sample_id
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"[stlearn_pipeline] Processing: {sample_id}")

   # ----------------------------
    # 1) Load AnnData & coordinates
    # ----------------------------
    adata = sc.read_h5ad(h5_file)

    # Visium scale factor: convert spot coords to 'hires' pixel coordinates
    lib = next(iter(adata.uns["spatial"].keys()))
    sf = adata.uns["spatial"][lib]["scalefactors"]
    adata.obsm["spatial"] = adata.obsm["spatial"] * sf["tissue_hires_scalef"]

    # stLearn expects 'imagecol' and 'imagerow' in .obs
    adata.obs["imagecol"] = adata.obsm["spatial"][:, 0]
    adata.obs["imagerow"] = adata.obsm["spatial"][:, 1]

    # -----------------------------------------
    # 2) Expression preprocessing (counts -> log)
    # -----------------------------------------
    st.pp.filter_genes(adata, min_cells=1)     # keep any gene seen in ≥1 spot
    st.pp.normalize_total(adata)               # library-size normalization
    st.pp.log1p(adata)                         # log-transform

    # -----------------------------
    # 3) Image tiling (hires level)
    # -----------------------------
    tiles_dir = outdir / "tiles"
    tiles_dir.mkdir(exist_ok=True)
    st.pp.tiling(
        adata,
        out_path=str(tiles_dir),
        crop_size=50,           # tile size in pixels on the hires image
        use_quality="hires"
    )

    # --------------------------------------------------------
    # 4) Morphology feature extraction (DL on image tiles)
    #    Populates adata.obsm['X_morphology']
    # --------------------------------------------------------
    st.pp.extract_feature(adata)
    np.save(outdir / "features.npy", np.asarray(adata.obsm["X_morphology"]))

    # -----------------------------------------------------
    # 5) PCA (expression) — used later or for diagnostics
    # -----------------------------------------------------
    st.em.run_pca(adata, n_comps=50)

    # ------------------------------------------------------
    # 6) SME normalization + clustering on normalized data
    # ------------------------------------------------------
    adata_sme = adata.copy()
    # SME on raw counts (creates 'raw_SME_normalized' in .obsm)
    st.spatial.SME.SME_normalize(adata_sme, use_data="raw")
    adata_sme.X = adata_sme.obsm["raw_SME_normalized"]

    # Scale & reduce (PCA is required before neighbors on adata_sme!)
    st.pp.scale(adata_sme)
    st.em.run_pca(adata_sme, n_comps=50)

    # Build graph & cluster (Leiden)
    st.pp.neighbors(adata_sme, n_neighbors=7, use_rep="X_pca")
    st.tl.clustering.leiden(adata_sme, key_added="leiden_sme")

    # ----------------------------
    # 7) Save result AnnData
    # ----------------------------
    adata_sme.write_h5ad(outdir / "results.h5ad")

    print(f"[stlearn_pipeline] Done: {sample_id} → {outdir}")