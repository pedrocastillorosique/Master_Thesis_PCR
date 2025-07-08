import os
import torch
import scanpy as sc

from GraphST import GraphST

# Check cuda availability: 
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

# Determine number of clusters manually 
## it will be checked in anndata.obs['Region'] any annotation: 
n_clusters = 7
# Stablish path: 
file_fold = r"D:\Pitagoras\Spatia_seq\Dataset\PITAGORAS"
# Set radius to specify the number of neighbors considered during refinement
radius = 50
# Set clustering tool:
tool = 'mclust' # mclust, leiden, and louvain

for file_name in os.listdir(file_fold):

    if 'h5ad' in file_name:
        # reading data: 
        adata = sc.read_h5ad(os.path.join(file_fold, file_name))
        adata.var_names_make_unique()

        # define model
        model = GraphST.GraphST(adata, device=device)

        # train model
        adata = model.train()

        if adata.obs['Region'].isna().sum() > 0:
            n_clusters = len(adata.obs['Region'].unique())-1
        else: 
            n_clusters = len(adata.obs['Region'].unique())

        # clustering
        from GraphST.utils import clustering

        if tool == 'mclust':
            clustering(adata, n_clusters, radius=radius, method=tool, refinement=True) # For DLPFC dataset, we use optional refinement step.
        elif tool in ['leiden', 'louvain']:
            clustering(adata, n_clusters, radius=radius, method=tool, start=0.1, end=2.0, increment=0.01, refinement=False) 
        adata.obs[tool] = adata.obs['domain']

        # saving path:
        saving_path =  file_fold + f'/graphST/{tool}'
        os.makedirs(saving_path, exist_ok=True)
        # save anndata
        sc.write(saving_path + f'/{file_name}', adata)


    elif file_name in ["GIN-63-3T-4", "GIN-65-3T-2", "GIN-67-3T-6", "GIN-71-3T-1"]:
        # reading data: 
        adata = sc.read_visium(os.path.join(file_fold, file_name, 'outs'))
        adata.var_names_make_unique()

        # define model
        model = GraphST.GraphST(adata, device=device)

        # train model
        adata = model.train()

        if file_name == "GIN-63-3T-4":
            n_clusters = 8
        elif file_name == "GIN-65-3T-2":
            n_clusters = 5
        elif file_name == "GIN-67-3T-6":
            n_clusters = 7
        elif file_name == "GIN-71-3T-1":
            n_clusters = 8

        # clustering
        from GraphST.utils import clustering

        if tool == 'mclust':
            clustering(adata, n_clusters, radius=radius, method=tool, refinement=False) # For DLPFC dataset, we use optional refinement step.
        elif tool in ['leiden', 'louvain']:
            clustering(adata, n_clusters, radius=radius, method=tool, start=0.1, end=2.0, increment=0.01, refinement=False) 
        adata.obs[tool] = adata.obs['domain']

        # saving path:
        saving_path =  file_fold + f'/graphST/{tool}'
        os.makedirs(saving_path, exist_ok=True)
        # save anndata
        sc.write(saving_path + f'/{file_name}', adata)