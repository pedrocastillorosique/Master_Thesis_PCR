import os
import scanpy as sc
import random
import numpy as np
import torch

from utils import get_predicted_results, parameter_setting
import pandas as pd
import warnings
warnings.filterwarnings("ignore")

parser = parameter_setting()
args = parser.parse_args(args=['--epochs', '10', '--dataset', 'DLPFC', '--gene_preprocess', 'hvg'])

def seed_torch(seed):
    random.seed(seed)
    os.environ['PYTHONHASHSEED'] = str(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.backends.cudnn.benchmark = False
    torch.backends.cudnn.deterministic = True


for data_name in os.listdir(args.path):
    seed_torch(2025)

    if 'h5ad' in data_name:

        args.name = data_name # file name with the extension
        name = args.name.split('.')[0] # file name 

        adata = sc.read_h5ad(os.path.join(args.path, args.name)) # read anndata
        # Load embeddings 
        xg = np.load(f'{args.path}/embeddings/{name}/{name}_xg.npy')
        xi = np.load(f'{args.path}/embeddings/{name}/{name}_xi.npy')
        z = xg + 0.1*xi # concate embeddings 
        # Predict
        ari, pred_label = get_predicted_results(args.dataset, args.name, args.path, z)
        print("Ari value : ", ari)
        #Saving path
        save_path = f'{args.path}/conGI'
        os.makedirs(save_path, exist_ok=True)
        # Save classification
        adata.obs['conGI'] = pred_label
        sc.write(save_path + f'/{name}.h5ad', adata)