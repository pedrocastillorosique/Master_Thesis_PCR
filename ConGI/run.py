import argparse
import numpy as np
import random
import torch
import os
import warnings

from torch.utils.data import DataLoader
from dataset import Dataset
from model import SpaCLR, TrainerSpaCLR
from utils import parameter_setting

warnings.filterwarnings("ignore")

def seed_torch(seed):
    random.seed(seed)
    os.environ['PYTHONHASHSEED'] = str(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.backends.cudnn.benchmark = False
    torch.backends.cudnn.deterministic = True

def fit_SpaCLR(trainloader, epochs):
    "Function to avoid freeze_support() of the threads"
    trainer.fit(trainloader, epochs)


if __name__ == "__main__":
    torch.multiprocessing.freeze_support() # Avoid freezing of threads
    seed_torch(2025)   

    # Get parameters settings and adjust them
    parser = parameter_setting()
    args = parser.parse_args(args=['--epochs', '10', '--dataset', 'DLPFC', '--gene_preprocess', 'hvg'])

    # Iterate through each dataset
    for data_name in os.listdir(args.path):
        
        args.name = data_name
        name = args.name.split('.')[0]
        path = args.path
        gene_preprocess = args.gene_preprocess
        n_gene = args.n_gene
        last_dim = args.last_dim
        gene_dims=[n_gene, 2*last_dim]
        image_dims=[n_gene]
        lr = args.lr
        p_drop = args.p_drop
        batch_size = args.batch_size
        dataset = args.dataset
        epochs = args.epochs
        img_size = args.img_size
        device = args.device
        log_name = args.log_name
        num_workers = args.num_workers
        prob_mask = args.prob_mask
        pct_mask = args.pct_mask
        prob_noise = args.prob_noise
        pct_noise = args.pct_noise
        sigma_noise = args.sigma_noise
        prob_swap = args.prob_swap
        pct_swap = args.pct_swap

        # dataset
        trainset = Dataset(dataset, path, data_name, gene_preprocess=gene_preprocess, n_genes=n_gene,
                        prob_mask=prob_mask, pct_mask=pct_mask, prob_noise=prob_noise, pct_noise=pct_noise, sigma_noise=sigma_noise,
                        prob_swap=prob_swap, pct_swap=pct_swap, img_size=img_size, train=True)
        trainloader = DataLoader(trainset, batch_size=batch_size, shuffle=True, num_workers=num_workers, pin_memory=True)

        testset = Dataset(dataset, path, data_name, gene_preprocess=gene_preprocess, n_genes=n_gene,
                        prob_mask=prob_mask, pct_mask=pct_mask, prob_noise=prob_noise, pct_noise=pct_noise, sigma_noise=sigma_noise,
                        prob_swap=prob_swap, pct_swap=pct_swap, img_size=img_size, train=False)
        testloader = DataLoader(testset, batch_size=batch_size, shuffle=False, num_workers=num_workers, pin_memory=True)

        # model
        network = SpaCLR(gene_dims=gene_dims, image_dims=image_dims, p_drop=p_drop, n_pos=trainset.n_pos, backbone='resnet', projection_dims=[last_dim, last_dim])
        optimizer = torch.optim.AdamW(network.parameters(), lr=lr)
        
        # log
        save_name = f'{name}_{args.w_g2i}_{args.w_g2g}_{args.w_i2i}'
        log_dir = f'{args.path}/log/{name}/{save_name}'

        # train
        trainer = TrainerSpaCLR(args, 7, network, optimizer, log_dir, device=device)
        fit_SpaCLR(trainloader, epochs)

        # inference
        xg, xi, _ = trainer.valid(testloader)

        # create folders
        os.makedirs(f'{args.path}/embeddings/{name}', exist_ok=True)

        # save embeddings
        np.save(f'{args.path}/embeddings/{name}/{name}_xg.npy', xg)
        np.save(f'{args.path}/embeddings/{name}/{name}_xi.npy', xi)