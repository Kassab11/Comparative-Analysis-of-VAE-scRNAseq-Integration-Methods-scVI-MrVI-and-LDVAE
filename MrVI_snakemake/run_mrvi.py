import scanpy as sc
import pandas as pd
import scvi
import matplotlib.pyplot as plt
import numpy as np
import time
import torch
import scib
import time
import torch
import os
import argparse
import pymde
import numpy as np
from scvi.external import MRVI
import scanpy.external as external
from evaluate import evaluate_model
from utils import filter_data
from utils import dataset_utils


def run_mrvi(input_file, output_dir, batch_key, sample_key, label_key, use_rep, encoder_n_hidden, n_latent, encoder_n_layers, n_latent_u, max_epochs, HVG):
    
    batch_key, sample_key, label_key, use_rep, name = dataset_utils(input_file, "MRVI")
    
    adata_filtered = filter_data(input_file)
    adata_raw = adata_filtered.copy()
    
    if HVG:
        sc.pp.highly_variable_genes(
        adata_filtered,
        n_top_genes=5000,
        subset=True,
        flavor="seurat_v3",
        batch_key=batch_key,
    )
    
    MRVI.setup_anndata(adata_filtered, batch_key=batch_key, sample_key=sample_key)
        
    model = MRVI(
        adata_filtered,
        encoder_n_hidden=encoder_n_hidden,
        n_latent=n_latent,
        encoder_n_layers=encoder_n_layers,
        n_latent_u=n_latent_u
    )
    
    torch.cuda.empty_cache()
    torch.cuda.reset_peak_memory_stats()
    
    start_time = time.time()
    
    model.train(max_epochs=max_epochs, accelerator='gpu', devices=1)
    
    torch.cuda.synchronize()
    end_time = time.time()
    
    max_memory_reserved = torch.cuda.max_memory_reserved() / (1024 ** 3) 
    max_memory_allocated = torch.cuda.max_memory_allocated() / (1024 ** 3)
    
    time_req = end_time - start_time
    
    latent = model.get_latent_representation()
    adata_filtered.obsm[use_rep] = latent
    
    sc.tl.pca(adata_filtered)
    sc.pp.neighbors(adata_filtered, use_rep=use_rep)
    sc.tl.umap(adata_filtered)
    sc.tl.tsne(adata_filtered, use_rep=use_rep)


    adata_filtered.obsm["X_mde"] = pymde.preserve_neighbors(
        adata_filtered.obsm[use_rep],
        embedding_dim=2,
        constraint=pymde.Standardized(),
        repulsive_fraction=1.5,
        verbose=False,
        device="cuda"
    ).embed(verbose=False).cpu().numpy()
    
    results = evaluate_model(model, adata_raw, adata_filtered, batch_key, label_key, use_rep)
    
    
    #os.makedirs(output_dir, exist_ok=True)
    model.save(output_dir)
    
    with open(os.path.join(output_dir, 'metrics.txt'), "w") as f:
        f.write(f"model: \nencoder_n_hidden :{encoder_n_hidden},\n n_latent: {n_latent},\n encoder_n_layers: {encoder_n_layers},\n n_latent_u: {n_latent_u} \n")
        f.write(f"training time: {time_req: .2f} sec\n")
        f.write(f"memory reserved: {max_memory_reserved: .2f} GB\n")
        f.write(f"memory allocated: {max_memory_allocated: .2f} GB\n\n")
        f.write(results.to_string())
        
    output_h5ad = os.path.join(args.output_dir, f"{name}_w_embeds.h5ad")
    adata_filtered.write(output_h5ad)    
    
    
    print(f"model and metrics saved to {output_dir}")
    
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", required=True)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--batch_key", type=str, default=None, required=False)
    parser.add_argument("--sample_key", type=str, default=None, required=False)
    parser.add_argument("--label_key", type=str, default=None, required=False)
    parser.add_argument("--use_rep", type=str, default=None, required=False)
    parser.add_argument("--encoder_n_hidden", type=int, required=True)
    parser.add_argument("--n_latent", type=int, required=True)
    parser.add_argument("--encoder_n_layers", type=int, required=True)
    parser.add_argument("--n_latent_u", type=int, required=True)
    parser.add_argument("--max_epochs", type=int,  required=True)
    parser.add_argument("--HVG", default=False,  required=False)
    args = parser.parse_args()
    
    run_mrvi(
        args.input_file,
        args.output_dir,
        args.batch_key,
        args.sample_key,
        args.label_key,
        args.use_rep,
        args.encoder_n_hidden,
        args.n_latent,
        args.encoder_n_layers,
        args.n_latent_u,
        args.max_epochs,
        args.HVG
    )

    
    
