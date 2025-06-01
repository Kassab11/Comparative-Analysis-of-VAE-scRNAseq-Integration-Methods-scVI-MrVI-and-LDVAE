import pandas as pd
import scanpy as sc
import scib
from scanpy.external import tl as external_tl
import warnings
warnings.filterwarnings("ignore")

def evaluate_model(model, adata_raw, adata_filtered, batch_key, label_key, use_rep):
    
    latent = model.get_latent_representation()
    adata_filtered.obsm[use_rep] = latent
    
    results = {}
    
    
    try: 
        results["Batch ASW"] = scib.me.silhouette_batch(adata_filtered, batch_key=batch_key, label_key=label_key, embed=use_rep, scale=True)
    except Exception as e: results["Batch ASW"] = str(e)
    
    try: 
        results["PCR Batch"] = scib.me.pcr_comparison(adata_raw, adata_filtered, covariate=batch_key, embed=use_rep)
    except Exception as e: results["PCR Batch"] = str(e)
    
    try: 
        results["iLISI"] = scib.me.ilisi_graph(adata_filtered, batch_key=batch_key, type_="embed", use_rep=use_rep)
    except Exception as e: results["iLISI"] = str(e)
    
    try:
        adata_filtered2 = adata_filtered.copy()
        adata_filtered2.obs[label_key] = adata_filtered2.obs[label_key].astype("category")
        results["Graph Connectivity"] = scib.me.graph_connectivity(adata_filtered2, label_key=label_key)
    except Exception as e: results["Graph Connectivity"] = str(e)
    
    try:
        scib.me.cluster_optimal_resolution(adata_filtered, cluster_key="leiden", label_key=label_key)
        results["NMI"] = scib.me.nmi(adata_filtered, cluster_key="leiden", label_key=label_key)
    except Exception as e: results["NMI"] = str(e)
    
    try:
        results["ARI"] = scib.me.ari(adata_filtered, cluster_key="leiden", label_key=label_key)
    except Exception as e: results["ARI"] = str(e)
    
    try:
        results["Label ASW"] = scib.me.silhouette(adata_filtered, label_key=label_key, embed=use_rep)
    except Exception as e: results["label ASW"] = str(e)
    
    try:
        results["Isolated Label F1"] = scib.me.isolated_labels_f1(adata_filtered, batch_key=batch_key, embed=None, label_key=label_key)
    except Exception as e: result["Isolated Label F1"] = str(e)
    
    try:
        results["Isolated Label ASW"] = scib.me.isolated_labels_asw(adata_filtered, batch_key=batch_key, label_key=label_key, embed=use_rep)
    except Exception as e: results["Isolated Label ASW"] = str(e)
    
    try:
        results["cLISI"] = scib.me.clisi_graph(adata_filtered, label_key=label_key, type_="embed", use_rep=use_rep)
    except Exception as e: results["cLISI"] = str(e)
    
    try: 
        results["Trajectory Conservation"] = scib.me.trajectory_conservation(adata_raw, adata_filtered, label_key=label_key)
    except Exception as e: results["Trajectory Conservation"] = str(e)
        
    return pd.DataFrame.from_dict(results, orient="index", columns=["score"])
    
