import os
import scanpy as sc
import numpy as np


def dataset_utils(input_file, model):
    
    file_name = os.path.basename(input_file)
    name = os.path.splitext(file_name)[0]
    
    if name == "immune_cell_hum":
        batch_key = "batch"
        label_key = "cell_type"
        if model == "MRVI":
            use_rep = "X_mrVI"
            sample_key = "sample_ID"
        
    
    elif name == "imyoo_remission_biome_pilot_2024-05-01":
        batch_key = "library"
        label_key = "cell_type"
        if model == "MRVI":
            use_rep = "X_mrVI"
            sample_key = "participant_id"
    
    elif name == "imyoo_capillary_blood_samples_76535_pbmcs":
        batch_key = "run_lane_batch"
        label_key = "cell_type_level_3"
        if model == "MRVI":
            use_rep = "X_mrVI"
            sample_key = "Participant IDs"
            
    else:
        print("please choose an available option")
        return None
    
    return batch_key, sample_key, label_key, use_rep, name
    
    
def filter_data(input_file):

    adata = sc.read_h5ad(input_file)
    
    file_name = os.path.basename(input_file)
    name = os.path.splitext(file_name)[0]
    
    if name == "immune_cell_hum":
    
        del adata.obsm["X_pca"] 
        
        adata.var_names_make_unique()
        adata.X = adata.layers["counts"] 
        adata.var["mt"] = adata.var_names.str.startswith("MT-")
        adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
        adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")
        sc.pp.calculate_qc_metrics(
            adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True
        )

        filtered_batches = []
        batches = adata.obs["batch"].unique()
        
        for batch in batches:
            batch_data = adata[adata.obs["batch"] == batch].copy()

    
            min_genes_batch = batch_data.obs["n_genes_by_counts"].quantile(0.01)  
            max_counts_batch = batch_data.obs["total_counts"].quantile(0.99) 

            print(f"Batch {batch}: Removing cells with <{int(min_genes_batch)} genes and >{int(max_counts_batch)} counts")

    
            batch_data = batch_data[(batch_data.obs["n_genes_by_counts"] >= min_genes_batch) & 
                            (batch_data.obs["total_counts"] <= max_counts_batch)]
    
            filtered_batches.append(batch_data)
            
        adata_filtered = sc.concat(filtered_batches, join="inner")
        sc.pp.filter_genes(adata_filtered, min_cells=3)
        sc.pp.scrublet(adata_filtered, batch_key="batch")
        adata_filtered.X = np.rint(adata_filtered.X).astype(np.int64)
        
    elif name == "imyoo_remission_biome_pilot_2024-05-01":
    
        del adata.obsm["X_scVI"] 
        del adata.obsm["X_tsne"] 
        del adata.obsm["X_umap"] 
        
        adata = adata[adata.obs["cell_type"] != "Unlabeled"].copy()
        
        adata.var_names_make_unique()
        
        adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
        adata.var["ribo"] = adata.var_names.str.upper().str.startswith(("RPS", "RPL"))
        adata.var["hb"] = adata.var_names.str.upper().str.contains("^HB[^(P)]")
        
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True)

        filtered_batches = []
        batches = adata.obs["library"].unique()

        for batch in batches:
            batch_data = adata[adata.obs["library"] == batch].copy()

            min_genes = batch_data.obs["n_genes_by_counts"].quantile(0.01)  
            max_counts = batch_data.obs["total_counts"].quantile(0.99)

            print(f"Library {batch}: Removing cells with <{int(min_genes)} genes and >{int(max_counts)} counts")

            batch_data = batch_data[
                (batch_data.obs["n_genes_by_counts"] >= min_genes) &
                (batch_data.obs["total_counts"] <= max_counts)
            ]   
            filtered_batches.append(batch_data)
            
        adata_filtered = sc.concat(filtered_batches, join="inner")
        sc.pp.filter_genes(adata_filtered, min_cells=3)
        sc.pp.scrublet(adata_filtered, batch_key="library")
        
    elif name == "imyoo_capillary_blood_samples_76535_pbmcs":

        del adata.obsm["X_mde"] 
        del adata.obsm["X_scvi"] 
        del adata.obsm["X_tsne"] 
        del adata.obsm["X_umap"] 
        
        non_pbmc_types = ["Mast Cells", "asDC", "cDC2", "cDC3", "pDC", "tumorDC"]
        
        adata = adata[
            adata.obs["cell_type_level_3"].notna() & 
            (~adata.obs["cell_type_level_3"].isin(non_pbmc_types))
        ].copy()
        
        adata.var_names_make_unique()

        adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
        adata.var["ribo"] = adata.var_names.str.upper().str.startswith(("RPS", "RPL"))
        adata.var["hb"] = adata.var_names.str.upper().str.contains("^HB[^(P)]")
        
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True)
        
        filtered_batches = []
        batches = adata.obs["run_lane_batch"].dropna().unique() 
        
        for batch in batches:
            batch_data = adata[adata.obs["run_lane_batch"] == batch].copy()
            valid_cells = batch_data.obs["n_genes_by_counts"].notna() & batch_data.obs["total_counts"].notna()
            batch_data = batch_data[valid_cells].copy()
            
            if batch_data.n_obs < 10:
                print(f"Skipping Library {batch}: too few valid cells after NaN removal.")
                continue
                
            min_genes = batch_data.obs["n_genes_by_counts"].quantile(0.01)
            max_counts = batch_data.obs["total_counts"].quantile(0.99)
            
            print(f"Library {batch}: Removing cells with <{int(min_genes)} genes and >{int(max_counts)} counts")
            
            batch_data = batch_data[
                (batch_data.obs["n_genes_by_counts"] >= min_genes) &
                (batch_data.obs["total_counts"] <= max_counts)
            ].copy()
            
            if batch_data.n_obs > 0:
                filtered_batches.append(batch_data)
            else:
                print(f"Skipping Library {batch}: no cells left after QC filtering.")
                
        adata_filtered = sc.concat(filtered_batches, join="inner")
        sc.pp.filter_genes(adata_filtered, min_cells=3)
        sc.pp.scrublet(adata_filtered, batch_key="run_lane_batch")

        
    return adata_filtered
