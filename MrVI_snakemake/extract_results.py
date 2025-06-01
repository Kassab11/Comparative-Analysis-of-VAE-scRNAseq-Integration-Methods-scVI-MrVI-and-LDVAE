import os
import pandas as pd
import numpy as np

# === CONFIGURATION ===
root_dir = "/home/kassab/Desktop/CB803_project/MrVI_snakemake/results"
selected_datasets = ["imyoo_capillary_blood_samples_76535_pbmcs_full", "imyoo_capillary_blood_samples_76535_pbmcs_hvg"]
metric_filename = "metrics.txt"  # Be sure the filename is correct!
method = "MrVI"
output = "Embedding"
out_filename = "imyoo_capillary_blood_samples_76535_pbmcs"

metrics_keys = [
    "Batch ASW", "PCR Batch", "iLISI", "Graph Connectivity",
    "NMI", "ARI", "Label ASW", "Isolated Label F1",
    "Isolated Label ASW", "cLISI", "Trajectory Conservation"
]

columns = [
    "Method", "Hyperparameters", "Output", "Features", "Reference Name",
    "Training Time (s)", "GPU Memory Used (GB)", "Overall", "Overall Batch",
    "Batch ASW Score", "PCR Batch Score", "iLISI Score", "Graph Connectivity Score", "Overall BioConv",
    "NMI Score", "ARI Score", "Label ASW Score", "Isolated Label F1 Score",
    "Isolated Label ASW Score", "cLISI Score", "Trajectory Conservation Score"   
]

def safe_float(val):
    try:
        return float(val)
    except:
        return np.nan

rows = []

# === PROCESS EACH DATASET ===
for dataset in selected_datasets:
    dataset_path = os.path.join(root_dir, dataset)
    if not os.path.isdir(dataset_path):
        continue

    for subdir in os.listdir(dataset_path):
        metrics_path = os.path.join(dataset_path, subdir, metric_filename)
        if not os.path.isfile(metrics_path):
            continue

        features = "HVG" if dataset.endswith("hvg") else "Full"
        training_time = np.nan
        memory_used = np.nan
        metrics_values = {key: np.nan for key in metrics_keys}
        hyperparams = ""
        reference_name = ""

        with open(metrics_path, "r") as f:
            lines = f.readlines()
            params = {
                "encoder_n_hidden": None,
                "n_latent": None,
                "encoder_n_layers": None,
                "n_latent_u": None
            }

            for i, line in enumerate(lines):
                line = line.strip()

                if "encoder_n_hidden" in line:
                    for j in range(i, i + 4):
                        key_val = lines[j].strip().replace(",", "").split(":")
                        if len(key_val) == 2:
                            key = key_val[0].strip()
                            val = key_val[1].strip()
                            if key in params:
                                params[key] = val
                    hyperparams = f"n_hidden: {params['encoder_n_hidden']}, n_latent: {params['n_latent']}, n_layers: {params['encoder_n_layers']}, n_latent_u: {params['n_latent_u']}"
                    reference_name = f"h{params['encoder_n_hidden']}_l{params['n_latent']}_n{params['encoder_n_layers']}_u{params['n_latent_u']}"
                
                elif "training time" in line:
                    training_time = safe_float(line.split(":")[-1].replace("sec", "").strip())
                elif "memory reserved" in line:
                    memory_used = safe_float(line.split(":")[-1].replace("GB", "").strip())
                else:
                    for key in metrics_keys:
                        if line.startswith(key):
                            metrics_values[key] = safe_float(line.split()[-1])

        # Append only after parsing the whole file
        row = [
            method, hyperparams, output, features, reference_name,
            training_time, memory_used, "", "",
            metrics_values["Batch ASW"], metrics_values["PCR Batch"], metrics_values["iLISI"],
            metrics_values["Graph Connectivity"], "",
            metrics_values["NMI"], metrics_values["ARI"], metrics_values["Label ASW"],
            metrics_values["Isolated Label F1"], metrics_values["Isolated Label ASW"],
            metrics_values["cLISI"], metrics_values["Trajectory Conservation"]
        ]
        rows.append(row)

# === SAVE OUTPUT ===
df = pd.DataFrame(rows, columns=columns)
output_file = os.path.join(root_dir, f"{out_filename}.xlsx")
df.to_excel(output_file, index=False)
print(f"Saved: {output_file}")

