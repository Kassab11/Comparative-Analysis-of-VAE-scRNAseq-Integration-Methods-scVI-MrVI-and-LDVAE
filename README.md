# Comparative-Analysis-of-VAE-scRNAseq-Integration-Methods-scVI-MrVI-and-LDVAE


This repository contains code, data references, and results for a benchmarking study comparing three variational autoencoder (VAE) models—**scVI**, **MrVI**, and **LDVAE**—on the task of batch integration in single-cell RNA sequencing (scRNA-seq) data. The evaluation is conducted using the **Human Immune Cell Dataset** from the [OpenProblems project](https://openproblems.bio/datasets/openproblems_v1/immune_cells).

## Project Overview

The goal is to assess and compare integration performance using:
- **4 batch correction metrics**
- **7 biological conservation metrics**
- **UMAP and PCA visualizations**
- **Scalability (training time, GPU usage)**
- **Ease of use / usability**

## Project Structure

```bash
.
├── notebooks/                # Google Drive Link to Jupyter notebooks for training and evaluation (size>25MB)
├── scripts/                  # Python scripts (in progress)
├── results/                  # Evaluation results and plots
├── environment.yml           # Reproducible environment definition (conda)
└── README.md                 # This file
