# Comparative-Analysis-of-VAE-scRNAseq-Integration-Methods-scVI-MrVI-and-LDVAE


This repository contains code, data references, and results for a benchmarking study comparing three variational autoencoder (VAE) models—**scVI**, **MrVI**, and **LDVAE**—on the task of batch integration in single-cell RNA sequencing (scRNA-seq) data. The evaluation is conducted using the **Human Immune Cell Dataset** from the [OpenProblems project](https://openproblems.bio/datasets/openproblems_v1/immune_cells).

## Project Overview

The goal is to assess and compare integration performance using:
- **4 batch correction metrics**
- **7 biological conservation metrics**
- **UMAP and t-SNE visualizations**
- **Scalability (training time, GPU usage)**
- **Ease of use / usability**


## Main Dependencies
Package	Version	Description
Python ==	3.12.2	Programming language

scvi-tools	1.3.0	Variational inference for scRNA-seq

Scanpy ==	1.11.0	Single-cell data analysis

PyTorch ==	2.6.0+cu124	Deep learning framework

JAX == 0.4.35	Accelerated numerical computing (used for MrVI)

scIB ==	1.1.7	Integration benchmarking metrics

matplotlib ==	3.x	Plotting and visualization

Full list of packages: see environment.yml

## Project Structure

```bash
.
├── notebooks/                # Google Drive Link to Jupyter notebooks for training and evaluation (size>25MB)
├── scripts/                  # Python scripts (in progress)
├── results/                  # Evaluation results and plots
├── environment.yml           # Reproducible environment definition (conda)
└── README.md                 # This file
