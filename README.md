# Comparative-Analysis-of-VAE-scRNAseq-Integration-Methods-scVI-MrVI-and-LDVAE


This repository contains code, data references, and results for a benchmarking study comparing three variational autoencoder (VAE) models—**scVI**, **MrVI**, and **LDVAE**—for batch integration in single-cell RNA sequencing (scRNA-seq) data. The evaluation is conducted across three biologically and technically diverse datasets:

- **Human Immune Cell Dataset** from the [OpenProblems project](https://openproblems.bio/datasets/openproblems_v1/immune_cells)
- **Capillary Blood PBMC Dataset** from [Zenodo 8020792](https://doi.org/10.5281/zenodo.8020792)
- **Remission Biome Pilot Dataset** from [Zenodo 11100300](https://doi.org/10.5281/zenodo.11100300)

The study systematically explores model performance using multiple batch correction and biological conservation metrics, under varying hyperparameter configurations and gene selection strategies.

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

## Datasets

This study utilizes three publicly released single-cell RNA-seq datasets, each representing distinct biological and technical contexts:

### 1. Human Immune Cell Dataset
- **Source:** Open Problems in Single-Cell Analysis
- **Description:** Contains 33,506 cells and 12,303 genes across 5 studies and 10 batches. Samples are derived from peripheral blood and bone marrow, sequenced using 10X Genomics (v2/v3) and Smart-seq2 technologies.
- **DOI:** [10.1038/s41592-021-01336-8](https://doi.org/10.1038/s41592-021-01336-8)
- **Link:** [https://openproblems.bio/datasets/openproblems_v1/immune_cells](https://openproblems.bio/datasets/openproblems_v1/immune_cells)

### 2. Zenodo 8020792 - Capillary Blood PBMC Dataset
- **Source:** Zenodo
- **Description:** Comprises 76,535 cells and 36,601 genes across 14 batches. Samples are collected from capillary blood using 10X Genomics.
- **DOI:** [10.5281/zenodo.8020792](https://doi.org/10.5281/zenodo.8020792)
- **Access:** Request-based (requires Zenodo login)

### 3. Zenodo 11100300 - Remission Biome Pilot Dataset
- **Source:** Zenodo
- **Description:** Contains 55,260 cells and 36,601 genes from 4 batches. Samples are collected from patients undergoing antibiotic-induced remission, sequenced using 10X Genomics.
- **DOI:** [10.5281/zenodo.11100300](https://doi.org/10.5281/zenodo.11100300)
- **Access:** Request-based (requires Zenodo login)



## Project Structure

```bash
.
├── notebooks/                # Google Drive Link to Jupyter notebooks for training and evaluation (size>25MB)
├── scripts/                  # Python scripts (in progress)
├── results/                  # Evaluation results and plots
├── environment.yml           # Reproducible environment definition (conda)
└── README.md                 # This file
