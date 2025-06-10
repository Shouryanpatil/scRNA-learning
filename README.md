# scRNA-learning

## Single-cell RNA-seq (scRNA-seq) Analysis Pipeline using Seurat

Welcome to the scRNA-seq data analysis pipeline using the **Seurat** package in R. This repository contains three key R scripts that guide you through the end-to-end analysis of scRNA-seq data: from quality control to data integration and clustering.

## Repository Structure

```
├── quality_control.R                 # Preprocess, QC and filter scRNA-seq data
├── SCT_integration_analysis.R       # Normalize and integrate multiple samples using SCTransform
├── clustering.R                     # Perform clustering and cell type annotation
└── README.md                        # This file
```

---

## 1️] quality\_control.R

**Purpose:** Performs quality control on raw scRNA-seq data (10X Genomics format) and filters low-quality cells/genes.

### 📥 Input:

* Raw 10X data from two samples: `ctrl_raw_feature_bc_matrix/`, `stim_raw_feature_bc_matrix/`

### 📤 Output:

* `merged_filtered_seurat.RData`: Merged object before gene filtering
* `seurat_filtered.RData`: Final Seurat object with high-quality cells and genes

### 🔁 Key Steps:

* Read and merge raw data
* Calculate quality metrics: `log10GenesPerUMI`, `mitoRatio`
* Visualize metrics for QC decisions
* Filter cells by:

  * nUMI >= 500
  * nGene >= 250
  * log10GenesPerUMI > 0.8
  * mitoRatio < 0.2
* Remove lowly expressed genes (keep genes expressed in >= 10 cells)
* Save filtered Seurat object

---

## 2️] SCT\_integration\_analysis.R

**Purpose:** Normalizes, scores cell cycles, and integrates multiple samples using **SCTransform (SCT)** to minimize batch effects.

### 📥 Input:

* Filtered Seurat object from previous step: `seurat_filtered.RData`
* Cell cycle genes: `cycle.rda`

### 📤 Output:

* `integrated_seurat.rds`: Final integrated Seurat object
* `split_seurat.rds`: Intermediate objects
* PCA and UMAP visualizations

### 🔁 Key Steps:

* Normalize data and score for cell cycle (G2M, S)
* Run PCA and visualize
* Split samples and run SCTransform
* Identify integration features and anchors
* Run `IntegrateData` to remove batch effects
* Visualize integrated data (UMAP by sample)
* Save final integrated object

---

## 3️] clustering.R

**Purpose:** Performs clustering to identify distinct cell types/states, visualize data, and annotate clusters.

### 📥 Input:

* Integrated Seurat object: `integrated_seurat.rds`
* Annotation file (CSV with gene descriptions)

### 📤 Output:

* `seurat_labelled.rds`: Annotated Seurat object with cell types
* `sessionInfo_scrnaseq_Feb2023.txt`: Session info
* UMAP plots, bar plots, marker gene tables

### 🔁 Key Steps:

1. Run PCA and Elbow Plot to select dimensions
2. Identify nearest neighbors and clusters
3. Run UMAP for 2D visualization
4. Explore metadata by sample, cluster, cell cycle
5. Generate FeaturePlot, DotPlot
6. Identify marker genes for clusters
7. Annotate cell types (e.g., monocytes, B cells)
8. Remove unwanted/stressed clusters
9. Save outputs

---

## ✅ Overall Workflow Summary

```
Raw 10X Data ─▶ quality_control.R ─▶ Filtered Seurat Object
                     │
                     ▼
          SCT_integration_analysis.R ─▶ Integrated Seurat Object
                     │
                     ▼
                 clustering.R ─▶ Final Annotated Seurat Object & Results
```

##  Notes

* All objects are saved for reuse in `.RData` or `.rds` formats
* Visualizations and intermediate objects help with reproducibility and QC validation
* This pipeline is modular — each script can be rerun independently if input files are present

---
