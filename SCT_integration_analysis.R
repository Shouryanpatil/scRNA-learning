# Single-cell RNA-seq - normalization

# Load libraries
library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)

# Normalize the counts
seurat_phase <- NormalizeData(filtered_seurat)

# Load cell cycle markers

load("C:/Users/SHOURYAN/Documents/SingleCellRNA/single_cell_rnaseq/single_cell_rnaseq/data/cycle.rda")

# Score cells for cell cycle
seurat_phase <- CellCycleScoring(seurat_phase,
                                 g2m.features = g2m_genes,
                                 s.features = s_genes)

# View cell cycle scores and phases assigned to cells
View(seurat_phase@meta.data)

# Identify the most variable genes
seurat_phase <- FindVariableFeatures(seurat_phase,
                                     selection.method = "vst",
                                     nfeatures = 2000,
                                     verbose = FALSE)

# Scale the counts
seurat_phase <- ScaleData(seurat_phase)

# Identify the 15 most hightly variable genes
ranked_variable_genes <- VariableFeatures(seurat_phase)
top_genes <- ranked_variable_genes[1:15]

# Plot the average expression and variance of these genes
# With labels to indicate which genes are in the top 15
p <- VariableFeaturePlot(seurat_phase)
print(p)
LabelPoints(plot = p, points = top_genes, repel = TRUE)
LabelPoints(plot = p, points = top_genes, repel = TRUE, xnudge = 0, ynudge = 0)

# Perform PCA
seurat_phase <- RunPCA(seurat_phase)

# Plot the PCA colored by cell cycle phase
DimPlot(seurat_phase,
        reduction = "pca",
        group.by = "Phase",
        split.by = "Phase")

# Check quartile values
summary(seurat_phase@meta.data$mitoRatio)

# Turn mitoRatio into categorical factor vector based on quartile values
seurat_phase@meta.data$mitoFr <- cut(seurat_phase@meta.data$mitoRatio,
                                     breaks = c(-Inf, 0.0144, 0.0199, 0.0267, Inf),
                                     labels = c("Low","Medium","Medium high","High"))

# Plot PCA colored by mitoFr
DimPlot(seurat_phase, reduction = "pca", group.by = "mitoFr", pt.size = 1) +
  ggtitle("PCA Colored by Mitochondrial Expression Quartiles")

# Split seurat object by condition to perform cell cycle scoring and SCT on all samples
split_seurat <- SplitObject(seurat_phase, split.by = "sample")

ctrl_reps <- split_seurat[c("ctrl_1", "ctrl_2")]

#Set Memory Limit(for Large Datasets)
options(future.globals.maxSize = 4000* 1024^2)

# Run SCTransform on Each Sample(with Regression)
for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- SCTransform(
    split_seurat[[i]],
    vars.to.regress = c("mitoRatio"),
    vst.flavor = "v2"
  )
}
# Inspecting Assays
split_seurat$ctrl@assays

# Check assays in 'stim' sample (e.g., if named "stim")
split_seurat$stim@assays

for (i in names(split_seurat)) {
  cat("Sample:", i, "\n")
  print(Assays(split_seurat[[i]]))
  cat("\n")
}

# First 10 features in the SCT assay
head(rownames(split_seurat$ctrl[["SCT"]]), 10)

# Top 10 variable features
head(VariableFeatures(split_seurat$ctrl), 10)

# Save the split seurat object
saveRDS(split_seurat, file = "C:/Users/SHOURYAN/Documents/SingleCellRNA/single_cell_rnaseq/single_cell_rnaseq/data/split_seurat.rds")

#Load the split seurat object into the environment
split_seurat <- readRDS("C:/Users/SHOURYAN/Documents/SingleCellRNA/single_cell_rnaseq/single_cell_rnaseq/data/split_seurat.rds")




# Select the most variable features to use for integration
integ_features <- SelectIntegrationFeatures(object.list = split_seurat,
                                            nfeatures = 3000)

# Prepare the SCT list object for integration
split_seurat <- PrepSCTIntegration(object.list = split_seurat,
                                   anchor.features = integ_features)

# Find best buddies = - can take a while to run
integ_anchors <- FindIntegrationAnchors(object.list = split_seurat,
                                        normalization.method = "SCT",
                                        anchor.features = integ_features)

# Integrate across conditions
seurat_integrated <- IntegrateData(anchorset = integ_anchors,
                                   normalization.method = "SCT")

# Run PCA
seurat_integrated <- RunPCA(object = seurat_integrated)

# Plot PCA
PCAPlot(seurat_integrated,
        split.by = "sample")

# Set seed
set.seed(123456)

# Run UMAP
seurat_integrated <- RunUMAP(seurat_integrated, 
                             dims = 1:40,
                             reduction = "pca")

# Plot UMAP                             
DimPlot(seurat_integrated)

# Plot UMAP split by sample
DimPlot(seurat_integrated,
        split.by = "sample")

# Save integrated seurat object 
saveRDS(seurat_integrated, "C:/Users/SHOURYAN/Documents/SingleCellRNA/single_cell_rnaseq/single_cell_rnaseq/results/integrated_seurat.rds")

# To save:
save.image("my_environment.RData")
# To load:
load("my_environment.RData")
