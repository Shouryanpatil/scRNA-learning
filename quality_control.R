# Load libraries
library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(ggplot2)
#install.packages("Seurat")
#install.packages("remotes")
#remotes::install_github("satijalab/seurat", ref = "release/5.0.1")  # Adjust version if needed
#install.packages("Matrix")
#packageVersion("Matrix")
#remove.packages("Seurat")
#remove.packages("SeuratObject")

library(Seurat)

# How to read in 10X data for a singl sample (output is a sparse matrix)
ctrl_counts <- Read10X(data.dir = "C:/Users/SHOURYAN/Documents/SingleCellRNA/single_cell_rnaseq/single_cell_rnaseq/data/ctrl_raw_feature_bc_matrix")

# Turn count matrix into a Seurat object (Output is a Seurat object)
ctrl <- CreateSeuratObject(counts = ctrl_counts,
                           min.features = 100)

# Explore the metadata
head(ctrl@meta.data)

Sample_names <- c("ctrl", "stim")

# Empty list to populate seurat object for each sample
list_seurat <- list()


for (sample in Sample_names) {
  #Path to data directory
  data_dir <- paste0("C:/Users/SHOURYAN/Documents/SingleCellRNA/single_cell_rnaseq/single_cell_rnaseq/data/", sample, "_raw_feature_bc_matrix")
  
  if (dir.exists(data_dir)) {
    # Create a Seurat object for each sample
    seurat_data <- Read10X(data.dir = data_dir)
    seurat_obj <- CreateSeuratObject(counts = seurat_data,
                                     min.features = 100,
                                     project = sample)
    # Save seurat object to list
    list_seurat[[sample]] <- seurat_obj
  } else {
    message(paste("❌ Directory does not exist:", data_dir))
  }
}

# Create a merged Seurat object
merged_seurat <- merge(x = list_seurat[["ctrl"]],
                       y = list_seurat[["stim"]],
                       add.cell.ids = c("ctrl", "stim"))

# Concatenate the count matrices of both sample together
merged_seurat <- merge(list_seurat[["ctrl"]], y = list_seurat[["stim"]],
                       add.cell.ids = c("ctrl", "stim"),
                       project = "merged")
head(merged_seurat@meta.data)
tail(merged_seurat@meta.data)

# Explore merged metadata
View(merged_seurat@meta.data)

# Add number of genes per UMI for each cell to metadata
merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)

# Compute percent mito ratio
merged_seurat$mitoRatio <- PercentageFeatureSet(merged_seurat, pattern = "^MT-") / 100

# Create metadata dataframe
metadata <- merged_seurat@meta.data

# Add cell IDs to metadata
metadata$cells <- rownames(metadata)

library(stringr)
# Create sample column
metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "^ctrl_"))] <- "ctrl"
metadata$sample[which(str_detect(metadata$cells, "^stim_"))] <- "stim"

# Rename columns 
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident, 
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)

# Add metadata back to Seurat object
merged_seurat@meta.data <- metadata

# Create .RData object to Load at any time
save(merged_seurat, file = "C:/Users/SHOURYAN/Documents/SingleCellRNA/single_cell_rnaseq/single_cell_rnaseq/data/merged_filtered_seurat.RData")

# Visualaize the number of cell counts per sample
metadata %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")

# Visualize the number UMIs/transcripts per cell
metadata %>%
  ggplot(aes(color = sample, x=nUMI, fill = sample)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)

# Visualize the distribution of genes detected per cell via histogram
metadata %>%
  ggplot(aes(color=sample, x=nGene, fill = sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  scale_x_log10() +
  geom_vline(xintercept = 300)

# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI (novelty score)
metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill = sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)

# Visualize the the distribution of mitochondrial gene expression detected per cell
metadata %>%
  ggplot(aes(color=sample, x=mitoRatio, fill = sample)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  geom_vline(xintercept = 0.2)

# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
metadata %>%
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) +
  geom_point() +
  scale_color_gradient(low = "gray90", high = "black") +
  stat_smooth(method = lm) +
  scale_x_log10() +
  scale_y_log10() +
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample)

# Filter out low quality cells using selected thresholds - these will change with experiment
filtered_seurat <- subset(x = merged_seurat,
                          subset = (nUMI >= 500) &
                            (nGene >= 250) &
                            (log10GenesPerUMI > 0.80) &
                            (mitoRatio < 0.20))

# Extract counts
counts <- GetAssayData(object = filtered_seurat, slot = "counts")

# Output a logical matrix specifying for each gene on whether or not there are more than zero counts per cell
nonzero <- counts > 0

# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10

# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]

# Reassign to filtered Seurat object 
filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)

# Exercises
# Extract metdata from the filtered Seurat object
metadata_clean <- filtered_seurat@meta.data

table(metadata_clean$orig.ident)

#Save Filtered Object
save(filtered_seurat, file = "C:/Users/SHOURYAN/Documents/SingleCellRNA/single_cell_rnaseq/single_cell_rnaseq/data/seurat_filtered.RData")


# To save:
save.image("my_environment.RData")
# To load:
load("my_environment.RData")
