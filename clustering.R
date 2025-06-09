# Single-cell RNA-seq - clustering

# Load libraries
library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)

# Explore heatmap of PCs
DimHeatmap(seurat_integrated,
           dims = 1:9,
           cells = 500,
           balanced = TRUE)

# Printing out the most variable genes driving PCs
print(x = seurat_integrated[["pca"]],
      dims = 1:10,
      nfeatures = 5)

# Plot the elbow plot
ElbowPlot(object = seurat_integrated,
          ndims = 40)

# Determine the K-nearest neighbor graph
seurat_integrated <- FindNeighbors(object = seurat_integrated,
                                   dims = 1:40)

# Determine the clusters for various resolutions
seurat_integrated <- FindClusters(object = seurat_integrated,
                                  resolution = c(0.4, 0.6, 0.8, 1.0, 1.4))

# Explore resolution
seurat_integrated@meta.data %>%
  View()

# Assign identity of clusters
Idents(object = seurat_integrated) <- "integrated_snn_res.0.8"

# Plot the UMAP
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 6)

# DON'T run this line:
# load(bzfile("data/additional_data/seurat_integrated.RData.bz2"))

# Assign identity of clusters
Idents(object = seurat_integrated) <- "integrated_snn_res.0.8"

# Plot the UMAP
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 6)

# Extract identity and sample information from seurat object to determine the number of cells per cluster per sample
n_cells <- FetchData(seurat_integrated,
                     vars = c("ident", "sample")) %>%
  dplyr::count(ident, sample)

# Barplot of number of cells per cluster by sample
ggplot(n_cells, aes(x=ident, y=n, fill = sample)) +
  geom_bar(position = position_dodge(), stat = "identity") +
    geom_text(aes(label = n), vjust = -.2, position = position_dodge(1))

# UMAP of cells in each cluster by sample
DimPlot(seurat_integrated,
        label = TRUE,
        split.by = "sample") + NoLegend()

# Barplot of proportion of cells in each cluster by sample
ggplot(seurat_integrated@meta.data) +
  geom_bar(aes(x=integrated_snn_res.0.8, fill=sample), position=position_fill())

# Explore wheter clusters segeregate by cell cycle phase 
DimPlot(seurat_integrated,
        label = TRUE,
        split.by = "Phase") +NoLegend()

# Determine metrics to plot presetn in seurat_integrate@meta.dat
metrics <- c("nUMI", "nGene", "S.Score", "G2M.Score", "mitoRatio")

FeaturePlot(seurat_integrated,
            reduction = "umap",
            features = metrics,
            pt.size = 0.4,
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE)
# Boxplot of nGene per cluster
ggplot(seurat_integrated@meta.data) +
  geom_boxplot(aes(x=integrated_snn_res.0.8, y=nGene, fill = integrated_snn_res.0.8)) +
  NoLegend()

# Definingg the information in the seurat object of interest 
columns <- c(paste0("PC_", 1:16),
             "ident",
             "UMAP_1", "UMAP_2")

# Extracting this data from the seurat object
pc_data <- FetchData(seurat_integrated,
                     vars = columns)
# Extract the UMAP coordinates for the first 10 cells
seurat_integrated@reductions$umap@cell.embeddings[1:10, 1:2]

# Defining the information in the seurat object of interest 
columns <- c(paste0("PC_", 1:16),
             "ident",
             "umap_1", "umap_2")

# Adding cluster label to center of cluster on UMAP
umap_label <- FetchData(seurat_integrated, 
                        vars = c("ident", "UMAP_1", "UMAP_2"))  %>%
  group_by(ident) %>%
  dplyr::summarise(x=mean(UMAP_1), y=mean(UMAP_2))

# Plotting a UMAP plot for each of the PCs
map(paste0("PC_", 1:16), function(pc){
  ggplot(pc_data, 
         aes(UMAP_1, UMAP_2)) +
    geom_point(aes_string(color=pc), 
               alpha = 0.7) +
    scale_color_gradient(guide = FALSE, 
                         low = "grey90", 
                         high = "blue")  +
    geom_text(data=umap_label, 
              aes(label=ident, x, y)) +
    ggtitle(pc)
}) %>% 
  plot_grid(plotlist = .)

# Examine PCA results
print(seurat_integrated[["pca"]], dims = 1:5, nfeatures = 5)

DimPlot(object = seurat_integrated,
        reduction = "umap",
        label = TRUE) + NoLegend()

# Select the RNA counts slot to be the default assay
DefaultAssay(seurat_integrated) <- "RNA"

# Normalize RNA data for visualization purposes
seurat_integrated <- NormalizeData(seurat_integrated, verbose = FALSE)

# CD14+ monocyte markers
FeaturePlot(seurat_integrated,
            reduction = "umap",
            features = c("CD14", "LYZ"),
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE)

# FCGR3A+ monocyte markers
FeaturePlot(seurat_integrated,
            reduction = "umap",
            features = c("FCGR3A", "MS4A7"),
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE)

# Macrophages
FeaturePlot(seurat_integrated,
            reduction = "umap",
            features = c("MARCO", "ITGAM", "ADGRE1"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)

# Conventional dendritic cell markers
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("FCER1A", "CST3"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)

# Plasmacytoid dendritic cell markers
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("IL3RA", "GZMB", "SERPINF1", "ITM2C"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)

# List of known celltype markers
markers <- list()
markers[["CD14+ monocytes"]] <- c("CD14", "LYZ")
markers[["FCGR3A+ monocyte"]] <- c("FCGR3A", "MS4A7")
markers[["Macrophages"]] <- c("MARCO", "ITGAM", "ADGRE1")
markers[["Conventional dendritic"]] <- c("FCER1A", "CST3")
markers[["Plasmacytoid dendritic"]] <- c("IL3RA", "GZMB", "SERPINF1", "ITM2C")

# Create dotplot based on RNA expression
DotPlot(seurat_integrated, markers, assay="RNA")

# To save:
save.image("my_environment2.RData")
# To load:
load("my_environment2.RData")



DefaultAssay(seurat_integrated) <- "RNA"

cluster0_conserved_markers <- FindConservedMarkers(seurat_integrated,
                                                   ident.1 = 0,
                                                   grouping.var = "sample",
                                                   only.pos = TRUE,
                                                   logfc.threshold = 0.25)

view(cluster0_conserved_markers)

annotations <- read.csv("C:/Users/SHOURYAN/Documents/SingleCellRNA/single_cell_rnaseq/single_cell_rnaseq/data/annotation.csv")

# Combine markers with gene descriptions
cluster0_ann_markers <- cluster0_conserved_markers %>%
  rownames_to_column(var = "gene") %>% 
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))

View(cluster0_ann_markers)

# Create function to get conserved markers for any given cluster
get_conserved <- function(cluster){
  FindConservedMarkers(seurat_integrated,
                       ident.1 = cluster,
                       grouping.var = "sample",
                       only.pos = TRUE) %>%
    rownames_to_column(var = "gene") %>%
    left_join(y = unique(annotations[, c("gene_name", "description")]),
              by = c("gene" = "gene_name")) %>%
    cbind(cluster_id = cluster, .)
}

# Iterate function across desired clusters
conserved_markers <- map_dfr(c(4,0,6,2), get_conserved)

# Extract top 10 markers per cluster
top10 <- conserved_markers %>% 
  mutate(avg_fc = (ctrl_avg_log2FC + stim_avg_log2FC) / 2) %>% 
  group_by(cluster_id) %>% 
  top_n(n = 10, wt = avg_fc)

# Visualize top 10 markers per cluster
View(top10)

# Plot interesting marker gene expression for cluster 4
FeaturePlot(object = seurat_integrated,
            features = c("HSPH1", "HSPE1", "DNAJB1"),
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE,
            repel = TRUE)

# Vln plot - cluster 4
VlnPlot(object = seurat_integrated,
        features = c("HSPH1", "HSPE1", "DNAJB1"))

# Determine differentiating markers for CD4+ T cell
cd4_tcells <- FindMarkers(seurat_integrated,
                          ident.1 = 2,
                          ident.2 = c(0,4,6))

# Add gene symbols to the DE table
cd4_tcells <- cd4_tcells %>%
  rownames_to_column(var = "gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))

# Reorder columns and sort by padj
cd4_tcells <- cd4_tcells[, c(1, 3:5, 2, 6:7)]

cd4_tcells <- cd4_tcells %>%
  dplyr::arrange(p_val_adj)

# View data
View(cd4_tcells)

# Rename all indentities
seurat_integrated <- RenameIdents(object = seurat_integrated,
                                  "0" = "Naive or memory CD4+ T cells",
                                  "1" = "CD14+ monocytes",
                                  "2" = "Activated T cells",
                                  "3" = "CD14+ monocytes",
                                  "4" = "Stressed cells / Unknown",
                                  "5" = "CD8+ T cells",
                                  "6" = "Naive or memory CD4+ T cells",
                                  "7" = "B cells",
                                  "8" = "NK cells",
                                  "9" = "CD8+ T cells",
                                  "10" = "FCGR3A+ monocytes",
                                  "11" = "B cells",
                                  "12" = "NK cells",
                                  "13" = "B cells",
                                  "14" = "Conventional dendritic cells",
                                  "15" = "Megakaryocytes",
                                  "16" = "Plasmacytoid dendritic cells")

# Plot the UMAP
DimPlot(object = seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 3,
        repel = TRUE)

# Remove the stressed or dying cells
seurat_subset_labeled <- subset(seurat_integrated,
                                idents = "Stressed cells / Unknown", invert = TRUE)

# Re-visualize the clusters
DimPlot(object = seurat_subset_labeled, 
        reduction = "umap", 
        label = TRUE,
        label.size = 3,
        repel = TRUE)

# Save final R object 
write_rds(seurat_integrated,
          file = "C:/Users/SHOURYAN/Documents/SingleCellRNA/single_cell_rnaseq/single_cell_rnaseq/results/seurat_labelled.rds")

# Create and save a text file with sessionInfo
sink("sessionInfo_scrnaseq_Feb2023.txt")
sessionInfo()
sink()

# Add celltype annotation as a column in meta.data
seurat_subset_labeled$celltype <- Idents(seurat_subset_labeled)

# Compute number of cells per celltype
n_cells <- FetchData(seurat_subset_labeled,
                     vars = c("celltype", "sample")) %>%
  dplyr::count(celltype, sample)

# Barplot of number of cells per celltype by sample
ggplot(n_cells, aes(x=celltype, y=n, fill = sample)) +
  geom_bar(position = position_dodge(), stat = "identity") +
  geom_text(aes(label = n), vjust = -.2, position = position_dodge(1))

# Subset seurat object to just B cells
seurat_b_cells <- subset(seurat_subset_labeled, subset = (celltype == "B cells"))

# Run a wilcox test to compare ctrl vs stim
b_markers <- FindMarkers(seurat_b_cells,
                         ident.1 = "ctrl",
                         ident.2 = "stim",
                         grouping.var = "sample",
                         only.pos = FALSE,
                         logfc.threshold = 0.25)

# Idents(seurat_b_cells)
# unique(Idents(seurat_b_cells))
# unique(Idents(seurat_b_cells))
# 
# Idents(seurat_b_cells) <- "condition"
# b_markers <- FindMarkers(seurat_b_cells,
#                          ident.1 = "ctrl",
#                          ident.2 = "stim",
#                          only.pos = FALSE,
#                          logfc.threshold = 0.25)

# unique(seurat_b_cells@meta.data$condition)
# table(seurat_b_cells@meta.data$condition)
# 
# Idents(seurat_b_cells) <- "condition"
# colnames(seurat_b_cells@meta.data)
# unique(seurat_b_cells@meta.data$orig.ident)
# unique(seurat_b_cells@meta.data$sample)
# 
# # Example: if 'orig.ident' is like 'ctrl_rep1', 'stim_rep2'
# seurat_b_cells$condition <- ifelse(grepl("ctrl", seurat_b_cells$orig.ident), "ctrl", "stim")
# 
# unique(seurat_b_cells@meta.data$sample)
# # If it shows "ctrl", "stim", "ctrl", "stim", you can just do:
# seurat_b_cells$condition <- seurat_b_cells$sample
# 
# Idents(seurat_b_cells) <- "condition"
# 
# b_markers <- FindMarkers(seurat_b_cells,
#                          ident.1 = "ctrl",
#                          ident.2 = "stim",
#                          only.pos = FALSE,
#                          logfc.threshold = 0.25)

BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
EnhancedVolcano(b_markers,
                row.names(b_markers),
                x="avg_log2FC",
                y="p_val_adj"
)

# Idents(seurat_b_cells)
# DimPlot(seurat_b_cells, reduction = "umap", label = TRUE)
# FeaturePlot(seurat_b_cells, features = c("CD19", "MS4A1", "CD79A"))
# markers <- FindMarkers(seurat_object, ident.1 = "B cells")
# head(markers)
# b_cells_only <- subset(seurat_object, idents = "B cells")
# table(Idents(seurat_object))
# 
# table(Idents(seurat_subset_labeled))
# table(Idents(seurat_data))
# table(Idents(seurat_integrated))
# table(Idents(seurat_obj))
# table(Idents(seurat_phase))
# table(Idents(seurat_b_cells))
# 
# seurat_object <- seurat_subset_labeled
# markers <- FindMarkers(seurat_object, ident.1 = "B cells")
# head(markers)
# b_cells_only <- subset(seurat_object, idents = "B cells")
# 
# 
# head(seurat_b_cells@meta.data)
# 
# b_markers <- FindMarkers(seurat_b_cells,
#                          ident.1 = "ctrl",
#                          ident.2 = "stim",
#                          group.by = "sample",
#                          only.pos = FALSE,
#                          logfc.threshold = 0.25)


# Bring in Seurat object
seurat <- readRDS("C:/Users/SHOURYAN/Documents/SingleCellRNA/single_cell_rnaseq/single_cell_rnaseq/results/integrated_seurat.rds")

# Extract raw counts and metadata to create SingleCellExperiment object
counts <- seurat@assays$RNA@counts 

metadata <- seurat@meta.data

# Set up metadata as desired for aggregation and DE analysis
metadata$cluster_id <- factor(seurat@active.ident)

# Create single cell experiment object
sce <- SingleCellExperiment(assays = list(counts = counts), 
                            colData = metadata)

# To save:
save.image("my_environment4.RData")
# To load:
load("my_environment4.RData")