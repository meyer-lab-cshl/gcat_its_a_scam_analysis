library(Seurat)
library(Signac)
library(pheatmap)
library(tidyverse)

# set the seed
set.seed(42)

# Load the dataset 
pbmc <- LoadSeuratRds("pbmc3k_joint_matrix.Rds")

# normalize RNA assay
DefaultAssay(pbmc) <- "RNA"
pbmc <- SCTransform(pbmc)
pbmc <- RunPCA(pbmc)

# normalize ATAC assaygt
DefaultAssay(pbmc) <- "ATAC"
pbmc <- FindTopFeatures(pbmc, min.cutoff = 5)
pbmc <- RunTFIDF(pbmc)
pbmc <- RunSVD(pbmc)
pbmc


# cluster
pbmc <- FindMultiModalNeighbors(pbmc, reduction.list = list("pca", "lsi"),
                                dims.list = list(1:50, 2:40))
pbmc <- RunUMAP(pbmc, nn.name = "weighted.nn", reduction.name = "wnn.umap")
pbmc <- FindClusters(pbmc, graph.name = "wsnn", algorithm = 3, resolution = 0.4)

# plot the clustering
pdf("pbmc3k_joint_cluster.pdf", height = 3, width = 3)
DimPlot(pbmc, reduction = "wnn.umap", label = TRUE) + NoLegend()
dev.off ()


# plot the marker genes
DefaultAssay(pbmc) <- "SCT"

azimuth_l1 <- list(
  Mono    = c("CTSS","FCN1","NEAT1","LYZ","PSAP","S100A9","AIF1","MNDA","SERPINA1","TYROBP"),
  CD4_T   = c("IL7R","MAL","LTB","CD4","LDHB","TPT1","TRAC","TMSB10","CD3D","CD3G"),
  CD8_T   = c("CD8B","CD8A","CD3D","TMSB10","HCST","CD3G","LINC02446","CTSW","CD3E","TRAC"),
  NK      = c("NKG7","KLRD1","TYROBP","GNLY","FCER1G","PRF1","CD247","KLRF1","CST7","GZMB"),
  B       = c("CD79A","RALGPS2","CD79B","MS4A1","BANK1","CD74","TNFRSF13C","HLA-DQA1","IGHM","MEF2C"),
  other_T = c("CD3D","TRDC","GZMK","KLRB1","NKG7","TRGC2","CST7","LYAR","KLRG1","GZMA"),
  DC      = c("CD74","HLA-DPA1","HLA-DPB1","HLA-DQA1","CCDC88A","HLA-DRA","HLA-DMA","CST3","HLA-DQB1","HLA-DRB1")
)


# calculate madule scores
pbmc <- AddModuleScore(pbmc, features = azimuth_l1, name = "l1_", seed = 1, nbin = 15)
head(pbmc)

# rename the module score metadata with the cell type names
score_cols <- paste0("l1_", seq_along(azimuth_l1))
new_cols   <- paste0(names(azimuth_l1), "_score")
names(pbmc@meta.data)[match(score_cols, names(pbmc@meta.data))] <- new_cols

# calculate the median module scores
scores_by_cluster <- aggregate(pbmc[[new_cols]], 
  by = list(cluster = Idents(pbmc)), FUN = median)
scores_by_cluster[new_cols] <- round(scores_by_cluster[new_cols], 2)
scores_by_cluster

# heatmap of median module scores
mat <- as.matrix(scores_by_cluster[, new_cols])
rownames(mat) <- scores_by_cluster$cluster
colnames(mat) <- sub("_score$", "", colnames(mat))

lim <- max(abs(mat))
breaks <- seq(-lim, lim, length.out = 101)

pdf("pbmc3k_joint_cluster_module_score.pdf", height = 3, width = 3)
pheatmap(mat,
         color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
         breaks = breaks,
         cluster_rows = FALSE, cluster_cols = FALSE,
         display_numbers = TRUE, number_format = "%.2f")
dev.off()


# manually label the cell types
l1_labels <- c("0" = "Mono",  "1" = "CD4_T", "2" = "CD4_T", "3" = "CD8_T",
               "4" = "B",     "5" = "CD8_T", "6" = "Mono",  "7" = "NK",
               "8" = "DC",    "9" = "DC")
pbmc$celltype.l1 <- unname(l1_labels[as.character(pbmc$wsnn_res.0.4)])
table(pbmc$wsnn_res.0.4, pbmc$celltype.l1)


# write the barcode celltypes
bc_celltype <- tibble(bc = colnames(pbmc), cell_type = pbmc$celltype.l1)
write_tsv(bc_celltype, "pbmc3k_joint_celltype.tsv", col_names = FALSE)

# plot the PBMC marker genes for reference 
l1_markers <- c("CD3E","CD4","IL7R","CD8A","CD8B","TRDC","KLRB1",
                "GNLY","NKG7","KLRD1",
                "MS4A1","CD79A","MZB1",
                "CD14","LYZ","FCGR3A",
                "FCER1A","CD1C","LILRA4","MKI67")
pdf("pbmc3k_joint_cluster_pbmc_markers.pdf", height = 4, width = 6)   
DotPlot(pbmc, features = l1_markers) + RotatedAxis()
dev.off()
