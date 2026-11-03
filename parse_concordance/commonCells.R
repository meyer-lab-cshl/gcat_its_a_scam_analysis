library(Seurat)
library(tidyverse)
library(UpSetR)

h1.parse <- LoadSeuratRds("../../../count_matrix/og/h1_og.Rds")
h2.parse <- LoadSeuratRds("../../../count_matrix/og/h2_og.Rds")
h1.star <- LoadSeuratRds("../../../count_matrix/h1/h1_n500_collapse.Rds")
h2.star <- LoadSeuratRds("../../../count_matrix/h2/h2_n500_collapse.Rds")

h1.cells.ids <- list(parse = colnames(h1.parse), star_scam = colnames(h1.star))

pdf("h1_common_cells.pdf", height=4, width = 4, onefile=F)
upset(fromList(h1.cells.ids), order.by="freq", mb.ratio = c(0.55, 0.45)) 
dev.off()

h2.cells.ids <- list(parse = colnames(h2.parse), star_scam = colnames(h2.star))

pdf("h2_common_cells.pdf", height=4, width = 4, onefile=F)
upset(fromList(h2.cells.ids), order.by="freq", mb.ratio = c(0.55, 0.45)) 
dev.off()
