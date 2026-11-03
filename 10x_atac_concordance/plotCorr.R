library(Seurat)
library(Signac)
library(cowplot)
library(tidyverse)

bwa <- LoadSeuratRds("bwa_matrix.Rds")
bwa

cr <- LoadSeuratRds("cr_matrix.Rds")
colnames(cr) <- substr(colnames(cr), 1, 16)
cr

common.peaks <- intersect(rownames(bwa), rownames(cr))
bwa <- bwa[rownames(bwa) %in% common.peaks,] 
bwa
cr <- cr[rownames(cr) %in% common.peaks,] 
cr

common.cells <- intersect(colnames(bwa), colnames(cr))

bwa.mat <- as.matrix(GetAssayData(bwa))
cr.mat <- as.matrix(GetAssayData(cr))

bwa.cr.cor <- c()
for (i in common.cells) {
  bwa.cr.cor <- c(bwa.cr.cor,
    cor(bwa.mat[, i], cr.mat[, i], method="spearman"))
}
length(bwa.cr.cor)

bwa.cr.cor <- tibble(corr = bwa.cr.cor)
bwa.cr.cor
median(bwa.cr.cor$corr)

ggplot(data = bwa.cr.cor) +
  geom_histogram(mapping = aes(corr), binwidth = 0.001) +
  labs(x = "Spearman's correlation", y = "Number of cells") +
  theme_cowplot(10)
ggsave("pbmc_granulocyte_sorted_3k_bwa_cr_spearman.pdf", height = 5, width = 5)
