library(Seurat)
library(ggpubr)
library(cowplot)
library(tidyverse)

h1.parse <- LoadSeuratRds("../../../count_matrix/og/h1_og.Rds")
h1.star <- LoadSeuratRds("../../../count_matrix/h1/h1_n500_collapse.Rds")

h1.common <- intersect(colnames(h1.parse), colnames(h1.star))
length(h1.common)

h1.common.genes <- intersect(rownames(h1.parse), rownames(h1.star))
length(h1.common.genes)

h1.parse
h1.parse <- h1.parse[, colnames(h1.parse) %in% h1.common]
h1.star
h1.star <- h1.star[, colnames(h1.star) %in% h1.common]

h1.star.mat <- as.matrix(GetAssayData(h1.star))
h1.parse.mat <- as.matrix(GetAssayData(h1.parse))

h1.star.mat <- h1.star.mat[h1.common.genes,]
h1.parse.mat <- h1.parse.mat[h1.common.genes,]
dim(h1.star.mat)
dim(h1.parse.mat)
sum(rownames(h1.star.mat) == rownames(h1.star.mat))



h1.cor <- c()
for (i in h1.common) {
  h1.cor <- c(h1.cor,
    cor(h1.star.mat[, i], h1.parse.mat[, i], method="spearman"))
}
length(h1.cor)

h1.cor <- tibble(corr = h1.cor)
h1.cor
median(h1.cor$corr)

ggplot(data = h1.cor) +
  geom_histogram(mapping = aes(corr), binwidth = 0.001) +
  labs(x = "Spearman's correlation", y = "Number of cells") +
  theme_cowplot(10)
ggsave("h1_gene_spearman.pdf", height = 5, width = 5)



h2.parse <- LoadSeuratRds("../../../count_matrix/og/h2_og.Rds")
h2.star <- LoadSeuratRds("../../../count_matrix/h2/h2_n500_collapse.Rds")

h2.common <- intersect(colnames(h2.parse), colnames(h2.star))
length(h2.common)

h2.common.genes <- intersect(rownames(h2.parse), rownames(h2.star))
length(h2.common.genes)


h2.parse
h2.parse <- h2.parse[, colnames(h2.parse) %in% h2.common]
h2.star
h2.star <- h2.star[, colnames(h2.star) %in% h2.common]

h2.star.mat <- as.matrix(GetAssayData(h2.star))
h2.parse.mat <- as.matrix(GetAssayData(h2.parse))

h2.star.mat <- h2.star.mat[h2.common.genes,]
h2.parse.mat <- h2.parse.mat[h2.common.genes,]
dim(h2.star.mat)
dim(h2.parse.mat)
sum(rownames(h2.star.mat) == rownames(h2.star.mat))



h2.cor <- c()
for (i in h2.common) {
  h2.cor <- c(h2.cor,
    cor(h2.star.mat[, i], h2.parse.mat[, i], method="spearman"))
}
length(h2.cor)

h2.cor <- tibble(corr = h2.cor)
h2.cor
median(h2.cor$corr)

ggplot(data = h2.cor) +
  geom_histogram(mapping = aes(corr), binwidth = 0.001) +
  labs(x = "Spearman's correlation", y = "Number of cells") +
  theme_cowplot(10)
ggsave("h2_gene_spearman.pdf", height = 5, width = 5)
