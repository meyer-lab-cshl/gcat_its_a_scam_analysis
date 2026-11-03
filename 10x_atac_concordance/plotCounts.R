library(Seurat)
library(Signac)
library(ggpubr)
library(cowplot)
library(tidyverse)

bwa <- LoadSeuratRds("bwa_matrix.Rds")
bwa

cr <- LoadSeuratRds("cr_matrix.Rds")
colnames(cr) <- substr(colnames(cr), 1, 16)
cr


bwa.c <- tibble(cells = colnames(bwa), 
  bwa.counts = colSums(GetAssayData(bwa)), 
  bwa.features = colSums(GetAssayData(bwa) > 0))
cr.c <- tibble(cells = colnames(cr), 
  cr.counts = colSums(GetAssayData(cr)), 
  cr.features = colSums(GetAssayData(cr) > 0))

bwa.cr <- bwa.c %>% left_join(cr.c)
bwa.cr

ggplot(data = bwa.cr, aes(x = bwa.counts, y = cr.counts)) + 
  geom_point(size = 0.5) + 
  stat_cor(method = "pearson", label.x = 0, label.y = 75000) +
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "BWA", y = "Cellranger") +
  theme_cowplot(10) 
ggsave("pbmc_granulocyte_sorted_3k_bwa_cr_nCounts.pdf", height = 4, width = 4)

ggplot(data = bwa.cr, aes(x = bwa.features, y = cr.features)) + 
  geom_point(size = 0.5) + 
  stat_cor(method = "pearson", label.x = 0, label.y = 30000) +
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "BWA", y = "Cellranger") +
  theme_cowplot(10) 
ggsave("pbmc_granulocyte_sorted_3k_bwa_cr_nFeatures.pdf", height = 4, width = 4)
