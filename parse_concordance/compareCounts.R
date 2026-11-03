library(Seurat)
library(ggpubr)
library(cowplot)
library(tidyverse)

h1.parse <- LoadSeuratRds("../../../count_matrix/og/h1_og.Rds")
h1.star <- LoadSeuratRds("../../../count_matrix/h1/h1_n500_collapse.Rds")

h1.common <- intersect(colnames(h1.parse), colnames(h1.star))
length(h1.common)

h1.parse
h1.parse <- h1.parse[, colnames(h1.parse) %in% h1.common]
h1.star
h1.star <- h1.star[, colnames(h1.star) %in% h1.common]

h1.star.counts <- tibble(cell_id = colnames(h1.star), 
  star_scam_counts = h1.star$nCount_RNA, 
  star_scam_features = h1.star$nFeature_RNA)
h1.star.counts
h1.parse.counts <- tibble(cell_id = colnames(h1.parse), 
  parse_counts = h1.parse$nCount_RNA, 
  parse_features = h1.parse$nFeature_RNA)
h1.parse.counts

h1.counts <- h1.star.counts %>% 
  left_join(h1.parse.counts)

ggplot(data = h1.counts, aes(x = star_scam_counts, y = parse_counts)) +
  geom_point(size = 0.2) +  
  stat_cor(method = "pearson") +
  geom_abline(intercept = 0, slope = 1) +
  labs(title = "Gene count", x = "STAR + its_a_scam", y = "Parse") +
  scale_x_log10() +
  scale_y_log10() +
  theme_cowplot(10) 
ggsave("h1_compare_counts.pdf", height = 4, width = 4)

ggplot(data = h1.counts, aes(x = star_scam_features, y = parse_features)) +
  geom_point(size = 0.2) +  
  stat_cor(method = "pearson") +
  geom_abline(intercept = 0, slope = 1) +
  labs(title = "Number of genes", x = "STAR + its_a_scam", y = "Parse") +
  scale_x_log10() +
  scale_y_log10() +
  theme_cowplot(10) 
ggsave("h1_compare_featuress.pdf", height = 4, width = 4)




h2.parse <- LoadSeuratRds("../../../count_matrix/og/h2_og.Rds")
h2.star <- LoadSeuratRds("../../../count_matrix/h2/h2_n500_collapse.Rds")

h2.common <- intersect(colnames(h2.parse), colnames(h2.star))
length(h2.common)

h2.parse
h2.parse <- h2.parse[, colnames(h2.parse) %in% h2.common]
h2.star
h2.star <- h2.star[, colnames(h2.star) %in% h2.common]

h2.star.counts <- tibble(cell_id = colnames(h2.star), 
  star_scam_counts = h2.star$nCount_RNA, 
  star_scam_features = h2.star$nFeature_RNA)
h2.star.counts
h2.parse.counts <- tibble(cell_id = colnames(h2.parse), 
  parse_counts = h2.parse$nCount_RNA, 
  parse_features = h2.parse$nFeature_RNA)
h2.parse.counts

h2.counts <- h2.star.counts %>% 
  left_join(h2.parse.counts)

ggplot(data = h2.counts, aes(x = star_scam_counts, y = parse_counts)) +
  geom_point(size = 0.2) +  
  stat_cor(method = "pearson") +
  geom_abline(intercept = 0, slope = 1) +
  labs(title = "Gene count", x = "STAR + its_a_scam", y = "Parse") +
  scale_x_log10() +
  scale_y_log10() +
  theme_cowplot(10) 
ggsave("h2_compare_counts.pdf", height = 4, width = 4)

ggplot(data = h2.counts, aes(x = star_scam_features, y = parse_features)) +
  geom_point(size = 0.2) +  
  stat_cor(method = "pearson") +
  geom_abline(intercept = 0, slope = 1) +
  labs(title = "Number of genes", x = "STAR + its_a_scam", y = "Parse") +
  scale_x_log10() +
  scale_y_log10() +
  theme_cowplot(10) 
ggsave("h2_compare_featuress.pdf", height = 4, width = 4)
