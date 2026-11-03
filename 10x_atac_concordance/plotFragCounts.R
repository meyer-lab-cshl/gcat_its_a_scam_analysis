library(Seurat)
library(Signac)
library(ggpubr)
library(cowplot)
library(tidyverse)

# compute the number of fragmets per barcode
bwa.frag.count <- CountFragments("../peaks/pbmc_granulocyte_sorted_3k_noM_fragments.tsv.gz")
print(dim(bwa.frag.count))
print(head(bwa.frag.count))

# get a list of barcodes that meet the min frags and their frag counts
bwa.counts <- as_tibble(bwa.frag.count[bwa.frag.count$frequency_count >= 1000, 
  c("CB", "frequency_count")])
names(bwa.counts) <- c("CB", "bwa_frequency_count")



# compute the number of fragmets per barcode
cr.frag.count <- CountFragments("../../../og/fragments/pbmc_granulocyte_sorted_3k_atac_fragments.tsv.gz")
print(dim(cr.frag.count))
print(head(cr.frag.count))

# get a list of barcodes that meet the min frags and their frag counts
cr.counts <- as_tibble(cr.frag.count[cr.frag.count$frequency_count >= 1000, 
  c("CB", "frequency_count")])
cr.counts$CB <- substr(cr.counts$CB, 1, 16)
names(cr.counts) <- c("CB", "cr_frequency_count")


# keep just the common barcoded
common.bc <- intersect(bwa.counts$CB, cr.counts$CB)
print(length(common.bc))
bwa.counts <- bwa.counts[bwa.counts$CB %in% common.bc, ]
cr.counts <- cr.counts[cr.counts$CB %in% common.bc, ]

# plot the frag counts
frag.counts <- bwa.counts %>% 
  left_join(cr.counts)

ggplot(data = frag.counts, aes(x = bwa_frequency_count, y = cr_frequency_count)) + 
  geom_point(size = 0.5) + 
  stat_cor(method = "pearson") +
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "BWA", y = "Cellranger") +
  theme_cowplot(10) 
ggsave("pbmc_granulocyte_sorted_3k_bwa_cr_raw_frag_count.pdf", height = 4, width = 4)
