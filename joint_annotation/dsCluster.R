library(tidyverse)
id <- read_tsv("pbmc3k_joint_celltype.tsv")
set.seed(42)
id.ds <- id[sample.int(nrow(id), 100),]
write_tsv(id.ds, "pbmc3k_joint_celltype_ds.tsv", col_names = FALSE)
