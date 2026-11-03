library(tidyverse)
library(cowplot)

a <- read_csv("genomic_step_vector_bench.csv")

e <- a %>% 
  filter(title == "GenomicStepVector::add  multi-chromosome random sorted overlapping") 

e$nchrom <- as.numeric(str_split_i(str_split_i(e$name, " ", 1), "=", 2))

ggplot(data = e) + 
  geom_boxplot(mapping = aes(x = as.factor(complexityN), 
    y = elapsed*10^6 / complexityN, color = as.factor(nchrom)), , outlier.size = 0.25) + 
    scale_y_continuous(transform = "log10") + 
    labs(x = "Number of elements", 
      y = "Time per element (us)", color = "Chroms") +
    theme_cowplot(10)
ggsave("multichrom_sorted_time.pdf", width = 4, height = 4)


e <- a %>% 
  filter(title == "GenomicStepVector::add  multi-chromosome random sorted overlapping 2") 

e$nchrom <- as.numeric(str_split_i(str_split_i(e$name, " ", 1), "=", 2))
e$nperchrom <- as.numeric(str_split_i(str_split_i(e$name, " ", 2), "=", 2))

ggplot(data = e) + 
  geom_boxplot(mapping = aes(x = as.factor(nperchrom), 
    y = elapsed*10^6 / complexityN, color = as.factor(nchrom)), , outlier.size = 0.25) + 
    scale_y_continuous(transform = "log10") + 
    labs(x = "Number of elements per chrom", 
      y = "Time per element (us)", color = "Chroms") +
    theme_cowplot(10)
ggsave("multichrom_sorted_2_time.pdf", width = 4, height = 4)
