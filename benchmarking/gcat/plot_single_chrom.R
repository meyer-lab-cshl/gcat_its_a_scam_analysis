library(tidyverse)
library(cowplot)
a <- read_csv("genomic_step_vector_bench.csv")

c <- a %>% 
  filter(title == "GenomicStepVector::add  single chromosome random overlapping") %>% 
  mutate(type = "random")

d <- a %>% 
  filter(title == "GenomicStepVector::add  single chromosome random sorted overlapping") %>% 
  mutate(type = "sorted")

d <- rbind(d, c)

ggplot(data = d) + 
  geom_boxplot(mapping=aes(x = as.factor(complexityN), 
    y = elapsed*10^6 / complexityN, color = type), , outlier.size = 0.25) + 
  scale_y_continuous(transform = "log10") + 
  labs(x = "Number of elements", y = "Time per element (us)", 
  color = "add type") +
  theme_cowplot(10)
ggsave("single_chrom_time.pdf", height = 4, width = 4)

