library(tidyverse)
library(cowplot)


a <- read_csv("genomic_step_vector_bench.csv")

b <- a %>% 
  filter(title == "GenomicStepVector::at single-chromosome unsorted query") %>% 
  mutate(type = "Fixed entries (1M); random query")
b$time_per_query <- b$elapsed*10^6 / b$complexityN

c <- a %>% 
  filter(title == "GenomicStepVector::at single-chromosome sorted query") %>% 
  mutate(type = "Fixed entries (1M); sorted query")
c$time_per_query <- c$elapsed*10^6 / b$complexityN

c <- rbind(b, c)

ggplot(data = c) + 
  geom_boxplot(mapping=aes(x = as.factor(complexityN), 
    y = time_per_query, color = type), outlier.size = 0.25) + 
  scale_y_continuous(transform = "log10") + 
  labs(x = "Number of queries", y = "Time per query (us)", 
  color = "Query type") +
  theme_cowplot(10)
ggsave("single_chrom_at_fixed_entry.pdf", height = 4, width = 6)


d <- a %>% 
  filter(title == "GenomicStepVector::at single-chromosome varying entries") %>% 
  mutate(type = "Fixed sorted queries (1M)")
d$time_per_query <- d$elapsed*10^6 / 1000000

ggplot(data = d) + 
  geom_boxplot(mapping=aes(x = as.factor(complexityN), 
    y = time_per_query, color = type), outlier.size = 0.25) + 
  scale_y_continuous(transform = "log10") + 
  labs(x = "Number of entries", y = "Time per query (us)", 
  color = "Query type") +
  theme_cowplot(10)
ggsave("single_chrom_at_fixed_queries.pdf", height = 4, width = 6)
