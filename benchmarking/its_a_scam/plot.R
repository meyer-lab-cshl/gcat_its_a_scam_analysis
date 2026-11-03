library(tidyverse)
library(cowplot)
library(hms)

t <- read_tsv("tss_cov_time.txt", col_names = FALSE)
names(t) <- c("program", "type", "time") 
t$program <- factor(t$program, levels = c("tss", "cov"))

t$time <- as_hms(paste0("00:", t$time))
ggplot(data = t) + 
  geom_boxplot(mapping = aes(x = program, y = time, color = type),
    outlier.size = 0.25) + 
  scale_y_time(labels = function(l) strftime(l, '%M:%S')) +
  labs(x = "Program", y = "Runtime (mm:ss)", 
    color = "Type") +
  theme_cowplot(10)
ggsave("tss_cov_runtime.pdf", height = 4, width = 4)


m <- read_tsv("tss_cov_mem.txt", col_names = FALSE)
names(m) <- c("program", "type", "memory") 
m$program <- factor(m$program, levels = c("tss", "cov"))
m$memory <- m$memory / 1024

ggplot(data = m) + 
  geom_boxplot(mapping = aes(x = program, y = memory, color = type),
    outlier.size = 0.25) + 
  labs(x = "Program", y = "Maximum memory (MB)", 
    color = "Type") +
  theme_cowplot(10)
ggsave("tss_cov_memory.pdf", height = 4, width = 4)
