library(plyr)
library(dplyr)
library(tidyr)
library(data.table)
library(reshape2)

infile <- '/home/kdaily/Projects/PCBC/pcbc_c4_analysis/code/pipelines/GOElite/output/final/results.tsv'

d <- fread(infile, data.table=FALSE)

# Clean up column names
colnames(d) <- gsub("[- ]", "_", colnames(d))

# Filter, select columns
d2 <- d %>%
  filter(Z_Score > 0, Number_Changed > 4, PermuteP < 0.05) %>%
  select(Gene_Set_Name, contrast, geneset, symbol, AdjustedP) %>%
  tidyr::unite(group, geneset, Gene_Set_Name, contrast, sep="_", remove = FALSE)

# cast into a wide matrix
d3 <- d2 %>% 
  select(contrast, geneset, Gene_Set_Name, AdjustedP) %>% 
  unique() %>% 
  dcast(., geneset + Gene_Set_Name ~ contrast, value.var = "AdjustedP")

# Any NA's get a p-value of 1
d3[is.na(d3)] <- 1

# Convert to log_10 for clustering
d3[, -(1:2)] <- -log10(d3[, -(1:2)])

cluster.cols <- hclust(dist(t(d3[, -(1:2)])))
plot(cluster.cols)

tofile <- function(x, dir) {
  group <- unique(x$group)
  group <- gsub("[]")
  fn <- paste0(group, ".tsv")
  write.table(x, file=fn, sep="\t", row.names=F, quote=F)
}

d_ply(d2, .(group), tofile, dir="/home/kdaily/Projects/PCBC/pcbc_c4_analysis/code/pipelines/GOElite/output/toppcluster/")
