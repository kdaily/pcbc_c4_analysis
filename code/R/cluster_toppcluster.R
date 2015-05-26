library(plyr)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(reshape2)

library(synapseClient)
synapseLogin()

# Requires file output from ToppCluster
d.up <- read_tsv("/home/kdaily/Projects/PCBC/pcbc_c4_analysis/code/pipelines/GOElite/output/toppcluster/results/up.tsv")
d.down <- read_tsv("/home/kdaily/Projects/PCBC/pcbc_c4_analysis/code/pipelines/GOElite/output/toppcluster/results/down.tsv")

idxcols <- c("Category", "ID", "Title (or Source)")
d.idx <- d.up %>%
  select(one_of(idxcols))

d2 <- d.up %>% 
  full_join(d.down, by=idxcols) %>% 
  dplyr::select(1:3, ends_with("logP"))

d.dist <- dist(t(d2[, -(1:3)])) # euclidean distances between the rows
fit <- cmdscale(d.dist, eig=TRUE, k=2) # k is the number of dim
fit # view results

# plot solution

x <- fit$points[, 1]
y <- fit$points[, 2]
labels <- rownames(fit$points)

plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", 
     main="Metric MDS", type="n")

text(x, y, labels = labels, cex=.7)