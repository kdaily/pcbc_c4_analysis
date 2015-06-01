library(plyr)
library(dplyr)
library(tidyr)
library(data.table)
library(stringr)
library(reshape2)

library(synapseClient)
synapseLogin()

# Requires file output from ToppCluster
# These give warnings on many rows
dUpObj <- synGet("syn4231394")
dDownObj <- synGet("syn4231385")

# tmp <- read.delim(getFileLocation(dUpObj), sep="\t", nrows=200)
# Use to get column classes
# sapply(tmp, class)

colTypes <- c(rep("character", 4),
              rep("numeric", sum(str_detect(colnames(tmp), "logP"))),
              rep("character", sum(str_detect(colnames(tmp), "_GeneSet$"))))

d.up <- fread(getFileLocation(dUpObj), colClasses=colTypes, data.table=F)
d.down <- fread(getFileLocation(dDownObj), colClasses=colTypes, data.table=F)

# Clean up columns
colnames(d.up)[1:4] <- c("Category", "ID", "Title", "VerboseID")
colnames(d.down)[1:4] <- c("Category", "ID", "Title", "VerboseID")

idxcols <- c("Category", "ID", "Title", "VerboseID")

d2 <- d.up %>%
  full_join(d.down, by=idxcols) %>% # join up and down together
  dplyr::select(VerboseID, ends_with("logP")) %>% # only get id and value columns
  filter(!is.na(VerboseID)) # some NAs to get rid of

# Clean up column names again to make them easier to read
colnames(d2) <- str_replace(colnames(d2), "_logP", "")

# Use only the diff state columns
d3 <- d2 %>% select(-contains("ale"), -contains("virus"), 
                    -contains("plasmid"), -contains("NotApplicable")) 

# Matrix for clustering
forclust <- as.data.frame(d3[, -1])
rownames(forclust) <- d3$VerboseID
forclust[is.na(forclust)] <- 0

# Only keep categories with a minimum number of significant differences
minDifferences <- 2
keep <- apply(forclust, 1, function(x) sum(x > 0) > minDifferences)
forclust2 <- forclust[keep, ]

# standard hierarchical clustering with euclidean distances
d.dist <- dist(forclust2)
cluster.cols <- hclust(d.dist)

# Not even useful
pdf("DiffState_ToppCluster_cluster.pdf", width=100, height=20)
plot(cluster.cols)
dev.off()

plot.data <- data.frame(forclust2, VerboseID=rownames(forclustBinary))
plot.data <- melt(plot.data, id.vars="VerboseID")
plot.data$direction <- "up"
plot.data$direction[str_detect(plot.data$variable, "down")] <- "down"

plot.data$VerboseID <- factor(plot.data$VerboseID, levels=cluster.cols$labels, ordered=TRUE)
plot.data$VerboseID2 <- as.numeric(plot.data$VerboseID)
p <- ggplot(plot.data, aes(x=VerboseID, y=variable, fill=value)) + geom_tile()
p <- p + facet_grid(direction ~ ., scale="free_y")
p <- p + scale_fill_gradient(low='white', high='black')
p <- p + theme(axis.text.x=element_blank())
p
