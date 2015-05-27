library(plyr)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(reshape2)
library(rGithubClient)

library(synapseClient)
synapseLogin()

# Read in results of differential expression
obj <- synGet("syn4216912")

d <- read_tsv(getFileLocation(obj))

# Separate out the gene names for each contrast
tolong <- function(x) {
  contr <- unique(x$Contrast.Names)
  genes <- str_split(x$Gene.Symbols, ",")
  d <- data.frame(gene=genes)
  colnames(d) <- c("gene")
  d
}

d2 <- ddply(d, .(Contrast.Names), tolong, .progress='text')

# Split between up and down-regulated
d2up <- d2 %>% 
  filter(str_detect(Contrast.Names, "up"), !is.na(gene), gene != "") %>%
  select(gene, Contrast.Names)

d2down <- d2 %>% 
  filter(str_detect(Contrast.Names, "down"), !is.na(gene), gene != "") %>%
  select(gene, Contrast.Names)

write.table(d2up, file="DiffnameShort_DiffGenes_Express_mixedEffects_up_toppcluster.csv",
            sep=",", row.names=FALSE, col.names=FALSE, quote=FALSE)

repo <- getRepo("kdaily/pcbc_c4_analysis", ref="branch", refName="goelite")
thisScript <- getPermlink(repo, repositoryPath="code/R/format_for_toppcluster.R")

fUp <- File("DiffnameShort_DiffGenes_Express_mixedEffects_up_toppcluster.csv", parentId='syn3471223')
fUp <- synStore(fUp, used=obj@properties$id, executed=thisScript)

write.table(d2down, file="DiffnameShort_DiffGenes_Express_mixedEffects_down_toppcluster.csv",
            sep=",", row.names=FALSE, col.names=FALSE, quote=FALSE)

fDown <- File("DiffnameShort_DiffGenes_Express_mixedEffects_up_toppcluster.csv", parentId='syn3471223')
fDown <- synStore(fDown, used=obj@properties$id, executed=thisScript)
