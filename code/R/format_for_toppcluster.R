library(plyr)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(reshape2)

library(synapseClient)
synapseLogin()

obj <- synGet("syn4216912")

d <- read_tsv(getFileLocation(obj))

tolong <- function(x) {
  contr <- unique(x$Contrast.Names)
  genes <- str_split(x$Gene.Symbols, ",")
  d <- data.frame(gene=genes)
  colnames(d) <- c("gene")
  d
}

d2 <- ddply(d, .(Contrast.Names), tolong, .progress='text')

d2.up <- d2 %>% 
  filter(str_detect(Contrast.Names, "up"), !is.na(gene), gene != "") %>%
  select(gene, Contrast.Names)

d2.down <- d2 %>% 
  filter(str_detect(Contrast.Names, "down"), !is.na(gene), gene != "") %>%
  select(gene, Contrast.Names)

write.table(d2.up, file="DiffnameShort_DiffGenes_Express_mixedEffects_up_toppcluster.csv",
            sep=",", row.names=FALSE, col.names=FALSE, quote=FALSE)

write.table(d2.down, file="DiffnameShort_DiffGenes_Express_mixedEffects_down_toppcluster.csv",
            sep=",", row.names=FALSE, col.names=FALSE, quote=FALSE)
