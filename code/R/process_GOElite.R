library(plyr)
library(dplyr)
library(tidyr)
library(data.table)
library(reshape2)
library(synapseClient)
library(rGithubClient)

synapseLogin()

thisRepo <- getRepo("kdaily/pcbc_c4_analysis", ref="branch", refName="goelite")
thisScript <- getPermlink(repo=thisRepo, repositoryPath="code/R/process_GOElite.R")

inObj <- synGet('syn4228803')

d <- fread(getFileLocation(inObj), data.table=FALSE)

# Clean up column names
colnames(d) <- gsub("[- ]", "_", colnames(d))

# Filter, select columns
d2 <- d %>%
  filter(Z_Score > 0, Number_Changed > 4, PermuteP < 0.05) %>%
  select(Gene_Set_Name, contrast, geneset, Z_Score, AdjustedP) %>%
  tidyr::unite(GS, Gene_Set_Name, geneset, sep="_")

d3 <- d2 %>% 
  dcast(GS ~ contrast, value.var="Z_Score")%>%
  select(-contains("ale"), -contains("plasmid"), -contains("virus"), 
         -contains("plasmid"), -contains("NotApplicable"))

colnames(d3)

# Any NA's get a p-value of 1
d3[is.na(d3)] <- 1

# Convert to log_10 for clustering
forclust <- -log10(d3[, -1])
rownames(forclust) <- d3$GS

cluster.cols <- hclust(dist(forclust), method = "ward.D2")

pdf("DiffState_GOElite_cluster.pdf", width=100, height=20)
plot(cluster.cols)
dev.off()

f <- File("DiffState_GOElite_cluster.pdf", parentId='syn4228846')
f <- synStore(f, used=inObj@properties$id, executed=thisScript)
