cwd <- getwd()
setwd("~/src/discordant/")
source("~/src/discordant/discordant.R")
setwd(cwd)

library(data.table)
library(plyr)
library(dplyr)
library(reshape2)
library(synapseClient)
synapseLogin()

# Get metadata
mrna_metdata_id <- "syn3156503"
mrna_metadata_obj <- synTableQuery(sprintf('SELECT * FROM %s', mrna_metdata_id))
mrna_metadata <- mrna_metadata_obj@values

mrna_metadata[mrna_metadata == 'N/A'] = NA

mirna_metdata_id <- "syn3219876"
mirna_metadata_obj <- synTableQuery(sprintf('SELECT * FROM %s', mirna_metdata_id))
mirna_metadata <- mirna_metadata_obj@values

mirna_metadata[mirna_metadata == 'N/A'] = NA

# Get mRNA data
mrna_id <- "syn5011097"
mrna_obj <- synGet(mrna_id)

mrna_mat <- fread(getFileLocation(mrna_obj), data.table=FALSE)

mirna_id <- "syn5014456"
mirna_obj <- synGet(mirna_id)
mirna_mat <- fread(getFileLocation(mirna_obj), data.table=FALSE)

# Filter metadata - will use UIDs from these to subset matrices
mrna_metadata_filtered <- mrna_metadata %>%
  filter(public, pass_qc, !exclude,
         UID %in% colnames(mrna_mat),
         Diffname_short != "",
         Cell_Type == "PSC",
         C4_Karyotype_Result != "abnormal")

mirna_metadata_filtered <- mirna_metadata %>%
  filter(public, pass_qc, !exclude,
         UID %in% colnames(mirna_mat),
         Diffname_short != "",
         Cell_Type == "PSC",
         C4_Karyotype_Result != "abnormal")

# Need to match up by biological sample between assays
biosampleInBoth <- intersect(mrna_metadata_filtered$biologicalSampleName, 
                             mirna_metadata_filtered$biologicalSampleName)

mrna_metadata_filtered <- mrna_metadata_filtered %>%
  filter(biologicalSampleName %in% biosampleInBoth)

mirna_metadata_filtered <- mirna_metadata_filtered %>%
  filter(biologicalSampleName %in% biosampleInBoth)


# Get differentially expressed transcripts and microRNAs ------------------
comparisons <- synTableQuery("SELECT * FROM syn4483642")@values

diffexpr_mrna_id <- "syn5013690"
DE.mRNA <- fread(getFileLocation(synGet(diffexpr_mrna_id)), data.table=FALSE) %>%
  mutate(Comparison=str_replace(Comparison, "_vs_", "_")) %>% 
  tidyr::separate(Comparison, into=c("dataRestriction", "variable1short", 
                                     "variable2short", "direction"), sep = "_+") %>%
  filter(logFC >= 2, adj.P.value <= 0.05, 
         variable1short %in% c("SC", "DE"),
         variable2short %in% c("SC", "DE"))

diffexpr_mirna_id <- "syn5014584"
DE.miRNA <- fread(getFileLocation(synGet(diffexpr_mirna_id)), data.table=FALSE) %>%
  mutate(Comparison=str_replace(Comparison, "_vs_", "_")) %>% 
  tidyr::separate(Comparison, into=c("dataRestriction", "variable1short", 
                                     "variable2short", "direction"), sep = "_+") %>%
  filter(logFC >= 2, adj.P.value <= 0.05, 
         variable1short %in% c("SC", "DE"),
         variable2short %in% c("SC", "DE"))

mrna_mat_use <- mrna_mat %>%
  dplyr::filter(GeneName %in% DE.mRNA$GeneSymbol) %>% 
  dplyr::select(GeneName, one_of(mrna_metadata_filtered$UID))

mirna_mat_use <- mirna_mat %>%
  dplyr::filter(GeneName %in% DE.miRNA$GeneSymbol) %>% 
  dplyr::select(GeneName, one_of(mirna_metadata_filtered$UID))

# Take the median across multiple biological samples per feature
mrna_mat_median <- melt(mrna_mat_use, id.vars="GeneName", 
                        variable.name = "UID", value.name = "expression") %>% 
  left_join(mrna_metadata_filtered %>% select(UID, biologicalSampleName)) %>% 
  group_by(GeneName, biologicalSampleName) %>% 
  summarize(median_expression=median(expression)) %>%
  dcast(GeneName ~ biologicalSampleName)

mirna_mat_median <- melt(mirna_mat_use, id.vars="GeneName", 
                         variable.name = "UID", value.name = "expression") %>% 
  left_join(mirna_metadata_filtered %>% select(UID, biologicalSampleName)) %>% 
  group_by(GeneName, biologicalSampleName) %>% 
  summarize(median_expression=median(expression)) %>% 
  dcast(GeneName ~ biologicalSampleName)

# Use diffstate as a test case - get metadata that is either SC or DE
mrna_metadata_use <- mrna_metadata_filtered %>%
  filter(Diffname_short %in% c("SC", "DE")) %>% 
  select(biologicalSampleName, Diffname_short)

mirna_metadata_use <- mirna_metadata_filtered %>%
  filter(Diffname_short %in% c("SC", "DE")) %>% 
  select(biologicalSampleName, Diffname_short)

metadata_use <- unique(rbind(mrna_metadata_use, mirna_metadata_use))

# get only the columns of data, no feature names
mrna_mat_median_use <- mrna_mat_median %>%
  select(one_of(metadata_use$biologicalSampleName))

mirna_mat_median_use <- mirna_mat_median %>%
  select(one_of(metadata_use$biologicalSampleName))

# Make the groups
refGroup <- "SC"
groups <- (metadata_use$Diffname_short == refGroup) + 1

# Run discordant on two omics datasets
vectors <- createVectors(mrna_mat_median_use, mirna_mat_median_use, groups = groups)

result <- discordantRun(vectors$v1, vectors$v2, 
                        mrna_mat_median_use, mirna_mat_median_use)

resultTable <- makeTable(result$discordPPMatrix,
                         mrna_mat_median_use, mirna_mat_median_use) 

resultTableFixed <- resultTable %>% 
  as.data.frame() %>% 
  rename(feature1=V1, feature2=featureRows, value=V3) %>% 
  mutate(value=as.numeric(as.character(value)),
         feature1=as.numeric(as.character(feature1)),
         feature2=as.numeric(as.character(feature2)),
         transcript=mrna_mat_median$GeneName[feature1],
         mirna=mirna_mat_median$GeneName[feature2]
         )

resultTableFiltered <- resultTableFixed %>% filter(value >= 0.95)

save(vectors, result, resultTable, resultTableFixed, 
     mrna_mat_median_use, mirna_mat_median_use,
     file="discordant_mRNA_miRNA__DiffExpr_DE_SC.Rdata")

folder <- synStore(Folder("Discordant", parentId="syn3237603"))
w <- synGetWiki(folder)

repo <- rGithubClient::getRepo("kdaily/pcbc_c4_analysis", ref="branch", refName="discordant_de_sc")
thisScript <- rGithubClient::getPermlink(repo, repositoryPath="code/R/discordant_mRNA_miRNA.R")

o <- File("discordant_mRNA_miRNA__DE_SC.Rdata", parentId=synGetProperty(folder, "id"))
annotations(o) <- list(fileType="RData", dataType=c("mRNA", "miRNA"), comparison="All__DE_vs_SC")
o <- synStore(o, used=c(mrna_metdata_id, mirna_metdata_id, mrna_id, mirna_id,
                        diffexpr_mrna_id, diffexpr_mirna_id),
              executed=thisScript)
