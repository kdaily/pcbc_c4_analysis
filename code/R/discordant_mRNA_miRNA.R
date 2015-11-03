setwd("~/src/discordant/")
source("~/src/discordant/discordant.R")

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
mrna_obj <- synGet("syn5011097")

mrna_mat <- fread(getFileLocation(mrna_obj), data.table=FALSE)

mirna_obj <- synGet("syn5014456")
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

mrna_mat <- mrna_mat %>%
  dplyr::select(GeneName, one_of(mrna_metadata_filtered$UID))

mirna_mat <- mirna_mat %>%
  dplyr::select(GeneName, one_of(mirna_metadata_filtered$UID))

# Take the median across multiple biological samples per feature
mrna_mat_median <- melt(mrna_mat, id.vars="GeneName", 
                        variable.name = "UID", value.name = "expression") %>% 
  left_join(mrna_metadata_filtered %>% select(UID, biologicalSampleName)) %>% 
  group_by(GeneName, biologicalSampleName) %>% 
  summarize(median_expression=median(expression)) %>%
  dcast(GeneName ~ biologicalSampleName)

mirna_mat_median <- melt(mirna_mat, id.vars="GeneName", 
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
