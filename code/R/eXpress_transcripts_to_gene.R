library(data.table)
library(dplyr)

expressFile <- synGet("syn3354728")
expressData <- tbl_df(fread(getFileLocation(expressFile), data.table=FALSE))

mappingFile <- synGet("syn3444900")
mappingData <- tbl_df(fread(getFileLocation(mappingFile), data.table=FALSE))
mappingData <- mappingData %>% rename(target_id=hg19.knownGene.name)

results <- expressData %>% 
  left_join(mappingData) %>%
  select(-1) %>%
  group_by(hg19.kgXref.geneSymbol) %>%
  summarise_each(funs(sum))

write.csv(results, "eXpress_tpm_geneSymbol.csv", row.names=FALSE)


