## Collapse eXpress transcript level values to gene level

library(data.table)
library(dplyr)

collapse_eXpress <- function(dataID, mappingData) {
  expressFile <- synGet(dataID)
  expressData <- tbl_df(fread(getFileLocation(expressFile), data.table=FALSE))
  
  expressData %>% 
    left_join(mappingData) %>%
    select(-1) %>%
    group_by(hg19.kgXref.geneSymbol) %>%
    summarise_each(funs(sum))
}

## Get the UCSC id to gene symbol mapping file
mappingFile <- synGet("syn3444900")
mappingData <- tbl_df(fread(getFileLocation(mappingFile), data.table=FALSE))
mappingData <- mappingData %>% rename(target_id=hg19.knownGene.name)

tpmID <- "syn3354728"
tpmData <- collapse_eXpress(tpmID, mappingData)
write.csv(tpmData, "eXpress_tpm_geneSymbol.csv", row.names=FALSE)

effCountID <- "syn3354716"
effCountData <- collapse_eXpress(effCountID, mappingData)
write.csv(effCountData, "eXpress_eff_count_geneSymbol.csv", row.names=FALSE)

estCountID <- "syn3354714"
estCountData <- collapse_eXpress(estCountID, mappingData)
write.csv(estCountData, "eXpress_est_count_geneSymbol.csv", row.names=FALSE)

fpkmCountID <- "syn3354736"
fpkmCountData <- collapse_eXpress(fpkmCountID, mappingData)
write.csv(fpkmCountData, "eXpress_fpkm_geneSymbol.csv", row.names=FALSE)
