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

## Set up provenance
repo <- getRepo("kdaily/pcbc_c4_analysis", ref="branch", refName="eXpress_synapse")
thisFile <- getPermlink(repo, "code/R/eXpress_transcripts_to_gene.R")

## Set up annotations
annots <- list(dataType="mRNA", fileType="matrix")

tpmID <- "syn3354728"
tpmData <- collapse_eXpress(tpmID, mappingData)
write.csv(tpmData, "eXpress_tpm_geneSymbol.csv", row.names=FALSE)
f <- File("eXpress_tpm_geneSymbol.csv", parentId="syn3354743")
synSetAnnotations(f) <- annots
generatedBy(f) <- Activity(used=c(tpmID, mappingFile), 
                           executed=thisFile,
                           name="Sum", description="Sum over gene symbol")
f <- synStore(f)

effCountID <- "syn3354716"
effCountData <- collapse_eXpress(effCountID, mappingData)
write.csv(effCountData, "eXpress_eff_count_geneSymbol.csv", row.names=FALSE)
f <- File("eXpress_eff_count_geneSymbol.csv", parentId="syn3354743")
synSetAnnotations(f) <- annots
generatedBy(resultFile) <- Activity(used=c(effCountID, mappingFile), 
                                    executed=thisFile,
                                    name="Sum", description="Sum over gene symbol")
f <- synStore(f)


estCountID <- "syn3354714"
estCountData <- collapse_eXpress(estCountID, mappingData)
write.csv(estCountData, "eXpress_est_count_geneSymbol.csv", row.names=FALSE)
f <- File("eXpress_est_count_geneSymbol.csv", parentId="syn3354743")
synSetAnnotations(f) <- annots
generatedBy(f) <- Activity(used=c(estCountID, mappingFile), 
                           executed=thisFile,
                           name="Sum", description="Sum over gene symbol")
f <- synStore(f)

fpkmID <- "syn3354736"
fpkmData <- collapse_eXpress(fpkmID, mappingData)
write.csv(fpkmData, "eXpress_fpkm_geneSymbol.csv", row.names=FALSE)
f <- File("eXpress_fpkm_geneSymbol.csv", parentId="syn3354743")
synSetAnnotations(f) <- annots
generatedBy(f) <- Activity(used=c(fpkmID, mappingFile), 
                           executed=thisFile,
                           name="Sum", description="Sum over gene symbol")
f <- synStore(f)