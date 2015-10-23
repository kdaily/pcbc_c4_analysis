## Collapse eXpress transcript level values to gene level
library(rGithubClient)
library(synapseClient)
library(data.table)
library(dplyr)

synapseLogin()

collapse_eXpress <- function(dataID, mappingData) {
  expressFile <- synGet(dataID)
  expressData <- fread(getFileLocation(expressFile), data.table=FALSE)
  
  expressData %>% 
    left_join(mappingData) %>%
    dplyr::select(-1) %>%
    group_by(hg19.kgXref.geneSymbol) %>%
    summarise_each(funs(sum))
}

## Get the UCSC id to gene symbol mapping file
mappingFile <- synGet("syn3444900")
mappingData <- tbl_df(fread(getFileLocation(mappingFile), data.table=FALSE))
mappingData <- mappingData %>% rename(target_id=hg19.knownGene.name)

## Set up provenance
repo <- getRepo("Sage-Bionetworks/pcbc_c4_analysis")
thisFile <- getPermlink(repo, "code/R/eXpress_transcripts_to_gene.R")

## Set up annotations
annots <- list(dataType="mRNA", fileType="genomicMatrix", expressionLevel="gene")

effCountID <- "syn5006129"
effCountData <- collapse_eXpress(effCountID, mappingData)
write.csv(effCountData, "eXpress_eff_count_geneSymbol.csv", row.names=FALSE)
f <- File("eXpress_eff_count_geneSymbol.csv", parentId="syn3354743")
synSetAnnotations(f) <- annots
generatedBy(f) <- Activity(used=c(effCountID, mappingFile), 
                                    executed=thisFile,
                                    name="Sum", description="Sum over gene symbol")
f <- synStore(f)
