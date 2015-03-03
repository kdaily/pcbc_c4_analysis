library(synapseClient)
library(plyr)

synapseLogin()

fastqDir <- "/external-data/DAT_114__PCBC_Data/miRNA/fastq/"
projectId <- "syn1773109"
queryString <- sprintf("select UID,id,parentId from file where projectId=='%s' and dataType=='miRNA' and fileType=='fastq'", projectId)
queryResults <- synQuery(queryString)
colnames(queryResults) <- gsub("file\\.", "", colnames(queryResults))

getFileInfo <- function(x) {
  o <- synGet(x$id, downloadFile=FALSE)
  
  data.frame(fileName=ifelse(is.null(o@fileHandle$fileName), NA, o@fileHandle$fileName), 
             filePath=ifelse(is.null(o@filePath), NA, o@filePath))
}

idWithName <- ddply(queryResults, .(id), getFileInfo, .progress="text")
idWithName <- transform(idWithName, filePath=gsub("file://", "", idWithName$filePath))

idWithNameToUpdate <- subset(idWithName, !is.na(filePath))
mergedResults <- merge(queryResults, idWithNameToUpdate, by="id")
write.csv(mergedResults, "mirna_upload.csv")

# ## BROKEN; WORKS IN PYTHON
# storeFile <- function(x, fastqDir) {
#   synFile <- synGet(x$id)
#   synFile@synapseStore <- TRUE
#   synFile@filePath <- x$filePath
#   synFile <- synStore(synFile)
# }
# 
# objs <- dlply(mergedResults, .(id), storeFile)
