## Get all the mRNA samples with original filenames

library(synapseClient)
library(plyr)

synapseLogin()

newq <- "select id,UID from file where projectId=='syn1773109' AND dataType=='mRNA' AND fileType=='bam'"
newRes <- synapseQuery(newq)
colnames(newRes) <- gsub(".*\\.", "", colnames(newRes))

getFileName <- function(x) {
  o <- synGet(x$id, downloadFile=FALSE)
  data.frame(fileName=o@fileHandle$fileName)
}

foo <- ddply(newRes, .(id), getFileName, .progress="text")
newResMerged <- merge(newRes, foo, by="id")
write.csv(newResMerged, file="sample_table.csv", row.names=FALSE, quote=FALSE)
