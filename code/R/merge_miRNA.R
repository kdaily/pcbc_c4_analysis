# Merge miRNA files from mirExpress

library(data.table)
library(plyr)
library(dplyr)
library(reshape2)
library(synapseClient)
library(rGithubClient)
synapseLogin()

# Query for the files
q <- "select name,id,UID FROM file WHERE projectId=='syn1773109' AND dataType=='miRNA' AND fileType=='expr'"
qr <- synQuery(q, blockSize=500)
res <- qr$collectAll()
synIds <- res$file.id

# Get the objects and download if necessary
objs <- mlply(synIds, synGet, downloadFile=TRUE, .progress='text')
names(objs) <- lapply(objs, function(x) annotations(x)$UID)

# Get the file names
listnames <- llply(objs, getFileLocation, .progress='text')

# Read all the files in, only specific columns
# Merge them into a long data frame
pp1 <- ldply(listnames, .id=NA, .progress='text',
             fread, header=FALSE, sep="\t", 
             data.table=FALSE)

pp1$UID <- names(listnames)[pp1$X1]
pp1$X1 <- NULL
colnames(pp1) <- c("target_id", "count", "UID")
pp1 <- pp1[, c("UID", "target_id", "count")]

## Cast into wide data frames, one per measurement
counts <- dcast(pp1, target_id ~ UID, value.var="count")
write.csv(estCounts, "miRNA_counts.csv", row.names=FALSE)
countsFile <- File("miRNA_counts.csv", parentId="syn3354743",
                      annotations=list(fileType="matrix", dataType="miRNA"))
countsFile <- synStore(countsFile, used=synIds)
