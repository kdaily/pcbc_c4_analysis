library(data.table)
library(plyr)
library(reshape2)
library(synapseClient)
library(rGithubClient)
synapseLogin()

# Query for the files
q <- "select id from file where parentId=='syn3270268' AND dataType=='mRNA' AND fileType=='eXpress'"
res  <- synQuery(q)
synIds <- res$file.id

# Get the objects and download if necessary
objs <- mlply(synIds, synGet, downloadFile=TRUE, .progress='text')
names(objs) <- lapply(objs, function(x) annotations(x)$UID)

# Get the file names
listnames <- llply(objs, getFileLocation, .progress='text')

# Read all the files in, only specific columns
# Merge them into a long data frame
pp1 <- ldply(listnames, .id=NA, .progress='text', .parallel=FALSE,
             fread, select=c("target_id", "est_counts",
                             "fpkm", "eff_counts", "tpm"),
             data.table=FALSE)

pp1$UID <- names(listnames)[pp1$X1]
pp1$X1 <- NULL
pp1 <- pp1[, c("UID", "target_id", "est_counts", 
               "fpkm", "eff_counts", "tpm")]

write.csv(pp1, "eXpress_merged.csv", row.names=FALSE)
mergedFile <- File("eXpress_merged.csv", parentId="syn3354743",
                   annotations=list(fileType="csv", dataType="mRNA", 
                                    expressionLevel="transcript"))
mergedFile <- synStore(mergedFile, used="syn3270268")

## Cast into wide data frames, one per measurement

estCounts <- dcast(pp1, target_id ~ UID, value.var="est_counts")
write.csv(estCounts, "eXpress_est_counts.csv", row.names=FALSE)
estCountsFile <- File("eXpress_est_counts.csv", parentId="syn3354743",
                      annotations=list(fileType="genomicMatrix", dataType="mRNA", expressionLevel="transcript"))
estCountsFile <- synStore(estCountsFile, used=mergedFile@properties$id)

effCounts <- dcast(pp1, target_id ~ UID, value.var="eff_counts")
write.csv(effCounts, "eXpress_eff_counts.csv", row.names=FALSE)
effCountsFile <- File("eXpress_eff_counts.csv", parentId="syn3354743",
                      annotations=list(fileType="genomicMatrix", dataType="mRNA", expressionLevel="transcript"))
effCountsFile <- synStore(effCountsFile, used=mergedFile@properties$id)

tpm <- dcast(pp1, target_id ~ UID, value.var="tpm")
write.csv(tpm, "eXpress_tpm.csv", row.names=FALSE)
tpmFile <- File("eXpress_tpm.csv", parentId="syn3354743",
                annotations=list(fileType="genomicMatrix", dataType="mRNA", expressionLevel="transcript"))
tpmFile <- synStore(tpmFile, used=mergedFile@properties$id)

fpkm <- dcast(pp1, target_id ~ UID, value.var="fpkm")
write.csv(fpkm, "eXpress_fpkm.csv", row.names=FALSE)
fpkmFile <- File("eXpress_fpkm.csv", parentId="syn3354743",
                 annotations=list(fileType="genomicMatrix", dataType="mRNA", expressionLevel="transcript"))
fpkmFile <- synStore(fpkmFile, used=mergedFile@properties$id)
