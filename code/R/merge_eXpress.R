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
objs <- mlply(synIds, synGet, downloadFile=TRUE)
names(objs) <- lapply(objs, function(x) annotations(x)$UID)

# Get the file names
listnames <- llply(objs, getFileLocation)

# Read all the files in, only specific columns
# Merge them into a long data frame
pp1 <- ldply(listnames, fread, select=c("target_id", "est_counts", "fpkm", "eff_counts", "tpm"),
             data.table=FALSE, .id=NA)

write.csv(pp1, "eXpress_merged.csv", row.names=FALSE)
mergedFile <- File("eXpress_merged.csv", parentId="syn3270268",
                   annotations=list(fileType="csv", dataType="mRNA"))
mergedFile <- synStore(mergedFile)

## Cast into wide data frames, one per measurement

estCounts <- dcast(pp1, target_id ~ .id, value.var="est_counts")
write.csv(estCounts, "eXpress_est_counts.csv", row.names=FALSE)
estCountsFile <- File("eXpress_est_counts.csv", parentId="syn3270268",
                      annotations=list(fileType="csv", dataType="mRNA"))
estCountsFile <- synStore(estCountsFile)

effCounts <- dcast(pp1, target_id ~ .id, value.var="eff_counts")
write.csv(effCounts, "eXpress_eff_counts.csv", row.names=FALSE)
effCountsFile <- File("eXpress_eff_counts.csv", parentId="syn3270268",
                      annotations=list(fileType="csv", dataType="mRNA"))
effCountsFile <- synStore(effCountsFile)

tpm <- dcast(pp1, target_id ~ .id, value.var="tpm")
write.csv(tpm, "eXpress_tpm.csv", row.names=FALSE)
tpmFile <- File("eXpress_tpm.csv", parentId="syn3270268",
                annotations=list(fileType="csv", dataType="mRNA"))
tpmFile <- synStore(tpmFile)

fpkm <- dcast(pp1, target_id ~ .id, value.var="fpkm")
write.csv(fpkm, "eXpress_fpkm.csv", row.names=FALSE)
fpkmFile <- File("eXpress_fpkm.csv", parentId="syn3270268",
                annotations=list(fileType="csv", dataType="mRNA"))
fpkmFile <- synStore(fpkmFile)
