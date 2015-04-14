# Merge miRNA files from mirExpress

library(data.table)
library(plyr)
library(dplyr)
library(reshape2)
library(synapseClient)
library(rGithubClient)
synapseLogin()

# Get the metadata
metaSchema <- synGet('syn3219876')
metaTbl <- synTableQuery(paste("SELECT * FROM", metaSchema@properties$id))
metadata <- metaTbl@values

# Query for the files
q <- "select id,UID FROM file WHERE projectId=='syn1773109' AND dataType=='miRNA' AND fileType=='expr'"
qr <- synQuery(q, blockSize=500)
res <- qr$collectAll()

metadataPublic <- metadata %>% filter(public)
resPublic <- res %>% filter(file.UID %in% metadataPublic$UID)
synIds <- resPublic$file.id

# Get the objects and download if necessary
objs <- mlply(synIds, synGet, downloadFile=TRUE, .progress='text')
names(objs) <- lapply(objs, function(x) annotations(x)$UID)

# Get the file names
listnames <- llply(objs, getFileLocation, .progress='text')

# Read all the files in, merge them into a long data frame
pp1 <- ldply(listnames, .id=NA, .progress='text',
             fread, header=FALSE, sep="\t", 
             data.table=FALSE)

pp1$UID <- names(listnames)[pp1$X1]
pp1$X1 <- NULL
colnames(pp1) <- c("mir", "count", "UID")
pp1 <- pp1[, c("UID", "mir", "count")]

## Cast into wide data frames, one per measurement
counts <- dcast(pp1, mir ~ UID, value.var="count")
write.csv(counts, "miRNA_counts.csv", row.names=FALSE)

repo <- getRepo("kdaily/pcbc_c4_analysis", ref="branch", refName="mergemirna")
script <- getPermlink(repo, repositoryPath="code/R/merge_miRNA_public.R")

countsFile <- File("miRNA_counts.csv", name="miRNA Count Matrix", parentId="syn3219792",
                   annotations=list(fileType="genomicMatrix", dataType="miRNA"))

countsFile <- synStore(countsFile, used="syn2247164", executed=script)
