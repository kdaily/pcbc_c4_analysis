# Merge miRNA files from mirExpress

library(data.table)
library(plyr)
library(dplyr)
library(reshape2)
library(synapseClient)
library(rGithubClient)
synapseLogin()

# Get the metadata
metaSchema <- synGet('syn3156503')
metaTbl <- synTableQuery(paste("SELECT * FROM", metaSchema@properties$id))
metadata <- metaTbl@values

# Query for the files
q <- "select id,UID FROM file WHERE benefactorId=='syn1773109' AND dataType=='mRNA' AND fileType=='fpkm' AND fileSubType=='genes'"
qr <- synQuery(q, blockSize=500)
res <- qr$collectAll()

q2 <- "select id,UID,dataType,fileType,fileSubType FROM file WHERE parentId=='syn2246521'"
qr2 <- synQuery(q2, blockSize=200)
res2 <- qr2$collectAll()

metadataPublic <- metadata %>% filter(public)
resPublic <- res %>% filter(file.UID %in% metadataPublic$UID)
synIds <- resPublic$file.id

metadata %>% 
  filter(!(UID %in% res2$file.UID)) %>% 
  filter(pass_qc, !exclude, public) %>% 
  dplyr::select(UID, pass_qc, exclude, public, Diffname_short)

# Get the objects and download if necessary
objs <- mlply(synIds, synGet, downloadFile=TRUE, .progress='text')
names(objs) <- lapply(objs, function(x) annotations(x)$UID)

# Get the file names
listnames <- llply(objs, getFileLocation, .progress='text')

# Read all the files in, merge them into a long data frame
pp1 <- ldply(listnames, .id=NA, .progress='text',
             fread, header=TRUE, sep="\t", 
             select=c("tracking_id", "FPKM"),
             data.table=FALSE, colClasses="character")

pp1 <- pp1 %>% rename(UID=X1) %>%
  mutate(UID=names(listnames)[UID])

## Cast into wide data frames, one per measurement
fpkms <- dcast(pp1, tracking_id ~ UID, value.var="FPKM")

  write.table(fpkms, "public_expressionCalls_summarized.tsv",
              row.names=FALSE, quote=FALSE, sep="\t")

repo <- getRepo("kdaily/pcbc_c4_analysis", ref="branch", refName="mergemrna")
script <- getPermlink(repo, repositoryPath="code/R/merge_cufflinks_public.R")

fpkmFile <- File("public_expressionCalls_summarized.tsv", 
                 name="mRNA Expression Matrix", parentId="syn3219792",
                 annotations=list(fileType="genomicMatrix", dataType="mRNA"))

fpkmFile <- synStore(countsFile, used="syn2246521", executed=script)
