## Merge GO-Elite results

library(plyr)
library(dplyr)
library(tidyr)
library(data.table)
library(synapseClient)

synapseLogin()

## Input files have been pre-processed to remove headers
inputDir <- "/external-data/DAT_114__PCBC_Data/mRNA/GO-Elite/GO_Elite_output/GO-Elite_results/CompleteResults/ORA/archived-20150520-045614/processed"

# Get list of files
fileList <- list.files(path = inputDir, pattern = '*.txt')
nameList <- gsub("\\.txt", "", fileList)

# Read them in
dataList <- mlply(fileList, 
                  function(x) fread(paste(inputDir, x, sep="/"), 
                                    data.table=FALSE))

# Convert to DF
dataDF <- ldply(dataList, I)

# Names never propogate for some reason...
dataDF$X1 <- nameList[dataDF$X1]

# Process them
dataDF2 <- dataDF %>%
  rename(comparison=X1) %>%
  separate(comparison, c("contrast", "geneset"), sep="-")

# Column names have spaces, synapse tables doesn't like it
colnames(dataDF2) <- gsub("[ -]", "_", colnames(dataDF2))

# Auto-generate table schema
tcresult <- as.tableColumns(dataDF2)
cols <- tcresult$tableColumns
fileHandleId <- tcresult$fileHandleId

# Need longer string cols
cols[[1]]@maximumSize <- as.integer(200)
cols[[2]]@maximumSize <- as.integer(200)
cols[[3]]@maximumSize <- as.integer(200)

# Store in Synapse
projectId <- "syn1773109"
schema <- TableSchema(name="GOElite_eXpress", parent=projectId, columns=cols)
table <- Table(schema, values=dataDF2)
table <- synStore(table)
