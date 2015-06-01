## GO-Elite results to Synapse table

library(plyr)
library(dplyr)
library(tidyr)
library(data.table)
library(synapseClient)
library(readr)
synapseLogin()

## Input files have been pre-processed in a pipeline
## https://github.com/kdaily/pcbc_c4_analysis/tree/goelite/code/pipelines/GOElite
id <- "syn4228803"
obj <- synGet(id)

d <- as.data.frame(read_tsv(getFileLocation(obj)))

# Column names have spaces, synapse tables doesn't like it
colnames(d) <- gsub("[ -]", "_", colnames(d))

# Auto-generate table schema
tcresult <- as.tableColumns(d)
cols <- tcresult$tableColumns
fileHandleId <- tcresult$fileHandleId

# Need longer string cols
cols[[1]]@maximumSize <- as.integer(200)
cols[[10]]@maximumSize <- as.integer(200)
cols[[11]]@maximumSize <- as.integer(200)

cols[[5]]@columnType <- "DOUBLE"
cols[[7]]@columnType <- "DOUBLE"
cols[[8]]@columnType <- "DOUBLE"
cols[[9]]@columnType <- "DOUBLE"

# Store in Synapse
projectId <- "syn1773109"
schema <- TableSchema(name="GOElite_eXpress", parent=projectId, columns=cols)
table <- Table(schema, values=d)
table <- synStore(table)
