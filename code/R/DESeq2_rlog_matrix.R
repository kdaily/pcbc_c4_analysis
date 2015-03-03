# Compute and write out DESeq processed regularized log matrix for mRNA-Seq data

library(synapseClient)
library(plyr)
library(dplyr)
library(DESeq2)
library(BiocParallel)
register(MulticoreParam(8))

synapseLogin()

countMatFile <- synGet("syn3164570")
countMat <- read.delim(countMatFile@filePath, row.names="Geneid", check.names=FALSE)
countMat <- as.matrix(countMat)

metadataTable <- synGet("syn3156503")
q <- sprintf("SELECT * FROM %s", metadataTable@properties$id)

metadata <- synTableQuery(q)@values

metadataFinal <- filter(metadata, UID %in% colnames(countMat))

rownames(metadataFinal) <- metadataFinal$UID
countMat <- countMat[, metadataFinal$UID]

## Read the data. The design formula influences the normalization, as dispersion
## estimates are performed within groups. However treating each sample
## separately should do something sensible for downstream analyses. According to
## the estimateDispersions help:
# "...substitute a design formula ~ 1 for the purposes of dispersion estimation.
# This treats the samples as replicates for the purpose of dispersion
# estimation. As mentioned in the DESeq paper: "While one may not want to draw
# strong conclusions from such an analysis, it may still be useful for
# exploration and hypothesis generation."

dds <- DESeqDataSetFromMatrix(countData=countMat, colData=metadataFinal, design=~1)

dds <- DESeq(dds, parallel=TRUE)

## regularized log transformation
rld <- rlog(dds, blind=TRUE, fast=TRUE)

finalmat <- cbind(data.frame(ID=rownames(rld)), as.data.frame(assay(rld)))
write.csv(finalmat, file="DESeq2_rlog.csv", row.names=FALSE)

