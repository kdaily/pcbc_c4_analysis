#!usr/bin/env Rscript

# Submission Script in R
# Clear R console screen output
cat("\014")

# Clear R workspace
setwd('/home/ec2-user/Work/Github/pcbc_c4_analysis/code/R')

# Load libraries
library(synapseClient)

# login to synapse
synapseLogin()

# Get differential methylation (ALL)
DMRNA.ALL_ID = 'syn4484232'
DMRNA.ALL = data.table::fread(synGet(DMRNA.ALL_ID)@filePath, data.table=F)

# Filter comparisons
Comparisons.ALL = split(DMRNA.ALL, factor(DMRNA.ALL$Comparison))
Comparisons.ALL = lapply(Comparisons.ALL, function(x){return(unique(x$nearestTx[abs(x$logFC) >= log2(1.3) & x$adj.P.value<=0.05]))})

# Get differential expression (EACH)
DMRNA.EACH_ID = 'syn4485993'
DMRNA.EACH = data.table::fread(synGet(DMRNA.EACH_ID)@filePath, data.table=F)

# Filter comparisons
Comparisons.EACH = split(DMRNA.EACH, factor(DMRNA.EACH$Comparison))
Comparisons.EACH = lapply(Comparisons.EACH, function(x){return(unique(x$nearestTx[abs(x$logFC) >= log2(1.2) & x$adj.P.value<=0.05]))})

Comparisons = c(Comparisons.ALL, Comparisons.EACH)

# Make directory and write shell scripts for running these files
system('mkdir sgeEnrichPCBC')
fp_all = file(paste('sgeEnrichPCBC/allSubmissions.sh'),'w+')    
cat('#!/bin/bash',file=fp_all,sep='\n')
close(fp_all)
for (id in names(Comparisons)){
  fp = file (paste('/home/ec2-user/Work/Github/pcbc_c4_analysis/code/R/sgeEnrichPCBC/SUB',id,sep='.'), "w+")
  cat('#!/bin/bash', 
      'sleep 30', 
      paste('Rscript /home/ec2-user/Work/Github/pcbc_c4_analysis/code/R/enrichment_mRNA_SGE.R',id), 
      file = fp,
      sep = '\n')
  close(fp)
  
  fp_all = file(paste('sgeEnrichPCBC/allSubmissions.sh'),'a+')    
  cat(paste('qsub','-cwd','-V',paste('/home/ec2-user/Work/Github/pcbc_c4_analysis/code/R/sgeEnrichPCBC/SUB',id,sep='.'),
            '-o',paste('/home/ec2-user/Work/Github/pcbc_c4_analysis/code/R/sgeEnrichPCBC/SUB',id,'o',sep='.'),
            '-e',paste('/home/ec2-user/Work/Github/pcbc_c4_analysis/code/R/sgeEnrichPCBC/SUB',id,'e',sep='.')),
      file=fp_all,
      sep='\n')
  close(fp_all)
}
