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

# Get differential miRNA (ALL)
DMIRNA.ALL_ID = 'syn4609631'
DMIRNA.ALL = data.table::fread(synGet(DMIRNA.ALL_ID)@filePath, data.table=F)

# Filter comparisons
Comparisons.ALL = split(DMIRNA.ALL, factor(DMIRNA.ALL$Comparison))

# Get differential expression (EACH)
DMIRNA.EACH_ID = 'syn4622430'
DMIRNA.EACH = data.table::fread(synGet(DMIRNA.EACH_ID)@filePath, data.table=F)

# Filter comparisons
Comparisons.EACH = split(DMIRNA.EACH, factor(DMIRNA.EACH$Comparison))

Comparisons = c(Comparisons.ALL, Comparisons.EACH)

# Make directory and write shell scripts for running these files
system('mkdir sgeEnrichMirna')
fp_all = file(paste('sgeEnrichMirna/allSubmissions.sh'),'w+')    
cat('#!/bin/bash',file=fp_all,sep='\n')
close(fp_all)
for (id in names(Comparisons)){
  fp = file (paste('/home/ec2-user/Work/Github/pcbc_c4_analysis/code/R/sgeEnrichMirna/SUB',id,sep='.'), "w+")
  cat('#!/bin/bash', 
      'sleep 30', 
      paste('Rscript /home/ec2-user/Work/Github/pcbc_c4_analysis/code/R/enrichment_miRNA_SGE.R',id), 
      file = fp,
      sep = '\n')
  close(fp)
  
  fp_all = file(paste('sgeEnrichMirna/allSubmissions.sh'),'a+')    
  cat(paste('qsub','-cwd','-V',paste('/home/ec2-user/Work/Github/pcbc_c4_analysis/code/R/sgeEnrichMirna/SUB',id,sep='.'),
            '-o',paste('/home/ec2-user/Work/Github/pcbc_c4_analysis/code/R/sgeEnrichMirna/SUB',id,'o',sep='.'),
            '-e',paste('/home/ec2-user/Work/Github/pcbc_c4_analysis/code/R/sgeEnrichMirna/SUB',id,'e',sep='.')),
      file=fp_all,
      sep='\n')
  close(fp_all)
}
