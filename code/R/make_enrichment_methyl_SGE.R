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
DMETHYL.ALL_ID = 'syn4527629'
DMETHYL.ALL = data.table::fread(synGet(DMETHYL.ALL_ID)@filePath, data.table=F)

# Filter comparisons
Comparisons.ALL = split(DMETHYL.ALL, factor(DMETHYL.ALL$Comparison))
Comparisons.ALL = lapply(Comparisons.ALL, function(x){return(unique(x$nearestTx[abs(x$logFC) >= log2(1.3) & x$adj.P.value<=0.05]))})

# Get differential expression (EACH)
DMETHYL.EACH_ID = 'syn4598861'
DMETHYL.EACH = data.table::fread(synGet(DMETHYL.EACH_ID)@filePath, data.table=F)

# Filter comparisons
Comparisons.EACH = split(DMETHYL.EACH, factor(DMETHYL.EACH$Comparison))
Comparisons.EACH = lapply(Comparisons.EACH, function(x){return(unique(x$nearestTx[abs(x$logFC) >= log2(1.2) & x$adj.P.value<=0.05]))})

Comparisons = c(Comparisons.ALL, Comparisons.EACH)

# Make directory and write shell scripts for running these files
system('mkdir sgeEnrichMethyl')
fp_all = file(paste('sgeEnrichMethyl/allSubmissions.sh'),'w+')    
cat('#!/bin/bash',file=fp_all,sep='\n')
close(fp_all)
for (id in names(Comparisons)){
  fp = file (paste('/home/ec2-user/Work/Github/pcbc_c4_analysis/code/R/sgeEnrichMethyl/SUB',id,sep='.'), "w+")
  cat('#!/bin/bash', 
      'sleep 30', 
      paste('Rscript /home/ec2-user/Work/Github/pcbc_c4_analysis/code/R/enrichment_methyl_SGE.R',id), 
      file = fp,
      sep = '\n')
  close(fp)
  
  fp_all = file(paste('sgeEnrichMethyl/allSubmissions.sh'),'a+')    
  cat(paste('qsub','-cwd','-V',paste('/home/ec2-user/Work/Github/pcbc_c4_analysis/code/R/sgeEnrichMethyl/SUB',id,sep='.'),
            '-o',paste('/home/ec2-user/Work/Github/pcbc_c4_analysis/code/R/sgeEnrichMethyl/SUB',id,'o',sep='.'),
            '-e',paste('/home/ec2-user/Work/Github/pcbc_c4_analysis/code/R/sgeEnrichMethyl/SUB',id,'e',sep='.')),
      file=fp_all,
      sep='\n')
  close(fp_all)
}
