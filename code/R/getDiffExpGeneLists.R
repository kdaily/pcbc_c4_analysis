library(dplyr)
library(synapseClient)
library(data.table)

synapseLogin()

FC = c(1,1.5,2)
PVAL = c(0.01,0.05,0.1)

for (fc in FC){
  for (pval in PVAL){
    # Get differential expression (ALL)
    DEXP.ALL_ID = 'syn4484232'
    ALL_USED_IDs = DEXP.ALL_ID
    DEXP.ALL = data.table::fread(synGet(DEXP.ALL_ID)@filePath, data.table=F)
    
    # Filter comparisons
    Comparisons.ALL = split(DEXP.ALL, factor(DEXP.ALL$Comparison))
    Comparisons.ALL = lapply(Comparisons.ALL, function(x){return(unique(x$GeneSymbol[abs(x$logFC) >= log2(fc) & x$adj.P.value<=pval]))})
    
    # Get differential expression (EACH)
    DEXP.EACH_ID = 'syn4485993'
    ALL_USED_IDs = c(ALL_USED_IDs, DEXP.EACH_ID)
    DEXP.EACH = data.table::fread(synGet(DEXP.EACH_ID)@filePath, data.table=F)
    
    # Filter comparisons
    Comparisons.EACH = split(DEXP.EACH, factor(DEXP.EACH$Comparison))
    Comparisons.EACH = lapply(Comparisons.EACH, function(x){return(unique(x$GeneSymbol[abs(x$logFC) >= log2(fc) & x$adj.P.value<=pval]))})
    
    Comparisons = c(Comparisons.ALL, Comparisons.EACH)
    
    # Write to table
    write.table(paste(c('ComparisonName','GeneList'),collapse='\t'), file=paste('mRNA',fc,pval,'tsv',sep='.'),sep='\t',quote=F,col.names=F,row.names=F)
    for (name in names(Comparisons)){
      write.table(paste(c(name,paste(Comparisons[[name]],collapse=',')),collapse='\t'), 
                  file=paste('mRNA',fc,pval,'tsv',sep='.'),
                  sep='\t', quote=F,col.names=F,row.names=F,append=T)
    }
    
    # Write to synapse
    OBJ <- File(paste('mRNA',fc,pval,'tsv',sep='.'),
                name = paste('mRNA','FC',fc,'PVAL',pval,'tsv',sep=' '),
                parentId = 'syn4640410')
    OBJ <- synStore(OBJ, 
                    used = ALL_USED_IDs, 
                    activityName = 'Subsetting differential expression genelists')    
  }
}

# methylation
for (fc in FC){
  for (pval in PVAL){
    # Get differential methylation (ALL)
    DMETHYL.ALL_ID = 'syn4527629'
    ALL_USED_IDs = DMETHYL.ALL_ID
    DMETHYL.ALL = data.table::fread(synGet(DMETHYL.ALL_ID)@filePath, data.table=F)
    
    # Filter comparisons
    Comparisons.ALL = split(DMETHYL.ALL, factor(DMETHYL.ALL$Comparison))
    Comparisons.ALL = lapply(Comparisons.ALL, function(x){return(unique(x$methProbeIDs[abs(x$logFC) >= log2(fc) & x$adj.P.value<=pval]))})
    
    # Get differential expression (EACH)
    DMETHYL.EACH_ID = 'syn4598861'
    ALL_USED_IDs = c(ALL_USED_IDs, DMETHYL.EACH_ID)
    DMETHYL.EACH = data.table::fread(synGet(DMETHYL.EACH_ID)@filePath, data.table=F)
    
    # Filter comparisons
    Comparisons.EACH = split(DMETHYL.EACH, factor(DMETHYL.EACH$Comparison))
    Comparisons.EACH = lapply(Comparisons.EACH, function(x){return(unique(x$methProbeIDs[abs(x$logFC) >= log2(fc) & x$adj.P.value<=pval]))})
    
    Comparisons = c(Comparisons.ALL, Comparisons.EACH)
    
    # Write to table
    write.table(paste(c('ComparisonName','methProbeIDs'),collapse='\t'), file=paste('methylation',fc,pval,'tsv',sep='.'),sep='\t',quote=F,col.names=F,row.names=F)
    for (name in names(Comparisons)){
      write.table(paste(c(name,paste(Comparisons[[name]],collapse=',')),collapse='\t'), 
                  file=paste('methylation',fc,pval,'tsv',sep='.'),
                  sep='\t', quote=F,col.names=F,row.names=F,append=T)
    }
    
    # Write to synapse
    OBJ <- File(paste('methylation',fc,pval,'tsv',sep='.'),
                name = paste('methylation','FC',fc,'PVAL',pval,'tsv',sep=' '),
                parentId = 'syn4640410')
    OBJ <- synStore(OBJ, 
                    used = ALL_USED_IDs, 
                    activityName = 'Subsetting differential expression probes list')    
  }
}

# miRNA
for (fc in FC){
  for (pval in PVAL){
    # Get differential miRNA (ALL)
    DMIRNA.ALL_ID = 'syn4609631'
    ALL_USED_IDs = DMIRNA.ALL_ID
    DMIRNA.ALL = data.table::fread(synGet(DMIRNA.ALL_ID)@filePath, data.table=F)
    
    # Filter comparisons
    Comparisons.ALL = split(DMIRNA.ALL, factor(DMIRNA.ALL$Comparison))
    Comparisons.ALL = lapply(Comparisons.ALL, function(x){return(unique(x$GeneSymbol[abs(x$logFC) >= log2(fc) & x$adj.P.value<=pval]))})
    
    # Get differential miRNA (EACH)
    DMIRNA.EACH_ID = 'syn4622430'
    ALL_USED_IDs = c(ALL_USED_IDs, DMIRNA.EACH_ID)
    DMIRNA.EACH = data.table::fread(synGet(DMIRNA.EACH_ID)@filePath, data.table=F)
    
    # Filter comparisons
    Comparisons.EACH = split(DMIRNA.EACH, factor(DMIRNA.EACH$Comparison))
    Comparisons.EACH = lapply(Comparisons.EACH, function(x){return(unique(x$GeneSymbol[abs(x$logFC) >= log2(fc) & x$adj.P.value<=pval]))})
    
    Comparisons = c(Comparisons.ALL, Comparisons.EACH)
    
    # Write to table
    write.table(paste(c('ComparisonName','GeneLists'),collapse='\t'), file=paste('miRNA',fc,pval,'tsv',sep='.'),sep='\t',quote=F,col.names=F,row.names=F)
    for (name in names(Comparisons)){
      write.table(paste(c(name,paste(Comparisons[[name]],collapse=',')),collapse='\t'), 
                  file=paste('miRNA',fc,pval,'tsv',sep='.'),
                  sep='\t', quote=F,col.names=F,row.names=F,append=T)
    }
    
    # Write to synapse
    OBJ <- File(paste('miRNA',fc,pval,'tsv',sep='.'),
                name = paste('miRNA','FC',fc,'PVAL',pval,'tsv',sep=' '),
                parentId = 'syn4640410')
    OBJ <- synStore(OBJ, 
                    used = ALL_USED_IDs, 
                    activityName = 'Subsetting differential expression miRNA list')    
  }
}
