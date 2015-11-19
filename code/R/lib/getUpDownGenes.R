getUpDownGenes <- function(PVAL,logFC, Genes, contName,PVAL_CUTOFF = 0.05, FC_CUTOFF = 1e-10){
  ind.pos <- which(PVAL<=PVAL_CUTOFF & logFC >= FC_CUTOFF)
  if (length(ind.pos) > 0){
    tmp.pos <- data.frame(GeneSymbol = Genes[ind.pos],
                          logFC = as.numeric(logFC[ind.pos]),
                          adj.P.value = as.numeric(PVAL[ind.pos]),
                          Comparison = paste0(contName,'__up'))
  } else {
    tmp.pos <- data.frame()
  }
  
  ind.neg <- which(PVAL<=PVAL_CUTOFF & logFC <= -FC_CUTOFF)
  if (length(ind.neg) > 0){
    tmp.neg <- data.frame(GeneSymbol = Genes[ind.neg],
                          logFC = as.numeric(logFC[ind.neg]),
                          adj.P.value = as.numeric(PVAL[ind.neg]),
                          Comparison = paste0(contName,'__down'))  
  } else {
    tmp.neg <- data.frame()
  }
  
  return(rbind(tmp.pos,tmp.neg))
}