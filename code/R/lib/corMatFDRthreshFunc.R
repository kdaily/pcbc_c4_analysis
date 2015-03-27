# Function for FDR thresholding of correlation matrix
corMatFDRthreshFunc <- function(cor_mat, indicesMask=NULL, MAX_FDR = 0.1) {    
  if (is.null(indicesMask)) {indicesMask = 1:nrow(cor_mat)}
  
  fdr = rep(1.0, nrow(cor_mat))
  fdr[indicesMask] = p.adjust(cor_mat$pvalue[indicesMask], method="fdr")
  
  return (fdr <= MAX_FDR)
}
