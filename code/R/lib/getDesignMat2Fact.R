getDesignMat2Fact <- function(designMat, FactorCovariates){
  ind <- grep(paste(FactorCovariates,collapse='|'),colnames(designMat))
  tmp <- as.data.frame(designMat)
  
  tmp[,ind] <- lapply(tmp[,ind],function(cols,FactorCovariates){cols <- factor(as.character(cols))})
  return(tmp)                      
}