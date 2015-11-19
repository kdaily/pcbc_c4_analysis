designMat2CovEffects <- function(COVARIATES_MAP, sigCovars.Effects, type='mean'){
  if (type == 'mean')
    sapply(COVARIATES_MAP,function(x,y = sigCovars.Effects){return(mean(y[names(y) %in% x]))})
  else if (type == 'sum') 
    sapply(COVARIATES_MAP,function(x,y = sigCovars.Effects){return(sum(y[names(y) %in% x]))})
}
