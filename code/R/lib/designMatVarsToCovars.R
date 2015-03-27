designMatVarsToCovars <- function(COVARIATES_MAP, designMatVars) {
  names(which(sapply(COVARIATES_MAP, function(v) any(v %in% designMatVars))))
}