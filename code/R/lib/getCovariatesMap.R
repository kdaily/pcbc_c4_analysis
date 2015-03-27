getCovariatesMap <- function(DESIGN_AND_COVARS) {
  FACTORS_AND_LEVELS = DESIGN_AND_COVARS$factorsLevels
  COVARIATES_MAP = sapply(1:length(FACTORS_AND_LEVELS), function(i) {paste(deparse(as.name(names(FACTORS_AND_LEVELS)[i]), backtick=TRUE), FACTORS_AND_LEVELS[[i]], sep="")}, simplify=FALSE)
  names(COVARIATES_MAP) = names(FACTORS_AND_LEVELS)
  
  numericCovars = DESIGN_AND_COVARS$numericCovars
  names(numericCovars) = numericCovars
  numericCovars = getSymbolicNamesList(numericCovars)
  COVARIATES_MAP = c(COVARIATES_MAP, numericCovars)
  
  # Re-order as originally given:
  COVARIATES_MAP = COVARIATES_MAP[DESIGN_AND_COVARS$covariates]
  
  return(COVARIATES_MAP)
}

