# Find Inter Class Correlation between factor and continuous covariates
# Inspired from http://stats.stackexchange.com/questions/108007/correlations-with-categorical-variables
getFactorContAssociationStatistics <- function(factorContNames,COVARIATES, na.action='remove', 
                                               alpha = 0.05){
  
  require(psych)
  
  if (na.action == "remove")
    COVARIATES = na.omit(COVARIATES[,factorContNames])
  
  stats = ICC(COVARIATES[,factorContNames], alpha = alpha)
  
  Pval = summary(aov(COVARIATES[,factorContNames[1]]~COVARIATES[,factorContNames[2]]))[[1]][["Pr(>F)"]][1]
  
  
  return(c(Estimate = stats$results['Single_raters_absolute','ICC'],
           Pval = Pval))
}
