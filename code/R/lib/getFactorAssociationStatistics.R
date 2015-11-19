getFactorAssociationStatistics <- function(factorNames,COVARIATES, na.action='remove'){
  
  require(vcd)
  
  if (na.action == "remove")
    COVARIATES = na.omit(COVARIATES[,factorNames])
  
  fac1 = as.factor(COVARIATES[,1])
  fac2 = as.factor(COVARIATES[,2])
  
  stats = assocstats(xtabs(~fac1+fac2))
  
  return(c(Estimate=stats$cramer,Pval=stats$chisq_tests['Pearson','P(> X^2)']))
}
