getAssociationStatistics <- function(COVARIATES,FactorCovariates,ContCovariates,PVAL=0.05){
  
  getFactorAssociationStatistics <- function(factorNames,COVARIATES, na.action='remove'){
    if (na.action == "remove")
      COVARIATES = na.omit(COVARIATES[,factorNames])
    fac1 = as.factor(COVARIATES[,1])
    fac2 = as.factor(COVARIATES[,2])
    
    stats = assocstats(xtabs(~fac1+fac2))
    
    return(c(Estimate=stats$cramer,Pval=stats$chisq_tests['Pearson','P(> X^2)']))
  }
  
  # Find association between factor covariates
  COVARIATES.FACTOR.CORRELATION = apply(expand.grid(FactorCovariates,FactorCovariates),1,getFactorAssociationStatistics,COVARIATES[,FactorCovariates])
  COVARIATES.FACTOR.CORRELATION.ESTIMATE <- matrix(COVARIATES.FACTOR.CORRELATION['Estimate',],nrow=length(FactorCovariates),ncol=length(FactorCovariates))
  COVARIATES.FACTOR.CORRELATION.PVAL <- matrix(COVARIATES.FACTOR.CORRELATION['Pval',],nrow=length(FactorCovariates),ncol=length(FactorCovariates))
  
  colnames(COVARIATES.FACTOR.CORRELATION.ESTIMATE) <- FactorCovariates
  rownames(COVARIATES.FACTOR.CORRELATION.ESTIMATE) <- FactorCovariates
  
  colnames(COVARIATES.FACTOR.CORRELATION.PVAL) <- FactorCovariates
  rownames(COVARIATES.FACTOR.CORRELATION.PVAL) <- FactorCovariates
    
  # Find association between factor covariates
  tmp <- apply(COVARIATES[,ContCovariates],2,as.numeric)
  COVARIATES.CONT.CORRELATION = corr.test(tmp,use = 'pairwise.complete.obs')
  COVARIATES.CONT.CORRELATION.ESTIMATE = COVARIATES.CONT.CORRELATION$r
  COVARIATES.CONT.CORRELATION.PVAL = COVARIATES.CONT.CORRELATION$p
  
  # Find ICC between factor and continuous covariates
  getFactorContAssociationStatistics <- function(factorContNames,COVARIATES, na.action='remove'){
    if (na.action == "remove")
      COVARIATES = na.omit(COVARIATES[,factorContNames])
    
    stats = ICC(COVARIATES)
    
    return(c(Estimate = stats$results['Single_raters_absolute','ICC'],
             Pval = stats$results['Single_raters_absolute','p']))
  }
  
  # Find association between factor covariates
  COVARIATES.FACTORCONT.CORRELATION = apply(expand.grid(FactorCovariates,ContCovariates),1,getFactorContAssociationStatistics,COVARIATES[,c(FactorCovariates,ContCovariates)])
  COVARIATES.FACTORCONT.CORRELATION.ESTIMATE <- matrix(COVARIATES.FACTORCONT.CORRELATION['Estimate',],nrow=length(FactorCovariates),ncol=length(ContCovariates))
  COVARIATES.FACTORCONT.CORRELATION.PVAL <- matrix(COVARIATES.FACTORCONT.CORRELATION['Pval',],nrow=length(FactorCovariates),ncol=length(ContCovariates))
  
  colnames(COVARIATES.FACTORCONT.CORRELATION.ESTIMATE) <- ContCovariates
  rownames(COVARIATES.FACTORCONT.CORRELATION.ESTIMATE) <- FactorCovariates
  
  colnames(COVARIATES.FACTORCONT.CORRELATION.PVAL) <- ContCovariates
  rownames(COVARIATES.FACTORCONT.CORRELATION.PVAL) <- FactorCovariates
  
  # Combine all estimates that are significant
  COVARIATES.CORRELATION.ESTIMATE = rbind(cbind(COVARIATES.FACTOR.CORRELATION.ESTIMATE,COVARIATES.FACTORCONT.CORRELATION.ESTIMATE),
                                          cbind(t(COVARIATES.FACTORCONT.CORRELATION.ESTIMATE),COVARIATES.CONT.CORRELATION.ESTIMATE))
  
  COVARIATES.CORRELATION.PVAL = rbind(cbind(COVARIATES.FACTOR.CORRELATION.PVAL,COVARIATES.FACTORCONT.CORRELATION.PVAL),
                                          cbind(t(COVARIATES.FACTORCONT.CORRELATION.PVAL),COVARIATES.CONT.CORRELATION.PVAL))
  
  # plot heatmap
  tmp <- COVARIATES.CORRELATION.ESTIMATE
  tmp[COVARIATES.CORRELATION.PVAL>PVAL] <- 0
  p <- ggheatmap(abs(tmp),hm.colours=brewer.pal(9,'Reds'))
  
  return(list(ESTIMATE = COVARIATES.CORRELATION.ESTIMATE, PVAL = COVARIATES.CORRELATION.PVAL, plot=p))  
}