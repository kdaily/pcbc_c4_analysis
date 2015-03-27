# Function to optain residual matrix (modified from covairates pipeline of Menachem Former)
calcResiduals <- function(geneBySampleValues, samplesByCovariates, factorCovariates = NULL, varsToAddBackIn=NULL, sampleWeights=NULL) {  
  
  # Convert factor covariates data frame to numeric matrix
  if (!is.null(factorCovariates)){
    colNames = intersect(colnames(samplesByCovariates),factorCovariates)
    samplesByCovariates[,colNames] = apply(samplesByCovariates[,colNames,drop=F],2, function(cols){unclass(factor(cols))})
  }    
  samplesByCovariates = as.matrix(samplesByCovariates)
  
  ##############################################################################
  #### If sampleWeights are given as matrix use calcResiduals in a for loop ####
  ##############################################################################
  if (is.matrix(sampleWeights)) {
    residualizedMat = matrix(NA, nrow=nrow(geneBySampleValues), ncol=ncol(geneBySampleValues), dimnames=dimnames(geneBySampleValues))
    for (gInd in 1:nrow(geneBySampleValues)) {
      gRow = calcResiduals(geneBySampleValues[gInd, , drop=FALSE], samplesByCovariates, varsToAddBackIn, sampleWeights[gInd, ])
      residualizedMat[gInd, ] = gRow
    }
    return(residualizedMat)
  }
  #################################################################################
  
  # If lest square model is needed (uncomment the following line)
  # result.lm = lsfit(x=samplesByCovariates, y=t(geneBySampleValues), wt=sampleWeights, intercept=FALSE)
  
  # Formula of "y ~ 0 + x" means no intercept:
  result.lm = lm(t(geneBySampleValues) ~ 0 + samplesByCovariates, weights=sampleWeights)
  covarNames = colnames(samplesByCovariates)
  
  coef = result.lm$coefficients
  isMatrixForm = is.matrix(coef)
  if (isMatrixForm) {
    rownames(coef) = covarNames
  }
  else {
    names(coef) = covarNames
  }
  
  allVarsToAddBack = '(Intercept)'
  if (!is.null(varsToAddBackIn)) {
    allVarsToAddBack = c(allVarsToAddBack, varsToAddBackIn)
  }
  allVarsToAddBack = intersect(allVarsToAddBack, covarNames)
  
  residualizedMat = result.lm$residuals
  for (v in allVarsToAddBack) {
    if (isMatrixForm) {
      multCoef = coef[v, , drop=FALSE]
    }
    else {
      multCoef = coef[v]
    }
    residualizedMat = residualizedMat + samplesByCovariates[, v, drop=FALSE] %*% multCoef
  }
  
  residualizedMat = t(residualizedMat)
  rownames(residualizedMat) = rownames(geneBySampleValues)
  colnames(residualizedMat) = colnames(geneBySampleValues)
  
  return(residualizedMat)
}
