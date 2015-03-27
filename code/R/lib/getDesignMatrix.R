# Function to optain desing matrix (modified from covairates pipeline of Menachem Former)
getDesignMatrix <- function(covariatesDataFrame, FACTOR_COVARIATE_NAMES, RELEVELS=list()) {
  ROWNAMES = rownames(covariatesDataFrame)
  COLNAMES = colnames(covariatesDataFrame)
  FACTOR_COVARIATE_NAMES = setdiff(FACTOR_COVARIATE_NAMES, FACTOR_COVARIATE_NAMES[!(FACTOR_COVARIATE_NAMES %in% colnames(covariatesDataFrame))])
  NUMERIC_COVARIATE_NAMES = setdiff(COLNAMES, FACTOR_COVARIATE_NAMES)
  
  # Ensure the factors are in fact of type factor, and the quantitative variables are numeric:
  covariatesDataFrame = as.data.frame( lapply(colnames(covariatesDataFrame), function(column) {if (column %in% FACTOR_COVARIATE_NAMES) {fac = as.factor(covariatesDataFrame[, column]); if (column %in% names(RELEVELS)) {fac = relevel(fac, ref=RELEVELS[[column]])}; return(fac)} else {return(as.numeric(covariatesDataFrame[, column]))}}) )
  rownames(covariatesDataFrame) = ROWNAMES
  colnames(covariatesDataFrame) = COLNAMES
  
  contra = NULL
  MAX_NUM_CATS = Inf
  catData = covariatesDataFrame[, FACTOR_COVARIATE_NAMES, drop=FALSE]
  if (ncol(catData) > 0) {
    numCats = sapply(colnames(catData), function(col) nlevels(factor(catData[, col])))
    EXCLUDE_CATEGORICAL_COLS = names(numCats)[numCats <= 1 | numCats > MAX_NUM_CATS]
    if (!is.null(EXCLUDE_CATEGORICAL_COLS) && length(EXCLUDE_CATEGORICAL_COLS) > 0) {
      warning(paste("Excluding categorical variables with less than 2", ifelse(is.infinite(MAX_NUM_CATS), "", paste(" or more than ", MAX_NUM_CATS, sep="")), " categories: ", paste(paste("'", EXCLUDE_CATEGORICAL_COLS, "'", sep=""), collapse=", "), sep=""))
      FACTOR_COVARIATE_NAMES = setdiff(FACTOR_COVARIATE_NAMES, EXCLUDE_CATEGORICAL_COLS)
      covariatesDataFrame = covariatesDataFrame[, !(colnames(covariatesDataFrame) %in% EXCLUDE_CATEGORICAL_COLS), drop=FALSE]
    }
    
    # Inspired by http://stackoverflow.com/questions/4560459/all-levels-of-a-factor-in-a-model-matrix-in-r
    #
    # And, already ensured above that covariatesDataFrame[, FACTOR_COVARIATE_NAMES] satisfies:
    # 1) fac is of type factor.
    # 2) fac is releveled as designated in RELEVELS.
    contra = lapply(FACTOR_COVARIATE_NAMES, function(column) {fac = covariatesDataFrame[, column]; fac = contrasts(fac);})
    names(contra) = FACTOR_COVARIATE_NAMES
  }
  
  # Inspired by http://stackoverflow.com/questions/5616210/model-matrix-with-na-action-null :
  current.na.action = getOption('na.action')
  # Model matrix will now include "NA":
  options(na.action='na.pass')
  
  # NOTE: this includes an '(Intercept)' column:
  design = model.matrix(~ ., data=covariatesDataFrame, contrasts.arg=contra)
  rownames(design) = rownames(covariatesDataFrame)
  
  options(na.action=current.na.action)
  
  return(list(design=design, covariates=COLNAMES, factorsLevels=sapply(contra, colnames, simplify=FALSE), numericCovars=NUMERIC_COVARIATE_NAMES, covariatesDataFrame=covariatesDataFrame))
}
