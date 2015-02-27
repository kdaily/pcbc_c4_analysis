# Utility Functions

# Function to find PC and significant covariates
findPCandSigCovariates <- function(GenesBySamples,SamplesByCovariates,minPCsdev){
  # Find PCA 
  pca.res = prcomp(t(GenesBySamples),center=T,scale.=T,retx=T)
  
  # Consider PCs with sdev >= minPCsdev
  npca = max(1,which(pca.res$sdev>=minPCsdev))
  
  SamplesByPCvals = pca.res$x[, 1:npca, drop=FALSE]
  PCsdev = pca.res$sdev[1:npca]
  colnames(SamplesByPCvals) = paste(colnames(SamplesByPCvals), " (", sprintf("%.2f", PCsdev[1:npca]), "%)", sep="")
  
  # Find correlation b/w PCs and Covariates
  CorrEstimate = matrix(0,dim(SamplesByPCvals)[2],dim(SamplesByCovariates)[2])
  CorrPval = matrix(0,dim(SamplesByPCvals)[2],dim(SamplesByCovariates)[2])
  for (i in 1:dim(SamplesByPCvals)[2])
    for (j in 1:dim(SamplesByCovariates)[2]){
      tmp = cor.test(SamplesByPCvals[,i], unclass(factor(SamplesByCovariates[,j])), 
                     use='pairwise.complete.obs', method='pearson', adjust="none")
      CorrPval[i,j] = tmp$p.value
      CorrEstimate[i,j] = tmp$estimate
    }
  return(list(SamplesByPCvals=SamplesByPCvals,PCsdev=PCsdev,
              CorrPval=CorrPval,CorrEstimate=CorrEstimate))
}

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

# Function to convert counts to cpm and filter cpm counts matrix
getGeneFilteredGeneExprMatrix <- function(genesBySamplesCounts,ONLY_USE_GENES=NULL,
                                          MIN_GENE_CPM=1,
                                          MIN_SAMPLE_PERCENT_WITH_MIN_GENE_CPM=0.5,
                                          EDGER_NORMALIZATION.calcNormFactors.method = "none", # Choices are: "TMM", "RLE", "upperquartile", "none" [see edgeR::calcNormFactors()]
                                          EDGER_NORMALIZATION.keep.lib.sizes = TRUE) {
  if (!is.null(ONLY_USE_GENES)) {
    useGenes = colnames(genesBySamplesCounts)
    useGenes = useGenes[useGenes %in% ONLY_USE_GENES]
    genesBySamplesCounts = genesBySamplesCounts[, useGenes]
    writeLines(paste("\nLimiting expression data to ", length(useGenes), " genes specified by the ONLY_USE_GENES parameter.", sep=""))
  }
  
  # Make edgeR object
  MATRIX.ALL_GENES = DGEList(counts=genesBySamplesCounts, genes=rownames(genesBySamplesCounts))
  
  # Keep genes with at least MIN_GENE_CPM count-per-million reads (cpm) in at least (MIN_SAMPLE_PERCENT_WITH_MIN_GENE_CPM)% of the samples:
  MATRIX.ALL_GENES.CPM = cpm(MATRIX.ALL_GENES)
  MATRIX.ALL_GENES.CPM[is.nan(MATRIX.ALL_GENES.CPM)] = 0
  fracSamplesWithMinCPM = rowMeans(MATRIX.ALL_GENES.CPM >= MIN_GENE_CPM)
  isNonLowExpr = fracSamplesWithMinCPM >= MIN_SAMPLE_PERCENT_WITH_MIN_GENE_CPM
  
  MATRIX.NON_LOW_GENES = MATRIX.ALL_GENES[isNonLowExpr, ,keep.lib.sizes=EDGER_NORMALIZATION.keep.lib.sizes]
  MATRIX.NON_LOW_GENES = calcNormFactors(MATRIX.NON_LOW_GENES, method=EDGER_NORMALIZATION.calcNormFactors.method)
  
  writeLines(paste("\nWill normalize expression counts for ", nrow(MATRIX.NON_LOW_GENES), " genes (those with a minimum of ", MIN_GENE_CPM, " CPM in at least ", sprintf("%.2f", 100 * MIN_SAMPLE_PERCENT_WITH_MIN_GENE_CPM), "% of the ", ncol(MATRIX.NON_LOW_GENES), " samples).", sep=""))
  
  FRACTION_BIN_WIDTH = 0.02
  plotFracSamplesWithMinCPM = data.frame(GeneFeature=names(fracSamplesWithMinCPM), fracSamplesWithMinCPM=as.numeric(fracSamplesWithMinCPM))
  gRes = ggplot(plotFracSamplesWithMinCPM, aes(x=fracSamplesWithMinCPM))
  gRes = gRes + geom_vline(xintercept=MIN_SAMPLE_PERCENT_WITH_MIN_GENE_CPM, linetype="solid", col="red")
  gRes = gRes + geom_histogram(color="black", fill="white", binwidth=FRACTION_BIN_WIDTH) #+ scale_x_log10()
  gRes = gRes + xlab(paste("Fraction of samples with at least ", MIN_GENE_CPM, " CPM", sep="")) + ylab("# of genes")
  
  return(list(filteredExprMatrix=MATRIX.NON_LOW_GENES, plotHist=gRes))
}

# Function for FDR thresholding of correlation matrix
corMatFDRthreshFunc <- function(cor_mat, indicesMask=NULL, MAX_FDR = 0.1) {    
  if (is.null(indicesMask)) {indicesMask = 1:nrow(cor_mat)}
  
  fdr = rep(1.0, nrow(cor_mat))
  fdr[indicesMask] = p.adjust(cor_mat$pvalue[indicesMask], method="fdr")
  
  return (fdr <= MAX_FDR)
}

# Function to plot correlation matrix
plotCorWithCompare <- function(plotCor, title=NULL, MARK_CORRELATIONS_NAME=NULL, markColumnsAsMissing=NULL) {
  # Mark the X-axis labels based on user-defined pattern:
  markMissingEntries = rep(FALSE, nrow(plotCor))
  
  use.mark.x = FALSE
  if (!is.null(markColumnsAsMissing)) {
    COVAR_NAMES = as.character(unique(levels(plotCor$COVAR)))
    mark.x = rep(FALSE, length(COVAR_NAMES))
    names(mark.x) = COVAR_NAMES
    mark.x[markColumnsAsMissing] = TRUE
    
    plot.x.labels = levels(plotCor$COVAR)[levels(plotCor$COVAR) %in% plotCor$COVAR]
    use.mark.x = as.character(mark.x[plot.x.labels])
    use.mark.x[is.na(use.mark.x)] = FALSE
    
    markMissingEntries[plotCor$COVAR %in% markColumnsAsMissing] = TRUE
  }
  plotCor$markMissingEntries = markMissingEntries
  
  use.face.x = ifelse(use.mark.x, "bold.italic", "plain")
  use.color.x = ifelse(use.mark.x, "darkgray", "black")
  
  plotSingle = FALSE
  plot_aes = aes(COVAR, COMPARE, fill=r, alpha=as.factor(markMissingEntries))
  if (length(unique(plotCor$r)) <= 1) {
    plot_aes = aes(COVAR, COMPARE, fill=factor(r))
    plotSingle = TRUE
  }
  
  # Reverse the Y-axis:
  plotCor$COMPARE = factor(plotCor$COMPARE, levels=rev(levels(plotCor$COMPARE)))
  
  alphaVals = c('TRUE'=0.85, 'FALSE'=1)
  gRes = ggplot(plotCor, plot_aes) + geom_tile() + scale_alpha_manual(values=alphaVals, guide="none")
  
  gRes = gRes + xlab("") + ylab("")
  gRes = gRes + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=6, face=use.face.x, color=use.color.x), axis.text.y=element_text(size=8, color="black", hjust=0))
  gRes = gRes + theme(panel.grid.major.x=element_line(color="black", linetype="dashed")) + theme(panel.grid.minor.x=element_blank())
  gRes = gRes + theme(panel.grid.major.y=element_blank()) + theme(panel.grid.minor.y=element_blank())
  gRes = gRes + theme(panel.background=element_rect(fill="white"))
  
  if (!plotSingle) {
    gRes = gRes + scale_fill_gradient2(low="blue", high="red")
  }
  
  if (!is.null(plotCor$markSignificantCorrelations) || !is.null(plotCor$markPotentialSignificantCorrelations)) {
    markSizes = c()
    markShapes = c()
    
    if (!is.null(plotCor$markSignificantCorrelations)) {
      sigName = MARK_CORRELATIONS_NAME
      
      useAES = modifyList( aes(size=ifelse(markSignificantCorrelations, "markSig", "noMarkSig")), aes_string(shape=paste("as.factor('", sigName, "')", sep="")) )
      gRes = gRes + geom_point(useAES, na.rm=TRUE)
      markSizes["markSig"] = 6
      markSizes["noMarkSig"] = NA
      
      markShapes[sigName] = utf8ToInt('*')
    }
    if (!is.null(plotCor$markPotentialSignificantCorrelations)) {
      potSigName = paste(MARK_CORRELATIONS_NAME, " (incomplete)", sep="")
      
      useAES = modifyList( aes(size=ifelse(markPotentialSignificantCorrelations, "markPotSig", "noMarkPotSig")), aes_string(shape=paste("as.factor('", potSigName, "')", sep="")) )
      gRes = gRes + geom_point(useAES, na.rm=TRUE)
      markSizes["markPotSig"] = 3
      markSizes["noMarkPotSig"] = NA
      
      markShapes[potSigName] = 21 # open circles
    }
    gRes = gRes + scale_size_manual(values=markSizes, guide="none") + scale_shape_manual(values=markShapes, guide="legend", name='')
    gRes = gRes + guides(shape=guide_legend(override.aes=list(size=6)))
  }
  
  if (!is.null(title)) {
    gRes = gRes + labs(title=title)
  }
  
  return(gRes)
}

# Function to calculate correlation and plot
calcCompleteCorAndPlot <- function(COMPARE_data, COVAR_data, correlationType, title, 
                                   PLOT_ALL_COVARS=FALSE, EXCLUDE_VARS_FROM_FDR=NULL, MAX_FDR = 0.1) {
  all_cor = corr.test(COMPARE_data, COVAR_data, use='pairwise.complete.obs', method=correlationType, adjust="none")
  all_cor_vals = all_cor$r
  all_cor_p = all_cor$p
  
  cor_mat = melt(all_cor_p, varnames=c("COMPARE", "COVAR"))
  colnames(cor_mat)[colnames(cor_mat) == "value"] = "pvalue"
  
  cor_mat$COMPARE = factor(cor_mat$COMPARE, levels=rownames(all_cor_p))
  cor_mat$COVAR = factor(cor_mat$COVAR, levels=colnames(all_cor_p))
  
  cor_mat$r = melt(all_cor_vals)$value
  
  calcFDRrows = rep(TRUE, nrow(cor_mat))
  markColumnsAsMissing = NULL
  if (!is.null(EXCLUDE_VARS_FROM_FDR)) {
    calcFDRrows = !(cor_mat$COVAR %in% EXCLUDE_VARS_FROM_FDR)
    markColumnsAsMissing = intersect(colnames(COVAR_data), EXCLUDE_VARS_FROM_FDR)
  }
  
  # Entries that pass the threshold of "significance":  
  markSignificantCorrelations = corMatFDRthreshFunc(cor_mat, indicesMask=calcFDRrows, MAX_FDR = 0.1)
  significantCorrelatedCovars = sort(unique(cor_mat$COVAR[markSignificantCorrelations]))
  
  markPotentialSignificantCorrelations = corMatFDRthreshFunc(cor_mat)
  # Specially mark only those incomplete covariates that would be significant in the context of all covariates:
  markPotentialSignificantCorrelations = markPotentialSignificantCorrelations & !calcFDRrows
  
  plotRows = 1:nrow(cor_mat)
  if (!PLOT_ALL_COVARS) {
    # Plot all correlations for:
    # 1) Covariates with at least one significant correlation
    # 2) Excluded covariates
    plotRows = (cor_mat$COVAR %in% significantCorrelatedCovars) | !calcFDRrows
  }
  plotCor = na.omit(cor_mat[plotRows, ])
  
  for (markCor in c("markSignificantCorrelations", "markPotentialSignificantCorrelations")) {
    useMarkCor = get(markCor)[plotRows]
    if (markCor != "markPotentialSignificantCorrelations" || length(which(useMarkCor)) > 0) {
      plotCor[, markCor] = useMarkCor[ setdiff(1:length(useMarkCor), as.numeric(attr(plotCor, "na.action"))) ]
    }
  }
  
  plot = plotCorWithCompare(plotCor, title, paste("FDR <= ", MAX_FDR, sep=""), markColumnsAsMissing)
  
  return(list(plot=plot, significantCovars=as.character(significantCorrelatedCovars)))
}

# Function to run principal component analysis
runPCA <- function(genesBySamples, SCALE_DATA_FOR_PCA = TRUE, MIN_PVE_PCT_PC = 1.0) {
  
  # estimate variance in data by PC:
  pca.res <- prcomp(t(genesBySamples), center=TRUE, scale.=SCALE_DATA_FOR_PCA, retx=TRUE)
  
  # examine how much variance is explained by PCs, and consider those with PVE >= (MIN_PVE_PCT_PC %):
  pc.var <- pca.res$sdev^2
  pve <- 100 * (pc.var / sum(pc.var))  
  npca <- max(1,length(which(pve >= MIN_PVE_PCT_PC)))
  
  samplePCvals <- pca.res$x[, 1:npca, drop=FALSE]
  
  list(samplePCvals=samplePCvals, pve=pve)
}

# Function to run principal component analysis and plot correlations
runPCAandPlotCorrelations <- function(genesBySamples, samplesByCovariates, dataName, isKeyPlot=FALSE, 
                                      SCALE_DATA_FOR_PCA = TRUE, MIN_PVE_PCT_PC = 1.0, CORRELATION_TYPE = "pearson",
                                      ALSO_PLOT_ALL_COVARS_VS_PCA = TRUE, MAX_NUM_LEVELS_PER_COVAR = 10) {
  title = paste(ifelse(SCALE_DATA_FOR_PCA, "S", "Un-s"), "caled ", dataName, " ", " data in PCA; PVE >= ", MIN_PVE_PCT_PC, "%; ", CORRELATION_TYPE, " correlations ", sep="")
  writeLines(paste("\nRunning PCA and calculating correlations for:\n", title, sep=""))
  
  pcaRes <- runPCA(genesBySamples=genesBySamples, SCALE_DATA_FOR_PCA=SCALE_DATA_FOR_PCA, 
                   MIN_PVE_PCT_PC=MIN_PVE_PCT_PC)
  
  samplePCvals <- pcaRes$samplePCvals
  pve <- pcaRes$pve
  
  npca <- ncol(samplePCvals)
  
  colnames(samplePCvals) = paste(colnames(samplePCvals), " (", sprintf("%.2f", pve[1:npca]), "%)", sep="")
  
  # Find covariates without any missing data
  samplesByFullCovariates = samplesByCovariates[, which(apply(samplesByCovariates, 2, 
                                                              function(dat) all(!is.na(dat))))]
  EXCLUDE_VARS_FROM_FDR = setdiff(colnames(samplesByCovariates), colnames(samplesByFullCovariates))
  
  add_PC_res = list()
  significantCovars = c()
  
  LOOP_PLOT_ALL_COVARS = FALSE
  if (ALSO_PLOT_ALL_COVARS_VS_PCA) { LOOP_PLOT_ALL_COVARS = unique(c(LOOP_PLOT_ALL_COVARS, TRUE)) }
  
  for (PLOT_ALL_COVARS in LOOP_PLOT_ALL_COVARS) {
    corrRes = calcCompleteCorAndPlot(samplePCvals, samplesByCovariates, CORRELATION_TYPE, title, PLOT_ALL_COVARS, EXCLUDE_VARS_FROM_FDR)
    add_PC_res[[length(add_PC_res)+1]] = list(plotData=corrRes$plot, isKeyPlot=(isKeyPlot && !PLOT_ALL_COVARS))
    if (!PLOT_ALL_COVARS) {
      significantCovars = corrRes$significantCovars
    }
  }
  
  return(list(significantCovars=significantCovars, PC_res=add_PC_res))
}
