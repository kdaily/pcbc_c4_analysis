# Function to run principal component analysis and plot correlations
runPCAandPlotCorrelations <- function(genesBySamples, samplesByCovariates, dataName, isKeyPlot=FALSE, 
                                      SCALE_DATA_FOR_PCA = TRUE, MIN_PVE_PCT_PC = 1.0, CORRELATION_TYPE = "pearson",
                                      ALSO_PLOT_ALL_COVARS_VS_PCA = TRUE, MAX_NUM_LEVELS_PER_COVAR = 50) {
      
  title = paste(ifelse(SCALE_DATA_FOR_PCA, "S", "Un-s"), "caled ", dataName, " ", " data in PCA; PVE >= ", MIN_PVE_PCT_PC, "%; ", CORRELATION_TYPE, " correlations ", sep="")
  writeLines(paste("\nRunning PCA and calculating correlations for:\n", title, sep=""))
 
  pcaRes <- runPCA(genesBySamples=genesBySamples, 
                   SCALE_DATA_FOR_PCA=SCALE_DATA_FOR_PCA, 
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
    corrRes = calcCompleteCorAndPlot(samplePCvals, 
                                     samplesByCovariates, 
                                     CORRELATION_TYPE, 
                                     title, 
                                     WEIGHTS = pve[1:dim(samplePCvals)[2]],
                                     PLOT_ALL_COVARS, 
                                     EXCLUDE_VARS_FROM_FDR)
    add_PC_res[[length(add_PC_res)+1]] = list(plotData=corrRes$plot, isKeyPlot=(isKeyPlot && !PLOT_ALL_COVARS))
    if (!PLOT_ALL_COVARS) {
      significantCovars = corrRes$significantCovars
      Effects.significantCovars = corrRes$Effects.significantCovars
    }
  }
  
  return(list(significantCovars=significantCovars, PC_res=add_PC_res, Effects.significantCovars = Effects.significantCovars))
}
