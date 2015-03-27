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