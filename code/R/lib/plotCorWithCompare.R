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
  gRes = gRes + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=14, face=use.face.x, color=use.color.x), axis.text.y=element_text(size=14, color="black", hjust=0))
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
