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
