get450KProbeMapping <- function(probeIDs, genome='hg19'){
 require(FDb.InfiniumMethylation.hg19)
  require(dplyr)
  require(plyr)
  
  hm450 <- getPlatform(platform = 'HM450', genome = genome)
  probes <- hm450[probeIDs]
    
  TSS = getNearestTSS(probes)
  TSS = rownameToFirstColumn(TSS,'methProbeIDs')
  TSS = dplyr::select(TSS, methProbeIDs, distance, nearestGeneSymbol, nearestTranscript)
  setnames(TSS,c("distance", "nearestGeneSymbol", "nearestTranscript"),
           c("distanceToTSS", "nearestTSS", "nearestTSS.ID"))
  
  Tx = getNearestTranscript(probes)
  Tx = rownameToFirstColumn(Tx, 'methProbeIDs')
  Tx = dplyr::select(Tx, methProbeIDs, distance, nearestGeneSymbol, nearestTranscript)
  setnames(Tx,c("distance", "nearestGeneSymbol", "nearestTranscript"),
           c("distanceToTx", "nearestTx", "nearestTx.ID"))
  
  hm450 = rownameToFirstColumn(hm450, 'methProbeIDs')
  
  Annotation = join_all(list(hm450,TSS,Tx), by = 'methProbeIDs', match = 'all')
  
  return(list(Annotation = Annotation))
}