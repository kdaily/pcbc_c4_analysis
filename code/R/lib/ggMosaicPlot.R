# Function to mosaic plot two factor variables
ggMosaicPlot <- function(var1, var2, var1Name = 'var1', var2Name = 'var2'){
  library(ggplot2)
  levVar1 <- length(levels(var1))
  levVar2 <- length(levels(var2))
  
  jointTable <- prop.table(table(var1, var2))
  plotData <- as.data.frame(jointTable)
  colnames(plotData)[1:2] <- c(var1Name,var2Name)
  
  plotData$marginVar1 <- prop.table(table(var1))
  plotData$var2Height <- plotData$Freq / plotData$marginVar1
  plotData$var1Center <- c(0, cumsum(plotData$marginVar1)[1:levVar1 -1]) +  plotData$marginVar1 / 2
  
  p <- ggplot(plotData, aes(var1Center, var2Height))
  p <- p + geom_bar(stat = "identity", aes(width = marginVar1, fill = var2), col = "Black")
  p <- p + geom_text(aes(label = as.character(var1), x = var1Center, y = 1.05, angle = 30))
  p <- p + xlab(var1Name) + ylab(var2Name)
  p
}
