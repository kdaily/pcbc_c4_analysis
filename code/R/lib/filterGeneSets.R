filterGeneSets <- function(GeneLists, # List of lists
                           genesInBackground, # background set of genes
                           minSize = 10,
                           maxSize = 1000){
 GeneLists = lapply(GeneLists, 
                    function(x, genesInBackground){
                      x = lapply(x, 
                                 function(x, genesInBackground){
                                   return(intersect(x, genesInBackground))
                                 },
                                 genesInBackground)
                      return(x)
                    }, 
                    genesInBackground)
 
 GeneLists = lapply(GeneLists, 
                    function(x, minSize, maxSize){
                      len = sapply(x, length)
                      x = x[len>minSize && len<maxSize]
                      return(x)
                    },
                    minSize,
                    maxSize)
 len = sapply(GeneLists, length)
 GeneLists = GeneLists[len != 0]
 
 return(GeneLists)
}