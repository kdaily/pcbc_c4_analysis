getGeneSetsFromSynapse <- function(FolderID){
  
  downloadGeneSets <- function(id){
    # Get all the gene sets and read as lists
    fp = file(synGet(id)@filePath)
    GeneList = strsplit(readLines(fp),'\t')
    close(fp)
    
    # Get all the names of the gene list
    names(GeneList) = sapply(GeneList, function(x){return(x[1])})
    
    # Remove the first two elements of all lists
    GeneList = lapply(GeneList, function(x){return(x = x[-(1:2)])})
    
    return(GeneList)
  }

  All.Files = synQuery(paste0('select name,id from file where parentId == "',FolderID,'"'))
  GeneLists = lapply(All.Files$file.id, downloadGeneSets)
 
  names(GeneLists) = tools::file_path_sans_ext(All.Files$file.name)

  return(GeneLists)
}
