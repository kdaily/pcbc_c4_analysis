#!/usr/bin/env Rscript

# Function to perform enrichment analysis od modules (from synapse as RData file)
# Get arguments from comman line
args = commandArgs(TRUE)

# Clear R console screen output
cat("\014")

# Clear R workspace
setwd('/home/ec2-user/Work/Github/pcbc_c4_analysis/code/R')
############################################################################################################

############################################################################################################
#### Libraries ####
library(synapseClient)
library(dplyr)
library(WGCNA)
library(tools)
library(stringr)
library(igraph)
library(data.table)
library(biomaRt)
library(plyr)

# Needs the dev branch
library(rGithubClient)

# login to synapse
synapseLogin()
############################################################################################################

############################################################################################################
#### Github commit ####
# Get github links for provenance
thisFileName <- 'enrichment_mRNA_SGE.R'

# Github link
thisRepo <- getRepo(repository = "th1vairam/pcbc_c4_analysis", 
                    ref="branch", 
                    refName='enrich')

thisFile <- getPermlink(repository = thisRepo,
                        repositoryPath=paste0('code/R/', thisFileName))

# Synapse specific parameters
activityName = 'Enrichment Analysis'
activityDescription = 'Enrichment analysis based on enrichr genesets'
############################################################################################################

############################################################################################################
#### Function definitions ####
# Function to filter Gene Sets
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
                       x = x[len>minSize & len<maxSize]
                       return(x)
                     },
                     minSize,
                     maxSize)
  len = sapply(GeneLists, length)
  GeneLists = GeneLists[len != 0]
  
  return(GeneLists)
}

# Function to perform Fishers enrichment analysis
fisherEnrichment <- function(genesInSignificantSet, # A character vector of differentially expressed or some significant genes to test
                             genesInGeneSet, # A character vector of genes in gene set like GO annotations, pathways etc...
                             genesInBackground # Background genes that are 
){
  genesInSignificantSet = intersect(genesInSignificantSet, genesInBackground) # back ground filtering
  genesInNonSignificantSet = base::setdiff(genesInBackground, genesInSignificantSet)
  genesInGeneSet = intersect(genesInGeneSet, genesInBackground) # back ground filtering
  genesOutGeneSet = base::setdiff(genesInBackground,genesInGeneSet)
  
  pval = fisher.test(
    matrix(c(length(intersect(genesInGeneSet, genesInSignificantSet)),             
             length(intersect(genesInGeneSet, genesInNonSignificantSet)),
             length(intersect(genesOutGeneSet, genesInSignificantSet)),
             length(intersect(genesOutGeneSet, genesInNonSignificantSet))), 
           nrow=2, ncol=2),
    alternative="greater")
  OR = (length(intersect(genesInGeneSet, genesInSignificantSet)) * length(intersect(genesOutGeneSet, genesInNonSignificantSet))) / (length(intersect(genesInGeneSet, genesInNonSignificantSet)) * length(intersect(genesOutGeneSet, genesInSignificantSet)))
  return(data.frame(pval = pval$p.value,
                    ngenes = length(genesInGeneSet),
                    noverlap = length(intersect(genesInGeneSet, genesInSignificantSet)),
                    Odds.Ratio = OR,
                    Genes = paste(intersect(genesInGeneSet, genesInSignificantSet), collapse = '|')
  )
  )
}

# Function to convert rownames to first column of a df
rownameToFirstColumn <- function(DF,colname){
  DF <- as.data.frame(DF)
  DF[,colname] <- row.names(DF)
  DF <- DF[,c(dim(DF)[2],1:(dim(DF)[2]-1))]
  return(DF)
}
############################################################################################################

############################################################################################################
#### Get gene sets ####
# Download enrichr gene sets from synapse
GL_OBJ = synGet('syn4867851')
ALL_USED_IDs = GL_OBJ$properties$id
load(GL_OBJ@filePath)

gsets = c("Achilles_fitness_decrease", "Achilles_fitness_increase", "Allen_Brain_Atlas_down", "Allen_Brain_Atlas_up",
          "BioCarta", "CMAP_down", "CMAP_up", "Cancer_Cell_Line_Encyclopedia", "ChEA", "Cross_Species_Phenotype",
          "Disease_Signatures_from_GEO_down", "Disease_Signatures_from_GEO_up", "Drug_Perturbations_from_GEO",
          "ENCODE_Histone_Modifications_2013", "ESCAPE", "GO_Biological_Process", "GO_Cellular_Component", "GO_Molecular_Function",
          "GeneSigDB", "Genome_Browser_PWMs.1", "HMDB_Metabolites", "HomoloGene", "Human_Gene_Atlas", "KEGG_2015",
          "MGI_Mammalian_Phenotype", "MGI_Mammalian_Phenotype_Level_3", "MGI_Mammalian_Phenotype_Level_4", "MSigDB_Computational",
          "MSigDB_Oncogenic_Signatures", "Mouse_Gene_Atlas", "NCI-60_Cancer_Cell_Lines", "NCI-Nature", 
          "NURSA_Human_Endogenous_Complexome", "OMIM_Disease", "OMIM_Expanded", "PPI_Hub_Proteins", "Pfam_InterPro_Domains",
          "Phosphatase_Substrates_from_DEPOD", "Reactome", "SILAC_Phosphoproteomics", "TF-LOF_Expression_from_GEO", 
          "TargetScan_microRNA", "Tissue_Protein_Expression_from_Human_Proteome_Map", "Tissue_Protein_Expression_from_ProteomicsDB",
          "Transcription_Factor_PPIs", "Virus_Perturbations_from_GEO_down", "Virus_Perturbations_from_GEO_up", "WikiPathways_2015")
GeneSets = GeneSets = GeneSets[gsets]
############################################################################################################

############################################################################################################
#### Get differentialy expressed genes ####
# Get differential expression (ALL)
DEXP.ALL_ID = 'syn4484232'
ALL_USED_IDs = c(ALL_USED_IDs, DEXP.ALL_ID)
DEXP.ALL = data.table::fread(synGet(DEXP.ALL_ID)@filePath, data.table=F)

# Filter comparisons
Comparisons.ALL = split(DEXP.ALL, factor(DEXP.ALL$Comparison))
Comparisons.ALL = lapply(Comparisons.ALL, function(x){return(unique(x$GeneSymbol[abs(x$logFC) >= log2(1.5) & x$adj.P.value<=0.05]))})

# Get differential expression (EACH)
DEXP.EACH_ID = 'syn4485993'
ALL_USED_IDs = c(ALL_USED_IDs, DEXP.EACH_ID)
DEXP.EACH = data.table::fread(synGet(DEXP.EACH_ID)@filePath, data.table=F)

# Filter comparisons
Comparisons.EACH = split(DEXP.EACH, factor(DEXP.EACH$Comparison))
Comparisons.EACH = lapply(Comparisons.EACH, function(x){return(unique(x$GeneSymbol[abs(x$logFC) >= log2(1.2) & x$adj.P.value<=0.05]))})

Comparisons = c(Comparisons.ALL, Comparisons.EACH)
############################################################################################################

############################################################################################################
#### Background gene list ####
logCPM_ID = 'syn4483934'
ALL_USED_IDs = c(ALL_USED_IDs, logCPM_ID)
logCPM = data.table::fread(synGet(logCPM_ID)@filePath, data.table=F)
backGroundGenes = logCPM$GeneName
############################################################################################################

############################################################################################################
#### Filter gene list ####
GeneSets = filterGeneSets(GeneSets, backGroundGenes, minSize = 10, maxSize = 5000)
############################################################################################################

############################################################################################################
#### Perform enrichment analysis ####
# Perform enrichment analysis
enrichResults.mRNA = list()
comp = Comparisons[[args]]
tmp = lapply(GeneSets,
             function(x, comp, genesInBackground){
               tmp = ldply(lapply(x, fisherEnrichment, comp, genesInBackground))
               setnames(tmp, '.id', 'GeneSetName')
               return(tmp)
             },
             comp, backGroundGenes)

for (name1 in names(tmp))
  tmp[[name1]]$Category = name1

enrichResults.mRNA[[args]] = as.data.frame(rbindlist(tmp))
enrichResults.mRNA[[args]]$FDR = p.adjust(enrichResults.mRNA[[args]]$pval,'fdr')
enrichResults.mRNA[[args]]$adj.pval = p.adjust(enrichResults.mRNA[[args]]$pval,'bonferroni')
writeLines(paste0('Completed ',args))  


# Write results to file
for(name in names(enrichResults.mRNA))
  enrichResults.mRNA[[name]]$ComparisonName = name
enrichmentResults = as.data.frame(rbindlist(enrichResults.mRNA))
 
# enrichmentResults$ngenes = unlist(enrichmentResults$ngenes)
# enrichmentResults$noverlap = unlist(enrichmentResults$noverlap)
# enrichmentResults$pval = unlist(enrichmentResults$pval)
# enrichmentResults$FDR = unlist(enrichmentResults$FDR)
# enrichmentResults$adj.pval = unlist(enrichmentResults$adj.pval)
# enrichmentResults$Odds.Ratio = unlist(enrichmentResults$Odds.Ratio)

write.table(enrichmentResults, file = paste(gsub('[^[:alnum:]]','_',args),'enrichmentResults_mRNA.tsv',sep='.'), sep='\t', row.names=F)
collectGarbage()
############################################################################################################

############################################################################################################
#### Write to synapse ####
# Write results to synapse
algo = 'Fisher'
ENR_OBJ = File(paste(gsub('[^[:alnum:]]','_',args),'enrichmentResults_mRNA.tsv',sep='.'), name = paste(gsub('[^[:alnum:]]','_',args),'mRNA'), parentId = 'syn4905034')
ENR_OBJ = synStore(ENR_OBJ, 
                   executed = thisFile,
                   used = ALL_USED_IDs,
                   activityName = activityName,
                   activityDescription = activityDescription)
############################################################################################################

############################################################################################################
# Write completed files to synapse
write.table(ENR_OBJ$properties$id, file = 'CompletedEnrichmentIDs.txt', sep='\n', append=T, quote=F, col.names=F, row.names=F)
writeLines(paste('Completed',args,'and stored in',ENR_OBJ$properties$id))
############################################################################################################
