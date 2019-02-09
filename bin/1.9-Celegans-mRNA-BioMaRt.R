## Obtain all the C. elegans protein coding genes from Biomart

#+ Load bioconductor and packages ----
source("https://bioconductor.org/biocLite.R")
if (!require(biomaRt)) {biocLite("biomaRt"); library(biomaRt)} # Ensembl biomaRt to obtain C. elegans sequences.
if (!require(tidyverse)) {install.packages("tidyverse"); library(tidyverse)}

#'Select a BioMart database and celegans dataset. ----
#'Shortcut command designed from the following:
# listMarts() # check which BioMart web services are available
# ensembl=useMart("ENSEMBL_MART_ENSEMBL") # choose to query the Ensembl BioMart database
# listDatasets(ensembl) # look at which datasets are available in the selected BioMart
# ensembl = useDataset("celegans_gene_ensembl",mart=ensembl) # choose to use the celegans dataset
#+ Select a BioMart database and celegans dataset
ensembl = biomaRt::useMart("ENSEMBL_MART_ENSEMBL",dataset="celegans_gene_ensembl")

#' Obtain gene IDs to all protein coding genes in the C. elegans biomart ----
#' Shortcut command designed from the following:
# listFilters(ensembl) # browse filter options
# attributes = listAttributes(ensembl) # values we are interested in to retrieve
# attributePages(ensembl) # browse attributes categories
# attributes[attributes$page == 'sequences',] # browse attributes for sequences
mRNA.ids <- getBM( # gene ids to all protein coding genes in celegans biomart
  attributes = c("ensembl_gene_id","external_gene_name","go_id"), 
  filters = "biotype", # focus on protein-coding genes as only these are identified by CEL-seq
  values = "protein_coding", 
  mart = ensembl)

save(ensembl,mRNA.ids, file = "data/celegans-biomart-database-mRNAs.Rdata")