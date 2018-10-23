# Obtain promoters of C. elegans genes-----------------------------------------

# Load bioconductor and packages-----------------------------------------------
source("https://bioconductor.org/biocLite.R")
if (!require(tidyverse)) {install.packages("tidyverse"); library(tidyverse)}
if (!require(biomaRt)) {biocLite("biomaRt"); library(biomaRt)} 
# Ensembl biomaRt to obtain C. elegans sequences.

# Select a BioMart database and celegans dataset-------------------------------
ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="celegans_gene_ensembl")
# Shortcut command designed from the following:

## listMarts() # check which BioMart web services are available
## ensembl=useMart("ENSEMBL_MART_ENSEMBL") # choose to query the Ensembl BioMart database
## listDatasets(ensembl) # look at which datasets are available in the selected BioMart
## ensembl = useDataset("celegans_gene_ensembl",mart=ensembl) # choose to use the celegans dataset

# Obtain gene IDs to all protein coding genes in the C. elegans biomart--------
mRNA.ids <- getBM( # gene ids to all protein coding genes in celegans biomart
  attributes = c("ensembl_gene_id","external_gene_name"), 
  filters = "biotype", # focus on protein-coding genes as only these are identified by CEL-seq
  values = "protein_coding", 
  mart = ensembl)
# Shortcut command designed from the following:

# listFilters(ensembl) # browse filter options
# attributes = listAttributes(ensembl) # values we are interested in to retrieve
# attributePages() # browse attributes categories
# attributes[attributes$page == 'sequences',] # browse attributes for sequences

# Define the promoter region, relative to transcription start-site (TSS)-------
# 500 bp upstream, and 100 bp downstream to TSS. 
# Obtain TSS upstream and downstream sequences based on ids, seperately.

# sequences upstream to TSS of all protein coding genes in celegans biomart
mRNA.seqs.TSS.up <- biomaRt::getSequence(
  id=mRNA.ids$ensembl_gene_id, # list of gene identifiers
  type="ensembl_gene_id", # The type of identifier used
  seqType="gene_flank", # gives the flanking region of the gene excluding the UTRs
  upstream=500, # value for the upstream or downstream attribute 
  mart = ensembl)
# sequences downstream to TSS of all protein coding genes in celegans biomart
mRNA.seqs.TSS.dn <- biomaRt::getSequence(
  id=mRNA.ids$ensembl_gene_id, # list of gene identifiers
  type= "ensembl_gene_id", # The type of identifier used
  seqType="gene_flank", # gives the flanking region of the gene excluding the UTRs
  upstream=100, # value for the upstream or downstream attribute 
  mart = ensembl)

# Combine the upstream and downstream TSS sequences to one data-frame----------
# Concatenate the upstream and downstream TSS sequences to form the promoters.
mRNA.seqs.promoter <- merge(
  mRNA.seqs.TSS.up,mRNA.seqs.TSS.dn,
  by = 'ensembl_gene_id',
  suffixes = c(".TSS.up",".TSS.dn")) %>%
  as.tibble() %>% 
  mutate(promoter = paste(gene_flank.TSS.up,gene_flank.TSS.dn,sep="")) %>%
  merge(mRNA.ids,by = 'ensembl_gene_id')
  
#mRNA.seqs.promoter %<>% 
#  mutate(promoter = paste(gene_flank.TSS.up,gene_flank.TSS.dn,sep="")) %>%
#  merge(mRNA.ids,by = 'ensembl_gene_id')
    
  
# convert promoter sequences dataframe to DNAStringSet-------------------------
# DNAStringSet is a common input to the motif analysis packages. 
# Associate gene id to each promoter sequence
prs <- Biostrings::DNAStringSet(x = mRNA.seqs.promoter$promoter)
names(prs) <- mRNA.seqs.promoter$ensembl_gene_id
save(mRNA.seqs.promoter,prs,
     file = "data/Celegans-promoters.rdata")
