#+ Load bioconductor and packages ----
source("https://bioconductor.org/biocLite.R")
if (!require(biomaRt)) {biocLite("biomaRt"); library(biomaRt)} # Ensembl biomaRt to obtain C. elegans sequences.
if (!require(tidyverse)) {install.packages("tidyverse"); library(tidyverse)}

# Obtain promoters of C. elegans genes-----------------------------------------
load("data/celegans-biomart-database-mRNAs.Rdata")

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
save(mRNA.seqs.promoter,prs, file = "data/Celegans-promoters.rdata")