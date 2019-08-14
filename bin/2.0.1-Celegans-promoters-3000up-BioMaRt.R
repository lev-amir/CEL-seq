#+ Load bioconductor and packages ----
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
require(BiocManager)

if (!require(biomaRt)) {BiocManager::install("biomaRt"); library(biomaRt)} # Ensembl biomaRt to obtain C. elegans sequences.
if (!require(Biostrings)) {BiocManager::install("Biostrings"); library(Biostrings)}
if (!require(tidyverse)) {install.packages("tidyverse"); library(tidyverse)}

# Obtain promoters of C. elegans genes-----------------------------------------
load("data/celegans-biomart-database-mRNAs.Rdata")

# Define the promoter region, relative to transcription start-site (TSS)-------
# 3000 bp upstream, and 0 bp downstream to TSS. 
# Obtain TSS upstream and downstream sequences based on ids, seperately.

# sequences upstream to TSS of all protein coding genes in celegans biomart
mRNA.seqs.TSS.up <- biomaRt::getSequence(
  id=mRNA.ids$ensembl_gene_id, # list of gene identifiers
  type="ensembl_gene_id", # The type of identifier used
  seqType="gene_flank", # gives the flanking region of the gene excluding the UTRs
  upstream=3000, # value for the upstream or downstream attribute 
  mart = ensembl)

# convert promoter sequences dataframe to DNAStringSet-------------------------
# DNAStringSet is a common input to the motif analysis packages. 
# Associate gene id to each promoter sequence
prs <- Biostrings::DNAStringSet(x = mRNA.seqs.TSS.up$gene_flank)
names(prs) <- mRNA.seqs.TSS.up$ensembl_gene_id
save(mRNA.seqs.TSS.up,prs, file = "data/Celegans-promoters-3000up.rdata")
