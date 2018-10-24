#+ Load packages and previously generated data of transcription factors enrichment in each DEG set. ----
source("https://bioconductor.org/biocLite.R")
if (!require(tidyverse)) {install.packages("tidyverse"); library(tidyverse)}
if (!require(PWMEnrich)) {install.packages("PWMEnrich"); library(PWMEnrich)}
if (!require(biomaRt)) {biocLite("biomaRt"); library(biomaRt)} 
load("experiments/2018-09-02-TFs-motifs/data/TFs-enriched-DEG.Rdata")
#' Load previosly obtained C. elegansensembl database
load("data/celegans-biomart-database-mRNAs.Rdata")

#' Fetch C. elegnas transcription factors based on GO terms.
# GO term accession   GO term name
# GO:0003700          DNA-binding transcription factor activity
# GO:0006355          regulation of transcription, DNA-templated
TFs <- getBM( # gene ids to all protein coding genes in celegans biomart
  attributes = c("ensembl_gene_id","external_gene_name"), 
  filters = c("biotype","go"), 
  values = list(biotype = "protein_coding",
                go = c("GO:0006355","GO:0003700")), 
  mart = ensembl)

#' load the 

TFs


