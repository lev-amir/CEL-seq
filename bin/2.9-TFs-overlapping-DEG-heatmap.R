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

#' Create p-values matrix. Obtain the p-values of motifs per sequence.
#' - rows correspond to the different input sequences and the columns correspond to motifs.
#' **sequence.nobg** contains the raw affinity scores.
#' **sequence.bg** contains the corresponding P-values.
pvals.dn <- res.dn$sequence.bg
pvals.up <- res.dn$sequence.bg




#' load the 
load("data/geneExpression.rdata")



## Plot heatmap of p-values ver2
based on (Heatmap - Static and Interactive: Absolute Guide
)[http://www.sthda.com/english/articles/28-hierarchical-clustering-essentials/93-heatmap-static-and-interactive-absolute-guide/]
if (!require(dendextend)) {install.packages("dendextend"); library("dendextend")}



if (!require(ComplexHeatmap)) {biocLite("ComplexHeatmap"); library("ComplexHeatmap")}

# Clean df before heatmap
a <- BOTH %>% 
  filter(strain == "N2") %>% 
  dplyr::select(c("Gene_stable_ID","Gene_name","log2FoldChange","padj"))

row_dend = hclust(dist(pvals.mat.sig)) # row clustering
col_dend = hclust(dist(t(pvals.mat.sig))) # column clustering

hclust(dist(t(b)))

b <- as.matrix(a$log2FoldChange)

Heatmap(a)

Heatmap(pvals.mat.sig,
        name = "Motif presence (p-value < 0.05)", # title of legend
        row_title = "Genes downregulated by anc-1 RNAi",
        column_names_side = "top",
        row_names_gp = gpar(fontsize = 6), # Text size for row names
        cluster_rows = color_branches(row_dend, k = 8),
        cluster_columns = color_branches(col_dend, k = 3),
        col = c("grey","blue")) 



