
if (!require(PWMEnrich)) {install.packages("PWMEnrich"); library(PWMEnrich)}
load("experiments/2018-09-02-TFs-motifs/data/TFs-enriched-DEG.Rdata")


#' Create p-values matrix. Obtain the p-values of motifs per sequence.
#' - rows correspond to the different input sequences and the columns correspond to motifs.
#' **sequence.nobg** contains the raw affinity scores.
#' **sequence.bg** contains the corresponding P-values.
pvals.dn <- res.dn$sequence.bg
pvals.up <- res.up$sequence.bg
pvals <- rbind(pvals.dn,pvals.up) # Combine both up- and down-regulated genes
colnames(pvals) <- toupper(colnames(pvals)) # Change Transcription factor names to UPPER-CASE

#' Grouping variables
p <- data.frame(pvals = c(rownames(pvals.dn),rownames(pvals.up)),
                direction = c(rep("dn",length(rownames(pvals.dn))),rep("up",length(rownames(pvals.up)))))


#' Hierarchichal clustering the transcription factors based on their p-value score 
#' to the genes that are differentially regulated by anc-1 RNAi.
col_dend = hclust(dist(t(pvals)))

#' Package to color cluster trees
if (!require(dendextend)) {install.packages("dendextend"); library("dendextend")}

#' based on (Heatmap - Static and Interactive: Absolute Guide)[http://www.sthda.com/english/articles/28-hierarchical-clustering-essentials/93-heatmap-static-and-interactive-absolute-guide/] 
Heatmap(pvals,
        name = "Motif presence (p-value)", # title of legend
        split = p, # split by a vector specifying rowgroups
        row_title = "Genes regulated by anc-1 RNAi",
        column_names_side = "top",
        show_row_names = FALSE,
        show_row_dend = FALSE,
        row_names_gp = gpar(fontsize = 6), # Text size for row names
        column_names_gp = gpar(fontsize = 6), # Text size for column names
        cluster_columns = color_branches(col_dend, k = 3))


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



