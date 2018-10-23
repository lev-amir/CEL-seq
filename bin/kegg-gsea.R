### decided not to perform ###

source("https://bioconductor.org/biocLite.R")
biocLite("clusterProfiler") ; library(clusterProfiler)
biocLite("org.Ce.eg.db") ; library(org.Ce.eg.db)

a <- bitr(NT$Gene_stable_ID, fromType="WORMBASE", toType="ENTREZID", OrgDb="org.Ce.eg.db", drop=T)
a <- a %>% dplyr::rename(Gene_stable_ID = WORMBASE)
c <- left_join(NT,a, by = "Gene_stable_ID")
gsea <- c %>% arrange(desc(log2FoldChange)) %>% .$log2FoldChange
names(gsea) <- c %>% arrange(desc(log2FoldChange)) %>% .$ENTREZID

## KEGG Gene Set Enrichment Analysis ##
kk2 <- gseKEGG(geneList     = gsea,
                keyType      = 'ncbi-geneid',
                organism     = 'cel',
                nPerm        = 1000,
                minGSSize    = 40,
                pvalueCutoff = 0.05,
                verbose      = TRUE)
head(kk2)

## Results, due to minGSSize    = 40, other values remove some of the results ##
# ID                    Description setSize enrichmentScore       NES      pvalue   p.adjust    qvalues rank
# cel03010 cel03010                       Ribosome     125       0.6713444  2.237956 0.001106195 0.03429204 0.02794597 2744
# cel03015 cel03015      mRNA surveillance pathway      69      -0.5240880 -1.999406 0.005952381 0.04696970 0.03827751 3162
# cel00240 cel00240          Pyrimidine metabolism      72      -0.3799745 -1.447325 0.006289308 0.04696970 0.03827751 4231
# cel04120 cel04120 Ubiquitin mediated proteolysis      80      -0.4375823 -1.713697 0.006535948 0.04696970 0.03827751 3425
# cel03013 cel03013                  RNA transport     114      -0.3601193 -1.526163 0.008547009 0.04696970 0.03827751 4167
# cel03040 cel03040                    Spliceosome     111      -0.5156142 -2.145414 0.009090909 0.04696970 0.03827751 3558
# leading_edge
# cel03010 tags=64%, list=23%, signal=50%
# cel03015 tags=51%, list=26%, signal=38%
# cel00240 tags=54%, list=35%, signal=35%
# cel04120 tags=48%, list=28%, signal=34%
# cel03013 tags=50%, list=34%, signal=33%
# cel03040 tags=65%, list=29%, signal=47%