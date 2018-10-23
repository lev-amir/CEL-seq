### Gene Ontology (GO) term analysis ###
source("https://bioconductor.org/biocLite.R")
biocLite("clusterProfiler") ; library(clusterProfiler)
biocLite("biomaRt") ; library(biomaRt)
biocLite("GOstats") ; library(GOstats)
biocLite("topGO") ; library(topGO)
biocLite("org.Ce.eg.db") ; library(org.Ce.eg.db)
biocLite("AnnotationHub") ; library(AnnotationHub)
### clusterProfiler GO classification of gene lists ###

# AM.dn.GO <- groupGO(gene = both.sets$AM.dn,
#                     keytype  = 'ENSEMBL',
#                     OrgDb    = org.Ce.eg.db,
#                     ont      = "MF",
#                     level    = 2,
#                     readable = TRUE)
# head(AM.dn.GO)
# ## need to remove anc-1

AnnotationHub()
ahs2 <- query(ah, c("Caenorhabditis elegans","Ensembl","WBcel235"))
query(ah, "Caenorhabditis_elegans.WBcel235.88.gtf")
cel.ensembl <- ah[['AH53484']]
names(cel.ensembl)
cel.ensembl$gene_id
unique(cel.ensembl$type) # type == "gene"
unique(cel.ensembl$gene_biotype) # gene_biotype == "protein_coding"
cel.proteinCoding.genes <- cel.ensembl$gene_id[cel.ensembl$gene_biotype == "protein_coding" & cel.ensembl$type == "gene"]


### clusterProfiler GO over-representation test ###
AM.dn.enrichGO <- enrichGO(gene          = both.sets$AM.dn,
                            keytype       = 'ENSEMBL',
                            # universe      = cel.proteinCoding.genes,
                            OrgDb         = org.Ce.eg.db,
                            ont           = "BP",
                            pAdjustMethod = "BH",
                            #pvalueCutoff  = 0.05,
                            #qvalueCutoff  = 0.05,
                            readable      = TRUE,
                            pool = TRUE)
head(AM.dn.enrichGO, n=15)
dropGO() # to remove specified terms which are too general.
gofilter() # to remove all terms that are not in specified levels.
head(clusterProfiler::simplify(AM.dn.enrichGO)) # to remove redundant GO terms

## Analyze GOrilla data
GOrilla <- read_csv("data\\GOrilla\\GOrillia_combined_qvalueLessThan0.05.csv")

GOrilla <- GOrilla %>% mutate("-log q-value" = -log(`FDR q-value`))
GOrilla %>%
  filter(strain == "N2") %>%
  ggplot(aes(x = reorder(Description, `-log q-value`), y = `-log q-value`, group = class, fill = class)) + 
  geom_col() + 
  coord_flip() + 
  xlab("Gene Ontology (GO) description")


ensembl = useMart("ensembl",dataset="celegans_gene_ensembl")
filter <- listFilters(ensembl)
attributes <- listAttributes(ensembl,page="feature_page",what = c("name"))
# for (i in 1:length(attributes)) {
#   entrez <- getBM(
#     attributes = attributes[i],
#     filters = "ensembl_gene_id",
#     values = "WBGene00004475",
#     mart = ensembl
#   )
#   print(entrez)
# }

entrez <- getBM(
  attributes = c("agilent_gpl13914","ucsf_gpl9450"),
  filters = "ensembl_gene_id",
  values = both.pval.FC[["Gene_stable_ID"]],
  mart = ensembl
)
entrezJoined <- union(entrez$agilent_gpl13914,entrez$ucsf_gpl9450)
CELE_entrezJoined <- paste0('CELE_',entrezJoined)

# Create a matrix of gene log2 fold changes
gene_matrix <- both.pval.FC[["log2FoldChange"]]
# Add the entrezID's as names for each logFC entry
names(gene_matrix) <- CELE_entrezJoined

kegg_enrich <- enrichKEGG(gene = names(gene_matrix),
                          organism = 'celegans',
                          pvalueCutoff = 0.05, 
                          qvalueCutoff = 0.10)

# Plot results
barplot(kegg_enrich, 
        drop = TRUE, 
        showCategory = 10, 
        title = "KEGG Enrichment Pathways",
        font.size = 8)

biocLite("org.Ce.eg.db") ; library(org.Ce.eg.db)

go_enrich <- enrichGO(gene = both.pval.FC[["NCBI_gene_ID"]],
                      OrgDb = 'org.Ce.eg.db', 
                      readable = T,
                      ont = "CC",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)

# Plot results
barplot(go_enrich, 
        drop = TRUE, 
        showCategory = 10, 
        title = "GO Biological Pathways",
        font.size = 8)


### GSEA decided not to perform ###
source("https://bioconductor.org/biocLite.R")
biocLite("clusterProfiler") ; library(clusterProfiler)
biocLite("org.Ce.eg.db") ; library(org.Ce.eg.db)

enrichGO <- gseGO(geneList     = gsea.N2,
                  OrgDb        = org.Ce.eg.db,
                  keyType      = "ENSEMBL",
                  ont          = "BP",
                  nPerm        = 1000,
                  minGSSize    = 25,
                  maxGSSize    = 500,
                  pvalueCutoff = 0.05,
                  verbose      = TRUE)
enrichGO@result %>% select(2:9)
### no term enriched under specific pvalueCutoff... ###

gseGO(geneList     = gsea.N2, OrgDb = org.Ce.eg.db)
