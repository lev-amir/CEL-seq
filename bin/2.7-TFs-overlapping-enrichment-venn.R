#+ Load previously generated data of transcription factors enrichment in each DEG set.
load("experiments/2018-09-02-TFs-motifs/data/TFs-motif-enrichment-results")
load("experiments/2018-09-02-TFs-motifs/data/TFs-motif-enrichment-results-random")


#' Focus on the intersect of the two strains, with upregulated and downregulated genes analyzed separately. 
#+ Downregulated genes 
goi <- DEG$INTERSECT.dn # CEL-seq DEGs, query genes
length(goi) # 93 downregulated genes in both strains.
mRNA.seqs.promoter <- mRNA.seqs.promoter %>% 
  mutate(CELseq.hits = ensembl_gene_id %in% goi)

#+ Random sequences as a control
set.seed(168)
goi <- sample(mRNA.seqs.promoter$ensembl_gene_id,168)
mRNA.seqs.promoter <- mRNA.seqs.promoter %>% 
  mutate(CELseq.hits = ensembl_gene_id %in% goi)
useBigMemoryPWMEnrich(TRUE) # using large amount of memory
res = motifEnrichment(prs[mRNA.seqs.promoter$CELseq.hits], bg.custom)
useBigMemoryPWMEnrich(FALSE) # disable using large amount of memory
report <- groupReport(res,by.top.motifs=TRUE) # Report of all enrichment scores.
plot(report[report$p.value < 0.05], fontsize=7, id.fontsize=6)