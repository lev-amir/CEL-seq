#+ Load packages and previously generated data of transcription factors enrichment in each DEG set.
source("https://bioconductor.org/biocLite.R")
if (!require(tidyverse)) {install.packages("tidyverse"); library(tidyverse)}
if (!require(PWMEnrich)) {install.packages("PWMEnrich"); library(PWMEnrich)}
load("experiments/2018-09-02-TFs-motifs/data/background-celegans-TF.rdata")
load("data/Celegans-promoters.rdata")

#' Focus on the intersect of the two strains, with upregulated and downregulated genes analyzed separately. 

#+ Downregulated genes 
goi <- DEG$INTERSECT.dn # CEL-seq DEGs, query genes
length(goi) # 402 downregulated genes in both strains (w/ p-value < 0.1)
mRNA.seqs.promoter <- mRNA.seqs.promoter %>% mutate(CELseq.hits = ensembl_gene_id %in% goi)

#+ Motif enrichment analysis in gene-sets
useBigMemoryPWMEnrich(TRUE) # using large amount of memory
res.dn = motifEnrichment(prs[mRNA.seqs.promoter$CELseq.hits], bg.custom)
useBigMemoryPWMEnrich(FALSE) # disable using large amount of memory

#' Generate table of transcription factors that are enriched in the promoters of anc-1 downregulated genes.
#' Export the reports as tables
report.dn <- groupReport(res.dn,by.top.motifs=FALSE)
report.dn %>% 
  as.data.frame() %>% 
  mutate("enriched.motif" = report.dn$p.value < 0.05) %>%
  write_csv(path = "experiments/2018-09-02-TFs-motifs/results/TFs-motifs-enrichment-overlapping-downregulated.csv")

#' Export the reports as PDFs of the enriched transcription factors -- no need to perform
# pdf(file = "experiments/2018-09-02-TFs-motifs/results/TFs-motifs-enrichment-overlapping-downregulated.pdf",
#     height = 10.69,
#     width = 7.27,
#     paper = "a4") 
# plot(report.dn[report.dn$p.value < 0.05], fontsize=7, id.fontsize=6, header.fontsize=7)
# dev.off()

#+ Upregulated genes 
goi <- DEG$INTERSECT.up # CEL-seq DEGs, query genes
length(goi) # 438 downregulated genes in both strains (w/ p-value < 0.1)
mRNA.seqs.promoter <- mRNA.seqs.promoter %>% mutate(CELseq.hits = ensembl_gene_id %in% goi)

#+ Motif enrichment analysis in gene-sets
useBigMemoryPWMEnrich(TRUE) # using large amount of memory
res.up = motifEnrichment(prs[mRNA.seqs.promoter$CELseq.hits], bg.custom)
useBigMemoryPWMEnrich(FALSE) # disable using large amount of memory

#' Generate table of transcription factors that are enriched in the promoters of anc-1 downregulated genes.
#' Export the reports as tables
report.up <- groupReport(res.up,by.top.motifs=FALSE)
report.up %>% 
  as.data.frame() %>% 
  mutate("enriched.motif" = report.up$p.value < 0.05) %>%
  write_csv(path = "experiments/2018-09-02-TFs-motifs/results/TFs-motifs-enrichment-overlapping-upregulated.csv")

#' Export the reports as PDFs of the enriched transcription factors -- no need to perform
# pdf(file = "experiments/2018-09-02-TFs-motifs/results/TFs-motifs-enrichment-overlapping-upregulated.pdf",
#     height = 10.69,
#     width = 7.27,
#     paper = "a4") 
# plot(report.up[report.up$p.value < 0.05], fontsize=7, id.fontsize=6, header.fontsize=7)
# dev.off()

#' As a control to the above analyses - repeat it using random promoter sequences
#+ Select the same amount of random *C. elegans* promoter sequences as a control
set.seed(168)
goi <- sample(mRNA.seqs.promoter$ensembl_gene_id,length(goi))
mRNA.seqs.promoter <- mRNA.seqs.promoter %>% mutate(CELseq.hits = ensembl_gene_id %in% goi)

useBigMemoryPWMEnrich(TRUE) # using large amount of memory
res = motifEnrichment(prs[mRNA.seqs.promoter$CELseq.hits], bg.custom)
useBigMemoryPWMEnrich(FALSE) # disable using large amount of memory

report.random <- groupReport(res,by.top.motifs=FALSE) # Report of all enrichment scores.
plot(report.random[report.random$p.value < 0.05], fontsize=7, id.fontsize=6)

save("report.up","report.dn", "res.dn", "res.up",
     file = "experiments/2018-09-02-TFs-motifs/data/TFs-enriched-DEG.Rdata")
