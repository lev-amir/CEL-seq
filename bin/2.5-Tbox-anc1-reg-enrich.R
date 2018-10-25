#+ Load packages and previously generated data of T-box transcriptions factors enrichment in each DEG set.
source("https://bioconductor.org/biocLite.R")
if (!require(tidyverse)) {install.packages("tidyverse"); library(tidyverse)}
if (!require(PWMEnrich)) {install.packages("PWMEnrich"); library(PWMEnrich)}
load("data/Celegans-promoters.rdata")
load("experiments/2018-07-11-tbx-motifs/data/backgound-celegans-proteincoding-tbx.rdata")


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
  write_csv(path = "experiments/2018-07-11-tbx-motifs/results/Tbox-motifs-enrichment-overlapping-downregulated.csv")

#' Export the reports 
pdf(file = "experiments/2018-07-11-tbx-motifs/results/Tbox-motifs-enrichment-overlapping-downregulated.pdf",
    height = 1.5,
    width = 5)
plot(report.dn, fontsize=7, id.fontsize=6)
dev.off()

#+ Upegulated genes 
goi <- DEG$INTERSECT.up # CEL-seq DEGs, query genes
length(goi) # 402 downregulated genes in both strains (w/ p-value < 0.1)
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
  mutate("enriched.motif" = report.dn$p.value < 0.05) %>%
  write_csv(path = "experiments/2018-07-11-tbx-motifs/results/Tbox-motifs-enrichment-overlapping-upregulated.csv")

#' Export the reports 
pdf(file = "experiments/2018-07-11-tbx-motifs/results/Tbox-motifs-enrichment-overlapping-upregulated.pdf",
    height = 1.5,
    width = 5)
plot(report.up, fontsize=7, id.fontsize=6)
dev.off()

#+ Differentially regulated genes 
goi <- DEG$INTERSECT.reg # CEL-seq DEGs, query genes
length(goi) # 840 downregulated genes in both strains (w/ p-value < 0.1)
mRNA.seqs.promoter <- mRNA.seqs.promoter %>% mutate(CELseq.hits = ensembl_gene_id %in% goi)

#+ Motif enrichment analysis in gene-sets
useBigMemoryPWMEnrich(TRUE) # using large amount of memory
res.reg = motifEnrichment(prs[mRNA.seqs.promoter$CELseq.hits], bg.custom)
useBigMemoryPWMEnrich(FALSE) # disable using large amount of memory

#' Generate table of transcription factors that are enriched in the promoters of anc-1 downregulated genes.
#' Export the reports as tables
report.reg <- groupReport(res.reg,by.top.motifs=FALSE)
report.reg %>% 
  as.data.frame() %>% 
  mutate("enriched.motif" = report.reg$p.value < 0.05) %>%
  write_csv(path = "experiments/2018-07-11-tbx-motifs/results/Tbox-motifs-enrichment-overlapping-regulated.csv")

#' Export the reports as PDFs of the enriched transcription factors -- no need to perform
# pdf(file = "experiments/2018-09-02-TFs-motifs/results/TFs-motifs-enrichment-overlapping-downregulated.pdf",
#     height = 10.69,
#     width = 7.27,
#     paper = "a4") 
# plot(report.dn, fontsize=7, id.fontsize=6, header.fontsize=7)
# dev.off()

save("report.up","report.dn","report.reg", "res.dn", "res.up", "res.reg",
     file = "experiments/2018-07-11-tbx-motifs/data/Tbox-anc1-reg-enriched.Rdata")
