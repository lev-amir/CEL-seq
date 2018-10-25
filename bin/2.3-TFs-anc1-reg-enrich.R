#+ Load packages and previously generated data of transcription factors enrichment in each DEG set.
source("https://bioconductor.org/biocLite.R")
if (!require(tidyverse)) {install.packages("tidyverse"); library(tidyverse)}
if (!require(PWMEnrich)) {install.packages("PWMEnrich"); library(PWMEnrich)}
load("data/Celegans-promoters.rdata")
load("experiments/2018-10-25-TFs-anc1-reg-motifs/data/bkgd-TF-anc1-reg.rdata")

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
  write_csv(path = "experiments/2018-10-25-TFs-anc1-reg-motifs/results/TFs-motifs-enrichment-overlapping-downregulated.csv")

#' Export the reports as PDFs of the enriched transcription factors -- no need to perform
# pdf(file = "experiments/2018-09-02-TFs-motifs/results/TFs-motifs-enrichment-overlapping-downregulated.pdf",
#     height = 10.69,
#     width = 7.27,
#     paper = "a4") 
# plot(report.dn, fontsize=7, id.fontsize=6, header.fontsize=7)
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
  write_csv(path = "experiments/2018-10-25-TFs-anc1-reg-motifs/results/TFs-motifs-enrichment-overlapping-upregulated.csv")

#' Export the reports as PDFs of the enriched transcription factors -- no need to perform
# pdf(file = "experiments/2018-09-02-TFs-motifs/results/TFs-motifs-enrichment-overlapping-upregulated.pdf",
#     height = 10.69,
#     width = 7.27,
#     paper = "a4") 
# plot(report.up, fontsize=7, id.fontsize=6)
# dev.off()

#+ Differentialy regulated genes 
goi <- DEG$INTERSECT.reg # CEL-seq DEGs, query genes
length(goi) # 438 downregulated genes in both strains (w/ p-value < 0.1)
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
  write_csv(path = "experiments/2018-10-25-TFs-anc1-reg-motifs/results/TFs-motifs-enrichment-overlapping-regulated.csv")

#' Export the reports as PDFs of the enriched transcription factors -- no need to perform
# svg(filename = "experiments/2018-10-25-TFs-anc1-reg-motifs/results/TFs-motifs-enrichment-overlapping-regulated.svg",
#     width = 4, height = 1)
# plot(report.reg[report.reg$p.value < 0.1], fontsize=7, id.fontsize=6)
# dev.off()

#' As a control to the above analyses - repeat it using random promoter sequences
#+ Select the same amount of random *C. elegans* promoter sequences as a control
set.seed(512)
goi <- sample(mRNA.seqs.promoter$ensembl_gene_id,length(goi))
mRNA.seqs.promoter <- mRNA.seqs.promoter %>% mutate(CELseq.hits = ensembl_gene_id %in% goi)

useBigMemoryPWMEnrich(TRUE) # using large amount of memory
res.rnd = motifEnrichment(prs[mRNA.seqs.promoter$CELseq.hits], bg.custom)
useBigMemoryPWMEnrich(FALSE) # disable using large amount of memory

report.rnd <- groupReport(res.rnd,by.top.motifs=FALSE) # Report of all enrichment scores.
plot(report.rnd, fontsize=7, id.fontsize=6)

save("report.up","report.dn","report.reg", "res.dn", "res.up", "res.reg",
     file = "experiments/2018-10-25-TFs-anc1-reg-motifs/data/TFs-anc1-reg-enriched.Rdata")
