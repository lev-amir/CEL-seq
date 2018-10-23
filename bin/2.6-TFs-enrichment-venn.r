#+ Load previously generated data of transcription factors enrichment in each DEG set.
load("experiments/2018-09-02-TFs-motifs/data/TFs-motif-enrichment-results")
load("experiments/2018-09-02-TFs-motifs/data/TFs-motif-enrichment-results-random")

## Visualize p-values ----

#+ Generate a report and extract p-values for motif pair.
report <- lapply(res, groupReport,by.top.motifs=FALSE) # Report of all enrichment scores.
report.rand <- lapply(rand.res,groupReport,by.top.motifs=FALSE)
report.ps <- unlist(lapply(report, function(x) {x$p.value}))
report.rand.ps <- unlist(lapply(report.rand, function(x) {x$p.value}))

#+ Plot pooled p-values for each motif, to check if there's a difference between true data and randomly sampled data.
par(mfrow = c(1,2))
hist(report.ps,20)
hist(report.rand.ps,20)

#' Plot transcription factor motifs with statistically significantly (p-value < 0.05) enriched motifs. 
#' **Top motifs** indicates in what percentage of sequences was the corresponding motif among the 
#' 5% lowest p-value (out of all examined motifs p-values).
plot(report$N2.up[report$N2.up$p.value < 0.05], fontsize=7, id.fontsize=6)

x <- report

a <- x$N2.up[x$N2.up$p.value<0.05,]
b <- x$AM.up[x$AM.up$p.value<0.05,]
ab <- intersect(a$target,b$target)

c <- x$N2.dn[x$N2.dn$p.value<0.05,]
d <- x$AM.dn[x$AM.dn$p.value<0.05,]
cd <- intersect(c$target,d$target)

abcd <- intersect(ab,cd)
ab.only <- setdiff(ab,cd)
cd.only <- setdiff(cd,ab)
```

# Which transcription factors are transcriptionally regulated by anc-1 and have an enrichemnt of motifs in the promoters of ANC-1 regulated genes?
```{r}
motif.ids <- report$N2.up$id # all motifs tested
deg.id <- BOTH %>% filter(Gene_stable_ID %in% DEG$UNION.reg) %>% pull(Gene_name) # all differentially regulated genes

motif.CELseq <- intersect(motif.ids,deg.id) # differentially regulated motifs tested TFs


motif.sig <- unique(c(a$target,b$target,c$target,d$target)) # significanlly regulated motifs that appear in any sample
motif.sig.deg <- intersect(motif.CELseq,motif.sig) # significanlly regulated motifs that appear in any sample and whose expression is differentialy regulated

motif.sig.up <- intersect(motif.CELseq,ab)
motif.sig.dn <- intersect(motif.CELseq,cd)
motif.sig.deg[ !( motif.sig.deg %in% c(motif.sig.up,motif.sig.dn) ) ]
```

# Venn of transcription factors based on  
```{r}
# ---------------------------------------------------
venn.diagram(
x = list("In upregulated promoters" = ab, "In downregulated promoter" = cd),
filename = "experiments/2018-09-02-TFs-motifs/results/venn_TF_regulators.tiff",
resolution = 300,
imagetype = "tiff",
main = "Transcription factor motifs",
main.cex = 3,
main.fontfamily = "sans",
fill = c("red2", "green2"),
cex = 2.5 ,# size of the areas’ labels
fontfamily = "sans",# the fontfamily of the areas’ labels
cat.cex = 3 ,#the size of the category names
cat.pos = c(180+30,180-30),
cat.dist = +0.04,
cat.col = c("red2", "green2"),
cat.fontfamily = "sans",
print.mode = c("raw","percent"),
sub = "Transcription factors with motifs enriched in anc-1 regulated genes' promoters",
sub.cex = 1.5,
sub.fontfamily = "sans"
)
```

Export transcription factors lists according to their motifs enrichment in upregulated / downregulated genes' promoters.
```{r}
write(ab.only,"experiments/2018-09-02-TFs-motifs/results/TF-motifs-in-up.txt")
write(cd.only,"experiments/2018-09-02-TFs-motifs/results/TF-motifs-in-dn.txt")
write(abcd,"experiments/2018-09-02-TFs-motifs/results/TF-motifs-in-up-dn.txt")
```