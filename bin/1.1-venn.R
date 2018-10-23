source("https://bioconductor.org/biocLite.R")
install.packages("VennDiagram") ; library(VennDiagram)

# Load processed data----------------------------------------------------------
load("data/geneExpression.rdata")

## Hypergeometric distribution test to check if gene overlap is significant----
# Check the p-values of overlaps that are equal or greater than the measured overlap

# Total population are genes that appeard in the results of either AM140 or N2 CEL-seq
totalpop <- length( union(NT$Gene_stable_ID,AM$Gene_stable_ID) )

# Upregulated genes

# Genes that were upregulated in AM140 animals
sample1 <- length(DEG$AM.up)
# Genes that were upregulated in N2 animals
sample2 <- length(DEG$N2.up)
# Genes that were upregulated in both AM140 and N2 animals
overlap <- length( intersect(DEG$AM.up,DEG$N2.up) )

hyper.up <- sum( dhyper(x = overlap:sample1,
                        m = sample2,
                        n = totalpop - sample2,
                        k = sample1) )
print(hyper.up)

# Downregulated genes

# Genes that were upregulated in AM140 animals
sample1 <- length(DEG$AM.dn)
# Genes that were upregulated in N2 animals
sample2 <- length(DEG$N2.dn)
# Genes that were upregulated in both AM140 and N2 animals
overlap <- length( intersect(DEG$AM.dn,DEG$N2.dn) )

hyper.dn <- sum( dhyper(x = overlap:sample1,
                        m = sample2,
                        n = totalpop - sample2,
                        k = sample1) )
print(hyper.dn)

# Venn of up-regulated genes---------------------------------------------------
venn.diagram(
  x = list(AM140 = DEG$AM.up, N2 = DEG$N2.up),
  filename = "experiments/2018-04-24-venn/results/venn_anc1_upregulated.tiff",
  resolution = 300,
  imagetype = "tiff",
  main = "Upregulated (p < 0.1)",
  main.cex = 3,
  main.fontfamily = "sans",
  fill = c("blue", "red"),
  cex = 2.5 ,# size of the areas’ labels
  fontfamily = "sans",# the fontfamily of the areas’ labels
  cat.cex = 3 ,#the size of the category names
  cat.pos = c(180+30,180-30),
  cat.dist = +0.04,
  cat.col = c("blue", "red"),
  cat.fontfamily = "sans",
  print.mode = c("raw","percent"),
  sub = "p < 0.0001",
  sub.cex = 2,
  sub.fontfamily = "sans"
)

# Venn of down-regulated genes-------------------------------------------------
venn.diagram(
  x = list(AM140 = DEG$AM.dn, N2 = DEG$N2.dn),
  filename = "experiments/2018-04-24-venn/results/venn_anc1_downregulated.tiff",
  resolution = 300,
  imagetype = "tiff",
  main = "Downregulated (p < 0.1)",
  main.cex = 3,
  main.fontfamily = "sans",
  fill = c("blue", "red"),
  cex = 2.5 ,# size of the areas’ labels
  fontfamily = "sans",# the fontfamily of the areas’ labels
  cat.cex = 3 ,#the size of the category names
  cat.pos = c(180+30,180-30),
  cat.dist = +0.04,
  cat.col = c("blue", "red"),
  cat.fontfamily = "sans",
  print.mode = c("raw","percent"),
  sub = "p < 0.0001",
  sub.cex = 2,
  sub.fontfamily = "sans",
  #inverted = TRUE
)

# Get top 10 genes per overlap by adjusted p-values----------------------------
a <- NT$Gene_name[NT$Gene_stable_ID %in% intersect(DEG$AM.up,DEG$N2.up)][1:15]
write(a, "experiments/2018-04-24-venn/results/topGenes-up.txt")
a <- NT$Gene_name[NT$Gene_stable_ID %in% intersect(DEG$AM.dn,DEG$N2.dn)][1:15]
write(a, "experiments/2018-04-24-venn/results/topGenes-dn.txt")
