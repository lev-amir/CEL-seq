# fbxc-25
sequence = prs$WBGene00019591
# perform motif enrichment!
res = motifEnrichment(sequence, bg.custom)

report = sequenceReport(res, 1)
report


ids = c("tbx-38","tbx-43","tbx-39","tbx-40")
sel.pwms = bg.custom$pwms[ids]

scores = motifScores(sequence, sel.pwms, raw.scores=TRUE)

plotMotifScores(scores, cols=c("green", "red", "blue","yellow"))


BiocManager::install("BSgenome.Celegans.UCSC.ce11")
library(BSgenome.Celegans.UCSC.ce11)

# empirical distribution for the tbx-43 motif
bcd.ecdf = motifEcdf(sel.pwms$"tbx-43", organism=BSgenome.Celegans.UCSC.ce11, bg.seq=prs, quick=TRUE)
# find the score that is equivalent to the P-value of 1e-3
threshold.1e3 = log2(quantile(bcd.ecdf, 1 - 1e-3))
threshold.1e3
# replot only the bcd motif hits with the P-value cutoff of 1e-3 (0.001)
plotMotifScores(scores, cols="green", sel.motifs="bcd", cutoff=threshold.1e3)
# Convert scores into P-values
pvals = 1 - bcd.ecdf(scores[[1]][,"bcd"])
head(pvals)


