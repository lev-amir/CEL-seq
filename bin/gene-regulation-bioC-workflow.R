# http://master.bioconductor.org/packages/release/workflows/vignettes/generegulation/inst/doc/generegulation.html

source("https://bioconductor.org/biocLite.R")
biocLite(c("MotifDb",  "GenomicFeatures", 
           "TxDb.Celegans.UCSC.ce11.refGene",
           "org.Ce.eg.db", "BSgenome.Celegans.UCSC.ce11",
           "motifStack", "seqLogo"))

library(MotifDb)
library(S4Vectors)
library(seqLogo)
library(motifStack)
library(Biostrings)
library(GenomicFeatures)
library(org.Ce.eg.db); library(BSgenome.Celegans.UCSC.ce11); library(TxDb.Celegans.UCSC.ce11.refGene)

View(query(MotifDb,"Celegans"))
query(query(MotifDb,"Celegans"),"cisbp")

motifs.SwissRegulon <- as.list(query(query(query(MotifDb,"tbx"),"Hsapiens"),"SwissRegulon"))
motifs.jaspar2018 <- query(query(query(MotifDb,"tbx"),"Hsapiens"),"jaspar2018")
motifs.HOCOMOCOv10 <- query(query(query(MotifDb,"tbx"),"Hsapiens"),"HOCOMOCOv10")

###

TBX1.SwissRegulon <- motifs.SwissRegulon[[3]]
TBX1.jaspar2018 <- motifs.jaspar2018[[6]]
TBX1.HOCOMOCOv10 <- motifs.HOCOMOCOv10[[3]]

pfm.TBX1.SwissRegulon <- new("pfm", mat=TBX1.SwissRegulon,name="TBX1-SwissRegulon")
pfm.TBX1.jaspar2018 <- new("pfm", mat=TBX1.jaspar2018,name="TBX1-jaspar2018")
pfm.TBX1.HOCOMOCOv10 <- new("pfm", mat=TBX1.HOCOMOCOv10,name="TBX1-HOCOMOCOv10")

TBX1.all <- c(pfm.TBX1.SwissRegulon, pfm.TBX1.jaspar2018, pfm.TBX1.HOCOMOCOv10)
plotMotifLogoStack(DNAmotifAlignment(TBX1.all))

###

TBX2.SwissRegulon <- motifs.SwissRegulon[[6]]
TBX2.jaspar2018 <- motifs.jaspar2018[[1]]
TBX2.HOCOMOCOv10 <- motifs.HOCOMOCOv10[[6]]

pfm.TBX2.SwissRegulon <- new("pfm", mat=TBX2.SwissRegulon,name="TBX2-SwissRegulon")
pfm.TBX2.jaspar2018 <- new("pfm", mat=TBX2.jaspar2018,name="TBX2-jaspar2018")
pfm.TBX2.HOCOMOCOv10 <- new("pfm", mat=TBX2.HOCOMOCOv10,name="TBX2-HOCOMOCOv10")

TBX2.all <- c(pfm.TBX2.SwissRegulon, pfm.TBX2.jaspar2018, pfm.TBX2.HOCOMOCOv10)
plotMotifLogoStack(DNAmotifAlignment(TBX2.all))

###

TBX15.SwissRegulon <- motifs.SwissRegulon[[1]]
TBX15.jaspar2018 <- motifs.jaspar2018[[4]]
TBX15.HOCOMOCOv10 <- motifs.HOCOMOCOv10[[1]]

pfm.TBX15.SwissRegulon <- new("pfm", mat=TBX2.SwissRegulon,name="TBX15-SwissRegulon")
pfm.TBX15.jaspar2018 <- new("pfm", mat=TBX2.jaspar2018,name="TBX15-jaspar2018")
pfm.TBX15.HOCOMOCOv10 <- new("pfm", mat=TBX2.HOCOMOCOv10,name="TBX15-HOCOMOCOv10")

TBX15.all <- c(pfm.TBX15.SwissRegulon, pfm.TBX15.jaspar2018, pfm.TBX15.HOCOMOCOv10)
plotMotifLogoStack(DNAmotifAlignment(TBX15.all))

## HOCOMOCOv10 gives odd results.
## jaspar2018 and SwissRegulon give identical results. 
## I'll go with jaspar2018.

motifs.jaspar2018

pfm.TBX1 <- motifs.jaspar2018[[6]]
## MotifDb returns a position frequency matrix (a PFM) with all columns summing to 1.0, 
## but the Biostrings matchPWM method expects a position count matrix (a PCM) with integer values. 
## Transform the frequency matrix into a count matrix
pcm.TBX1 <- round(100 * pfm.TBX1)

## Create a list of the genes
genes <- c('anc-1','unc-84','zyg-12','sun-1')
## Lookup their systematic names (this is actually the Entrez ID)
orfs <- as.character(mget(genes, org.Ce.egALIAS2EG))

## Obtain the coordinates of the transcripts for the orfs.
## These coordinates, returned in a GRangesList object, specify the start location 
## (chromosome and base pair) of every known transcript for each gene. 
grl <- transcriptsBy(TxDb.Celegans.UCSC.ce11.refGene, by="gene") [orfs]

## GenomicFeatures::getPromoterSeq calculates and returns the DNA sequence of the promoter
promoter.seqs <- getPromoterSeq(grl, Celegans, upstream=1000, downstream=0)

## Note that some restructuring is needed for us to use the results of getPromoterSeqs as an argument to matchPWM.
## we need DNA strings without that overarching by-gene-name structure, we call unlist to strip off that structure, leaving us with the desired DNAStringSet
promoter.seqs <- unlist(promoter.seqs)

## search for a match of the motif to the first of these sequences, 
## we use a 90% “min.score” when we call matchPWM. This high minimum match score
matchPWM(pcm.TBX1, promoter.seqs[[1]], "90%")

## All of the matches in the promoters  may be found at once with this command
pwm.hits <- sapply(promoter.seqs, 
                   function(pseq) 
                     matchPWM(pcm.TBX1, pseq, min.score="90%"))

## summarize the motif hits for each of the three motifs


