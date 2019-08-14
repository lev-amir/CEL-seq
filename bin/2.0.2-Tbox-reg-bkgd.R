if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
require(BiocManager)

if (!require(biomaRt)) {BiocManager::install("biomaRt"); library(biomaRt)} # Ensembl biomaRt to obtain C. elegans sequences.
if (!require(Biostrings)) {BiocManager::install("Biostrings"); library(Biostrings)}
if (!require(PWMEnrich)) {BiocManager::install("PWMEnrich"); library(PWMEnrich)}
if (!require(tidyverse)) {install.packages("tidyverse"); library(tidyverse)}

load("data/celegans-biomart-database-mRNAs.Rdata")
load("data/geneExpression.rdata")
load("data/Celegans-promoters-3000up.rdata")

load("experiments/2018-09-02-TFs-motifs/data/TFs-PWM.rdata")
# Select T-Box transcription factors motifs 
PWM <- PWM[grep("tbx-",names(PWM))]
# convert PWM to PCM (partially arbitrariy by multiplying by 100), rounding, and converting the matrix type to "integer". This is the required input to the enrichment calculation algorithm.
PCM <- lapply(PWM,"*",100)
PCM <- lapply(PCM,round)
PCM <- lapply(PCM, FUN=function(x) type.convert(x))

## Define background set of sequences
#Using the *PWMEnrich* package. 
#The background sequences are calculated by:
#  1. The promoter regions I defined for all protein-coding genes in *C. elegans*.
#2. The known PCM motifs of *C. elegans* TFs.
#3. The promoter sequences are split into 100bp chunks and fitted.
#This process is very long and so it's best to save the output for future sessions.


useBigMemoryPWMEnrich(TRUE) # using large amount of memory
bg.custom <- makeBackground(motifs = PCM,
                            type = "logn",
                            bg.seq = prs,
                            bg.len=100,
                            bg.source="All promoters: TSS-up-500, TSS-dn-100, split into 100bp chunks")

dir.create("experiments/2019-08-14-tbx-motifs-3000bp")
dir.create("experiments/2019-08-14-tbx-motifs-3000bp/data")
save(bg.custom,
     file = "experiments/2019-08-14-tbx-motifs-3000bp/data/backgound-celegans-proteincoding-tbx.rdata")
useBigMemoryPWMEnrich(FALSE)
