---
title: "T-Box-TFs-background: Background set of *C. elegans* T-Box transcription factors binding motifs representation in protein coding genes promoters"
author: "Amir Levine"
date: "2018-10 25"
output: md_document
---

This script checks the enrichment of the currently known position-weight matrices of *C. elegans* T-box transcription factors in all protein coding genes promoters. This object will serve as the background for the enrichment of the PWMs in specified gene lists. The TFs that were selected have a known binding motif.

Load packages and previously generated data of transcription factors enrichment in each DEG set.
```{r}
source("https://bioconductor.org/biocLite.R")
if (!require(tidyverse)) {install.packages("tidyverse"); library(tidyverse)}
if (!require(biomaRt)) {biocLite("biomaRt"); library(biomaRt)} 
if (!require(PWMEnrich)) {install.packages("PWMEnrich"); library(PWMEnrich)}
```


Load previosly obtained *C. elegans* ensembl database and the differentially expression genes sorted data.
```{r}
load("data/celegans-biomart-database-mRNAs.Rdata")
load("data/geneExpression.rdata")
load("data/Celegans-promoters.rdata")
```

# Obtain transcription factors motifs
Known Position Weight Matrices for *C. elegans* transcription factors were downloaded from the (CIS-BP Database)[http://cisbp.ccbr.utoronto.ca] (downloaded only TF info and PWMs). `"experiments\\2018-09-02-TFs-motifs\\raw\\CisBP_2018_09_02_12_51_pm.log"` contains the metadata. `TF.PWM-clean.R` wrangles `"PWM.txt"` to an object suitable for downstream analysis.


Load the PWM object.
```{r}
load("experiments/2018-09-02-TFs-motifs/data/TFs-PWM.rdata")
```

Select T-Box transcription factors motifs 
```{r}
# Subset PWM list for TBX genes
PWM <- PWM[grep("tbx-",names(PWM))]
```

convert PWM to PCM (partially arbitrariy by multiplying by 100), rounding, and converting the matrix type to "integer". This is the required input to the enrichment calculation algorithm.
```{r}
PCM <- lapply(PWM,"*",100)
PCM <- lapply(PCM,round)
PCM <- lapply(PCM, FUN=function(x) type.convert(x))
# mean(sapply(PCM,storage.mode) == "integer") # verify that all matrices are stored as an "integer" (=1)
# storage.mode(PCM[[1]]) # Check just the first matrix.
```


## Define background set of sequences
Using the *PWMEnrich* package. 
The background sequences are calculated by:
1. The promoter regions I defined for all protein-coding genes in *C. elegans*.
2. The known PCM motifs of *C. elegans* TFs.
3. The promoter sequences are split into 100bp chunks and fitted.
This process is very long and so it's best to save the output for future sessions.
```{r}
useBigMemoryPWMEnrich(FALSE) # using large amount of memory
bg.custom <- makeBackground(motifs = PCM,
                            type = "logn",
                            bg.seq = prs,
                            bg.len=100,
                            bg.source="All promoters: TSS-up-500, TSS-dn-100, split into 100bp chunks")

dir.create("experiments/2018-07-11-tbx-motifs/data")
save(bg.custom,
     file = "experiments/2018-07-11-tbx-motifs/data/backgound-celegans-proteincoding-tbx.rdata")
```