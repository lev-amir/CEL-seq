source("https://bioconductor.org/biocLite.R")
if (!require(tidyverse)) {install.packages("tidyverse"); library(tidyverse)}
if (!require(biomaRt)) {biocLite("biomaRt"); library(biomaRt)}


## load Celegans dataset
ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="celegans_gene_ensembl")

## locate filters of interest
# filters = listFilters(ensembl)
# 215                       interpro                                     Interpro ID(s) [e.g. IPR000001]
# 219                           pfam                                    Pfam domain ID(s) [e.g. PF00001]
# 224                         pfscan                               PROSITE profiles ID(s) [e.g. PS01031]

## genes of interest as filtered by 'F-box' domain logged in different databases 
pfam <- # Pfam Family: F-box (PF00646)
  getBM(attributes=c('ensembl_gene_id'), 
        filters = 'pfam',
        values = 'PF00646',
        mart = ensembl)
pfscan <- # Prosite: F-box domain profile FBOX, PS50181
  getBM(attributes=c('ensembl_gene_id'), 
        filters = 'pfscan',
        values = 'PS50181',
        mart = ensembl)
interpro <- # Interpro F-box domain (IPR001810)
  getBM(attributes=c('ensembl_gene_id'), 
        filters = 'interpro',
        values = 'IPR001810',
        mart = ensembl)

## uniq F-box gene IDs from the three databases combined
fbox.genes <- Reduce(union, list(pfam$ensembl_gene_id, pfscan$ensembl_gene_id, interpro$ensembl_gene_id)) 

## locate attributes of interest
# attributes = listAttributes(ensembl)
# attributePages(ensembl) # list the different attribute catergories
# 1                   ensembl_gene_id                             Gene stable ID feature_page
# 2             ensembl_transcript_id                       Transcript stable ID feature_page
# 3                ensembl_peptide_id                          Protein stable ID feature_page
# 4                   ensembl_exon_id                             Exon stable ID feature_page
# 15               external_gene_name                                  Gene name feature_page
# 17         external_transcript_name                            Transcript name feature_page
# 21                     gene_biotype                                  Gene type feature_page
# 22               transcript_biotype                            Transcript type feature_page
# 52                    wormbase_gene                           WormBase Gene ID feature_page
# 53                wormbase_gseqname             WormBase Gene Sequence-name ID feature_page
# 54                   wormbase_locus                          WormBase Locus ID feature_page
# 55              wormbase_transcript                     WormBase Transcript ID feature_page
# 56                       wormpep_id                                 Wormpep ID feature_page
# 79                             pfam                             Pfam domain ID feature_page
# 91                      scanprosite                        PROSITE patterns ID feature_page
# 94                           pfscan                        PROSITE profiles ID feature_page
# 106                        interpro                                Interpro ID feature_page

## find attributes associated to homologes in humans
# homologs <- attributes[attributes$page == "homologs",]
# 890                 hsapiens_homolog_ensembl_gene                           Human gene stable ID homologs
# 891         hsapiens_homolog_associated_gene_name                                Human gene name homologs

# homologs.name <- homologs$name
# i <- grep(".*hsapien.*",tolower(homologs.name))
# homologs[i,]$name

## human homologs attributes of interest
human.filters <- c('hsapiens_homolog_ensembl_gene',
                   'hsapiens_homolog_associated_gene_name',
                   'hsapiens_homolog_orthology_type',
                   'hsapiens_homolog_perc_id',
                   'hsapiens_homolog_perc_id_r1',
                   'hsapiens_homolog_orthology_confidence')

## data-frame w/ Celegans F-BOX genes and information of their human homologs
fbox.human <- getBM(attributes = c('ensembl_gene_id','external_gene_name',human.filters), 
                  filters = 'ensembl_gene_id',
                  values = fbox.genes,
                  mart = ensembl)

## export data-frame to file
write.csv(fbox.human,"fbox.human.csv")

## genes of interest, f-boxes
fbox <- read_csv("fbox\\fbox.human.csv")
glimpse(fbox)

## Combine data-frames
NT.fbox <- NT %>%
  filter(Gene_name %in% fbox[["external_gene_name"]] | Transcript_stable_ID %in% fbox[["ensembl_gene_id"]]) %>% # select genes of interest
  left_join(fbox, by = c("Gene_name" = "external_gene_name"))# combine columns of two df  
  
## plot expression fold-change
NT.fbox %>% ggplot(aes(x = log2FoldChange)) + geom_histogram()

## gene that passed the p-value threshold
glimpse(NT.fbox)
NT.fbox %>% filter(padj < 0.05) %>% 
  select(Gene_name,
         Gene_description,
         log2FoldChange,
         hsapiens_homolog_ensembl_gene,
         hsapiens_homolog_associated_gene_name,
         hsapiens_homolog_orthology_confidence,
         GO_term_name) %>% 
  View()

## test if there's a significant enrichment of f-box genes that had a statistically significant log2FC
not.fbox.sig <- NT %>%
  filter( !(Gene_name %in% fbox[["external_gene_name"]] | Transcript_stable_ID %in% fbox[["ensembl_gene_id"]]) ) %>% # select genes of interest
  left_join(fbox, by = c("Gene_name" = "external_gene_name")) %>%
  filter(padj < 0.05) %>% nrow()

not.fbox.notsig <- NT %>%
  filter( !(Gene_name %in% fbox[["external_gene_name"]] | Transcript_stable_ID %in% fbox[["ensembl_gene_id"]]) ) %>% # select genes of interest
  left_join(fbox, by = c("Gene_name" = "external_gene_name")) %>%
  filter(padj > 0.05) %>% nrow()

yes.fbox.Sig <- NT %>%
  filter( (Gene_name %in% fbox[["external_gene_name"]] | Transcript_stable_ID %in% fbox[["ensembl_gene_id"]]) ) %>% # select genes of interest
  left_join(fbox, by = c("Gene_name" = "external_gene_name")) %>%
  filter(padj < 0.05) %>% nrow()

yes.fbox.notSig <- NT %>%
  filter( (Gene_name %in% fbox[["external_gene_name"]] | Transcript_stable_ID %in% fbox[["ensembl_gene_id"]]) ) %>% # select genes of interest
  left_join(fbox, by = c("Gene_name" = "external_gene_name")) %>%
  filter(padj > 0.05) %>% nrow()

fet = matrix(c(1,2,3,4), nrow = 2)
fet = matrix(c(yes.fbox.notSig,
               not.fbox.notsig,
               yes.fbox.Sig,
               not.fbox.sig), nrow = 2)
fisher.test(fet)

