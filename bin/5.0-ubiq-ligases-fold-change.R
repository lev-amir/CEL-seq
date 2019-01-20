if (!require(tidyverse)) {install.packages("tidyverse"); library(tidyverse)}
if (!require(viridis)) {install.packages("viridis"); library(viridis)}
if (!require(egg)) {install.packages("egg"); library(egg)}

# Load the dataset of differentially expressed genes `BOTH`
load("data/2018-11-27-geneExpression-pval-thresh-0.05.rdata")
load("data/celegans-biomart-database-mRNAs.Rdata") # Dataset of gene ID conversions.

path = "experiments/2018-07-04-ubiquitin-foldchange/results/"
today = Sys.Date()

## genes of interest, related to protein ubiquitination
ubiq <- read_tsv("experiments/2018-07-04-ubiquitin-foldchange//data/KEGG_Celegans_UbiqMediateProteolysis.tab")
ubiq$complex[is.na(ubiq$complex)] <- ubiq$enzyme[is.na(ubiq$complex)] # replace NA values w/ the name of of the enzyme (E1+E2)
ubiq <- left_join(ubiq,mRNA.ids,by=c("gene" = "external_gene_name")) # Add the Ensembl gene IDs to the ubiquitin dataset
ubiq <- left_join(ubiq,BOTH,by = c("ensembl_gene_id" = "Gene_stable_ID")) # Annotate the gene expression dataset with the ubiquitin labels

# N2 ----
ubiq.N2 <- ubiq %>% 
  filter(strain == 'N2') %>%
  arrange(complex,log2FoldChange) %>% # sort by log2FC and add the plot row labels
  mutate(gene = factor(.[["gene"]], levels=unique(.[["gene"]])))

## Plot E1 + E2 Ubiquitin Ligases expression fold-change
pltE1E2 <- 
  ubiq.N2 %>% 
  filter(enzyme != "E3 ubiquitin ligase") %>%
  ggplot(aes(x = gene,
             y = log2FoldChange,
             ymin = log2FoldChange - lfcSE,
             ymax = log2FoldChange + lfcSE,
             fill = complex)) + 
  geom_bar(stat = "identity", ## Create black outline bars
           position = "dodge",
           color = "black",
           size = 0.25,
           show.legend=FALSE) + 
  ylim(-3,3) +
  facet_grid(~ complex, ## Separate the bar plot by complex/enzyme
             scales="free_x", 
             space="free_x") +
  geom_text(aes(label=padj_stars, # draw significance stars
                y = log2FoldChange + sign(log2FoldChange) * (lfcSE + 0.5)), # position significance stars by edges of error bars
            angle = 90, # rotate significance stars
            vjust = "center",
            size = 3,
            nudge_x = 0.1) + ## Add minor adjustment to the stars
  geom_errorbar(width = 0.25, size = 0.25) + ## Add error bars
  xlab("") + ylab(expression(paste(log[2]," fold-change"))) + ## Remove X-label and add Y-label
  scale_fill_viridis_d(option = "C") + ## Change bar colors
  theme_light() + ## Change to white backgroud
  theme(text = element_text(size=7),
        axis.text.x = element_text(angle = 60, # Format x-axis, Angle
                                   hjust = 1, # bottom of axis
                                   face = "italic"), # italics
        strip.text.x = element_text(color = "black", face = "bold"),  # Format Facet title
        strip.background = element_rect(colour=NA,fill=NA), # Remove facet bkrng
        panel.border=element_rect(colour="grey50"))  # Format Facet remove outline

ggsave(filename = paste0(today,"-ubiquitin-foldchange-N2-E1E2-raw.pdf"),
       plot = pltE1E2,
       device = "pdf",
       path = path,
       width = 8.5,
       height = 4.5,
       units = "cm")

## Plot E3 Ubiquitin Ligases expression fold-change
pltE3 <- 
  ubiq.N2 %>% 
  filter(enzyme == "E3 ubiquitin ligase") %>%
  ggplot(aes(x = gene,
             y = log2FoldChange,
             ymin = log2FoldChange - lfcSE,
             ymax = log2FoldChange + lfcSE,
             fill = complex)) + 
  geom_bar(stat = "identity", ## Create black outline bars
           position = "dodge",
           color = "black",
           size = 0.25,
           show.legend=FALSE) + 
  ylim(-3,+3) + 
  facet_grid(~ complex, ## Separate the bar plot by complex/enzyme
             scales="free_x", 
             space="free_x") +
  geom_text(aes(label=padj_stars, # draw significance stars
                y = log2FoldChange + sign(log2FoldChange) * (lfcSE + 0.75)), # position significance stars by edges of error bars
            angle = 90, # rotate significance stars
            vjust = "center",
            size = 3,
            nudge_x = 0.2) + ## Add minor adjustment to the stars
  geom_errorbar(width = 0.25, size = 0.25) + ## Add error bars
  xlab("") + ylab(expression(paste(log[2]," fold-change"))) + ## Remove X-label and add Y-label
  scale_fill_viridis_d() + ## Change bar colors
  theme_light() + ## Change to white backgroud
  theme(text = element_text(size=7),
        axis.text.x = element_text(angle = 60, # Format x-axis, Angle
                                   hjust = 1, # bottom of axis
                                   face = "italic"), # italics
        strip.text.x = element_text(color = "black",angle = 0),  # Format Facet title
        strip.background = element_rect(colour=NA,fill=NA), # Remove facet bkrng
        panel.border=element_rect(colour="grey50"))  # Format Facet remove outline

ggsave(filename = paste0(today,"-ubiquitin-foldchange-N2-E3-raw.pdf"),
       plot = pltE3,
       device = "pdf",
       path = path,
       width = 16.7,
       height = 5,
       units = "cm")

## Align the two plots vertically, is a 1:2 height ratio.
labels = c("","E3 ubiquitin ligases")

cairo_pdf(filename = paste0(path,today,"-ubiquitin-foldchange-N2-raw.pdf"),
          width = 6.574803,
          height = 3.54331,
          pointsize = 7)
ggarrange(pltE1E2,pltE3,labels=labels,heights=2:3)
dev.off()


# AM140 ----
ubiq.AM <- ubiq %>% 
  filter(strain == 'AM140') %>%
  arrange(complex,log2FoldChange) %>% # sort by log2FC and add the plot row labels
  mutate(gene = factor(.[["gene"]], levels=unique(.[["gene"]])))

## Plot E1 + E2 Ubiquitin Ligases expression fold-change
pltE1E2 <- 
  ubiq.AM %>% 
  filter(enzyme != "E3 ubiquitin ligase") %>%
  ggplot(aes(x = gene,
             y = log2FoldChange,
             ymin = log2FoldChange - lfcSE,
             ymax = log2FoldChange + lfcSE,
             fill = complex)) + 
  geom_bar(stat = "identity", ## Create black outline bars
           position = "dodge",
           color = "black",
           size = 0.25,
           show.legend=FALSE) + 
  facet_grid(~ complex, ## Separate the bar plot by complex/enzyme
             scales="free_x", 
             space="free_x") +
  geom_text(aes(label=padj_stars, # draw significance stars
                y = log2FoldChange + sign(log2FoldChange) * (lfcSE + 0.5)), # position significance stars by edges of error bars
            angle = 90, # rotate significance stars
            vjust = "center",
            size = 3,
            nudge_x = 0.1) + ## Add minor adjustment to the stars
  geom_errorbar(width = 0.25, size = 0.25) + ## Add error bars
  xlab("") + ylab(expression(paste(log[2]," fold-change"))) + ## Remove X-label and add Y-label
  scale_fill_viridis_d(option = "C") + ## Change bar colors
  theme_light() + ## Change to white backgroud
  theme(text = element_text(size=7),
        axis.text.x = element_text(angle = 60, # Format x-axis, Angle
                                   hjust = 1, # bottom of axis
                                   face = "italic"), # italics
        strip.text.x = element_text(color = "black", face = "bold"),  # Format Facet title
        strip.background = element_rect(colour=NA,fill=NA), # Remove facet bkrng
        panel.border=element_rect(colour="grey50"))  # Format Facet remove outline

## Plot E3 Ubiquitin Ligases expression fold-change
pltE3 <- 
  ubiq.AM %>% 
  filter(enzyme == "E3 ubiquitin ligase") %>%
  ggplot(aes(x = gene,
             y = log2FoldChange,
             ymin = log2FoldChange - lfcSE,
             ymax = log2FoldChange + lfcSE,
             fill = complex)) + 
  geom_bar(stat = "identity", ## Create black outline bars
           position = "dodge",
           color = "black",
           size = 0.25,
           show.legend=FALSE) + 
  ylim(-3,+3) + 
  facet_grid(~ complex, ## Separate the bar plot by complex/enzyme
             scales="free_x", 
             space="free_x") +
  geom_text(aes(label=padj_stars, # draw significance stars
                y = log2FoldChange + sign(log2FoldChange) * (lfcSE + 0.75)), # position significance stars by edges of error bars
            angle = 90, # rotate significance stars
            vjust = "center",
            size = 3,
            nudge_x = 0.2) + ## Add minor adjustment to the stars
  geom_errorbar(width = 0.25, size = 0.25) + ## Add error bars
  xlab("") + ylab(expression(paste(log[2]," fold-change"))) + ## Remove X-label and add Y-label
  scale_fill_viridis_d() + ## Change bar colors
  theme_light() + ## Change to white backgroud
  theme(text = element_text(size=7),
        axis.text.x = element_text(angle = 60, # Format x-axis, Angle
                                   hjust = 1, # bottom of axis
                                   face = "italic"), # italics
        strip.text.x = element_text(color = "black",angle = 0),  # Format Facet title
        strip.background = element_rect(colour=NA,fill=NA), # Remove facet bkrng
        panel.border=element_rect(colour="grey50"))  # Format Facet remove outline

## Align the two plots vertically, is a 1:2 height ratio.
labels = c("","E3 ubiquitin ligases")

cairo_pdf(filename = paste0(path,today,"-ubiquitin-foldchange-AM140-raw.pdf"),
          width = 6.574803,
          height = 3.54331,
          pointsize = 7)
ggarrange(pltE1E2,pltE3,labels=labels,heights=2:3)
dev.off()
