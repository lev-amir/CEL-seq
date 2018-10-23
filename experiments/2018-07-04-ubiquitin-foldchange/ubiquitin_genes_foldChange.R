if (!require(tidyverse)) {install.packages("tidyverse"); library(tidyverse)}
if (!require(viridis)) {install.packages("viridis"); library(viridis)}

## genes of interest, related to protein ubiquitination
ubiq <- read_tsv("ubiquitin-foldchange\\KEGG_Celegans_UbiqMediateProteolysis.tab")
ubiq$complex[is.na(ubiq$complex)] <- "" # replace NA values w/ empty string
ubiq <- ubiq %>% mutate(description = paste(ubiq$complex,ubiq$enzyme))

## plot expression fold-change
NT %>%
  filter(Gene_name %in% ubiq[["gene"]] | Transcript_stable_ID %in% ubiq[["gene"]]) %>% # select genes of interest
  left_join(ubiq, by = c("Gene_name" = "gene")) %>% # combine columns of two df  
  arrange(description,log2FoldChange) %>% # sort by description first, then by log2FC
  mutate(Gene_name = factor(.[["Gene_name"]], levels=unique(.[["Gene_name"]]))) %>% # needed for plot
  ggplot(aes(x = Gene_name,
             y = log2FoldChange,
             ymin = log2FoldChange - lfcSE,
             ymax = log2FoldChange + lfcSE,
             fill = description)) + 
  geom_bar(position="dodge",
           stat = "identity", 
           color = "black", 
           size = 0.25,
           show.legend=FALSE) + 
  geom_errorbar(width = 0.25, size = 0.25) +
  facet_grid(. ~ description, 
             scales="free_x", 
             space="free_x") + 
  # facet_grid(strain ~ .) +  ## to plot both strains onto the same graph
  geom_text(aes(label=padj_stars, # draw significance stars
            y = log2FoldChange + sign(log2FoldChange) * (lfcSE + 1)), # position significance stars by edges of error bars
            angle = 90, # rotate significance stars
            vjust = "center",
            size = 3,
            nudge_x = .4) +
  # theme(legend.key.size = unit(x = 0.25,units = "cm")) +  # Format legend
  xlab("") + 
  ylab(expression(paste(log[2]," fold-change"))) + 
  scale_fill_viridis_d(option = "C") + 
  theme_light() + 
  theme(text = element_text(size=7),
        axis.text.x = element_text(angle = 60, # Format x-axis
                                   hjust = 1,
                                   face = "italic"),
        strip.text.x = element_text(color = "black"),  # Format Facet
        strip.background = element_rect(colour=NA,fill=NA),
        panel.border=element_rect(colour="grey50"))  # Format Facet
  

## Warning message: Removed 3 rows containing missing values (geom_text). 
## this is likely due to the analysis in which low intensity reads had their 'padj' set as NA.
## in new analysis this should not occur.

## dots size may be too large and need to be reduced to 0.1 in the above `geom_point(shape = 16,size = 2)`
ggsave("ubiquitin-foldchange\\results\\ubiquitin-foldchange-bars.pdf", 
       height = 5, 
       width = 16.7, 
       units = "cm")

ggsave("ubiquitin-foldchange\\results\\ubiquitin-foldchange-bars.svg", 
       height = 5, 
       width = 16.7, 
       units = "cm")

# ggsave("ubiquitin-foldchange\\results\\ubiquitin-foldchange-bars.png",
#        device = "png",
#        height = 5, 
#        width = 16.7, 
#        units = "cm",
#        dpi = 300)

