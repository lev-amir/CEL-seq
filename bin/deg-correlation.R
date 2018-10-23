## to export plots with transperency (alpha value) not as PDF
# install.packages("Cairo") ; library(Cairo)

## combine the two strain data-frames, removing genes that don't appear in both.
bothStrains <- inner_join(AM,NT, by = 'Gene_stable_ID') 
# add variable to indicate which gene are statistically significant in both strains
bothStrains <- bothStrains %>%
  mutate(significant_change = padj.x < 0.05 & padj.y < 0.05)
# keep only genes which are statistically significant in both strains, for correlation calc.
bothStrain.signif <- bothStrains %>% filter(significant_change == TRUE)

glimpse(bothStrains)

## correlations between the two strains

corBothStrains <- cor(bothStrains$log2FoldChange.x,bothStrains$log2FoldChange.y)
#cor.test(bothStrains$log2FoldChange.x,bothStrains$log2FoldChange.y)
corBothStrain.signif <- cor(bothStrain.signif$log2FoldChange.x,bothStrain.signif$log2FoldChange.y)
#cor.test(bothStrain.signif$log2FoldChange.x,bothStrain.signif$log2FoldChange.y)

## Scatter plot gene expression fold-change of the two strains

bothStrains %>%
  ggplot(aes(x = log2FoldChange.x,
             y = log2FoldChange.y,
             col = significant_change,
             alpha = significant_change)) +
  geom_point(shape = 16,size = 0.1) + 
  annotate(geom = "text", 
           x = -3, 
           y = 2.5, 
           label = paste("r =",signif(corBothStrains)),
           size = 1) +
  annotate(geom = "text", 
           x = -3, 
           y = 2.2, 
           label = paste("r_significant =",signif(corBothStrain.signif)),
           size = 1) +
  theme_classic() +
  ylab("N2 Fold-Change (log2)") +
  xlab("AM140 Fold-Change (log2)") + 
  scale_alpha_manual(values = c(0.1,0.4), na.value = 0.1) + 
  scale_color_manual(values = c("Black","red"),na.value = "Black") + 
  theme(legend.position = "none") + 
  theme(text = element_text(size=7))
ggsave("results\\strains correlation\\strains_correlation_scatterPlot.pdf", 
       height = 4.5, 
       width = 4.5, 
       units = "cm")
