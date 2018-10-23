source("https://bioconductor.org/biocLite.R")
biocLite("clusterProfiler") ; library(clusterProfiler)
biocLite("org.Ce.eg.db") ; library(org.Ce.eg.db)

## convert biological IDs using clusterProfiler ## Only genes down-regualted in both N2 and AM140
# keytypes(org.Ce.eg.db) # which ids are available to C. elegans
ids <- bitr(both.sets$both.dn, fromType="WORMBASE", toType="ENTREZID", OrgDb="org.Ce.eg.db")

## KEGG over-representation test ##
kk <- enrichKEGG(gene         = ids[['ENTREZID']],
                 keyType      = 'ncbi-geneid',
                 organism     = 'cel', 
                 pvalueCutoff = 0.05)
head(kk, n=15)

## plotting and pathway visualization ##

# barplot(kk,
#         showCategory=12,
#         drop=TRUE,
#         colorBy = 'p.adjust'
#         )
#
#dotplot(kk,showCategory=12)

##Visualize the pathways
# browseKEGG(kk,"cel04120")


## plotting by ggplot2

kk <- as.data.frame(kk)

factor(kk$p.adjust,levels = names(table(kk$p.adjust))[order(kk$p.adjust)])
order(kk$p.adjust)
ggplot(kk, aes(x = reorder(Description,p.adjust), y = Count)) +
  geom_bar(aes(fill = p.adjust),stat = 'identity') + 
  coord_flip() + 
  theme_classic() + 
  scale_fill_gradient(name = "P value",low = "red3",high = "orange") + 
  #ggtitle("KEGG pathways") + 
  ylab("Gene count") + 
  xlab("") + 
  theme(text = element_text(size=7))
ggsave("results\\KEGG_OVERREP_BOTH.DOWN.wmf", height = 4.5, width = 10, units = "cm")

## export data to dable
write.csv(as.data.frame(kk),'results\\KEGG_OVERREP_BOTH.DOWN.csv')

sessionInfo()
# R version 3.5.1 (2018-07-02)
# Platform: i386-w64-mingw32/i386 (32-bit)
# Running under: Windows >= 8 x64 (build 9200)
# 
# Matrix products: default
# 
# locale:
#   [1] LC_COLLATE=English_Israel.1252  LC_CTYPE=English_Israel.1252    LC_MONETARY=English_Israel.1252
# [4] LC_NUMERIC=C                    LC_TIME=English_Israel.1252    
# 
# attached base packages:
#   [1] stats4    parallel  grid      stats     graphics  grDevices utils     datasets  methods  
# [10] base     
# 
# other attached packages:
#   [1] org.Ce.eg.db_3.6.0    AnnotationDbi_1.42.1  Biobase_2.40.0        clusterProfiler_3.8.1
# [5] VennDiagram_1.6.20    futile.logger_1.4.3   bindrcpp_0.2.2        PWMEnrich_4.16.0     
# [9] Biostrings_2.48.0     XVector_0.20.0        IRanges_2.14.10       S4Vectors_0.18.3     
# [13] BiocGenerics_0.26.0   forcats_0.3.0         stringr_1.3.1         dplyr_0.7.6          
# [17] purrr_0.2.5           readr_1.1.1           tidyr_0.8.1           tibble_1.4.2         
# [21] ggplot2_3.0.0         tidyverse_1.2.1       BiocInstaller_1.30.0 
# 
# loaded via a namespace (and not attached):
#   [1] gridExtra_2.3        GOSemSim_2.6.2       stringi_1.1.7        evaluate_0.11       
# [5] memoise_1.1.0        haven_1.1.2          lattice_0.20-35      cli_1.0.0           
# [9] RSQLite_2.1.1        fansi_0.3.0          evd_2.3-3            blob_1.1.1          
# [13] DBI_1.0.0            data.table_1.11.4    rstudioapi_0.7       bindr_0.1.1         
# [17] units_0.6-1          nlme_3.1-137         zlibbioc_1.26.0      qvalue_2.12.0       
# [21] ggrepel_0.8.0        UpSetR_1.3.3         rprojroot_1.3-2      tools_3.5.1         
# [25] magrittr_1.5         Rcpp_0.12.18         xml2_1.2.0           gdata_2.18.0        
# [29] readxl_1.1.0         httr_1.3.1           rmarkdown_1.10       assertthat_0.2.0    
# [33] R6_2.2.2             fastmatch_1.1-0      munsell_0.5.0        cellranger_1.1.0    
# [37] gtools_3.8.1         digest_0.6.16        splines_3.5.1        cowplot_0.9.3       
# [41] DOSE_3.6.1           colorspace_1.3-2     ggraph_1.0.2         pkgconfig_2.0.2     
# [45] pillar_1.3.0         GO.db_3.6.0          formatR_1.5          tweenr_1.0.0        
# [49] rvcheck_0.1.1        plyr_1.8.4           farver_1.0           gtable_0.2.0        
# [53] tidyselect_0.2.4     rvest_0.3.2          reshape2_1.4.3       knitr_1.20          
# [57] viridisLite_0.3.0    futile.options_1.0.1 rlang_0.2.2          broom_0.5.0         
# [61] glue_1.3.0           Matrix_1.2-14        backports_1.1.2      scales_1.0.0        
# [65] ggridges_0.5.1       lubridate_1.7.4      lambda.r_1.2.3       modelr_0.1.2        
# [69] igraph_1.2.2         enrichplot_1.0.2     ggforce_0.1.3        fgsea_1.6.0         
# [73] hms_0.4.2            labeling_0.3         htmltools_0.3.6      yaml_2.2.0          
# [77] lazyeval_0.2.1       utf8_1.1.4           crayon_1.3.4         withr_2.1.2         
# [81] MASS_7.3-50          BiocParallel_1.14.1  bit64_0.9-7          seqLogo_1.46.0      
# [85] DO.db_2.9            viridis_0.5.1        jsonlite_1.5         compiler_3.5.1      
# [89] bit_1.1-14    