columns.names <- c("Gene_stable_ID","Gene_name", "strain")
hyper <- read_tsv("hypersensitive to protein aggregation induced paralysis via RNAi.tsv", col_names = F)
names(hyper) <- columns.names
resist <- read_tsv("resistant to protein aggregation induced paralysis via RNAi.tsv", col_names = F)
names(resist) <- columns.names

hyper.DEG <- semi_join(both,hyper,by = "Gene_stable_ID")
hyper.DEG %>% ggplot(aes(x = Gene_name, y = log2FoldChange, fill = padj)) + 
  geom_bar(stat = "identity") + 
  facet_grid(strain ~ .)
View(hyper.DEG)

resist.DEG <- semi_join(both,resist,by = "Gene_stable_ID")
resist.DEG %>% ggplot(aes(x = Gene_name, y = log2FoldChange, fill = padj)) + 
  geom_bar(stat = "identity") + 
  facet_grid(strain ~ .)
View(resist.DEG)
