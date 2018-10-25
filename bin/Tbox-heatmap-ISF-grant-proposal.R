row_dend = hclust(dist(pval.mat.sig)) # row clustering
col_dend = hclust(dist(t(pval.mat.sig))) # column clustering

pdf("tbx_heatmap_up_grant.pdf",width=4,height=10)
Heatmap(pval.mat.sig, 
        name = "Motif enrichment p-value", #title of legend
        column_title = "T-BOX motif enrichment", 
        row_title = "Genes upregulated by anc-1 RNAi",
        column_names_side = "top",
        row_names_gp = gpar(fontsize = 7), # Text size for row names
        column_names_gp = gpar(fontsize = 7), # Text size for column names
        cluster_rows = color_branches(row_dend, k = 2),
        #km = 8 ,# To split the dendrogram using k-means, divide into 4 groups. Must set.seed
        cluster_columns = color_branches(col_dend, k = 3),
        col = mycols)

dev.off()

dev.print(png, 'filename.png')
dev.copy2pdf('filename.pdf',out.type = "pdf")
dev.off()