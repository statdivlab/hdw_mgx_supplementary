# Figure 3 
library(pheatmap)
library(RColorBrewer)
card_heatmap <- readRDS("card_heatmap.RDS")
card_heatmap <- card_heatmap[ , c("C1","C2","C3","C4","C5","C6","D1","D2","D3","D4","D5","D6","D7","D8","D9","D10")]
pheatmap(card_heatmap, color = colorRampPalette(brewer.pal(9,"Blues")[3:9])(100), treeheight_row = 0, treeheight_col = 0, cluster_cols = F, cluster_rows = F, 
         cellwidth = 15, cellheight = 12, fontsize = 8, filename = "figure3.pdf",  na_col = "white")
