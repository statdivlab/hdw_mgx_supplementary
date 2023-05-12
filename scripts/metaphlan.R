# script running CLR tests at the species and phylum level 
library(data.table)

clr <- function(x) {
  (x %>% log) - (x %>% log %>% mean)
}

all_long_counts_phylum <- readRDS("all_long_counts_phylum.RDS")

clr_tests_counts_phylum <- all_long_counts_phylum %>%
  group_by(sample) %>% 
  mutate("clr" = clr(abundance_counts)) %>%
  group_by(clade_name, group) %>%
  summarise(mean = mean(clr), sd = sd(clr), n = n()) %>%
  filter(n > 2) %>% 
  pivot_wider(names_from=group, values_from=mean:n) %>%
  filter(n_dairy > 2, n_community > 2) %>%
  mutate(t = (mean_community - mean_dairy)/sqrt(sd_community^2/n_community + sd_dairy^2/n_dairy),
         df = n_community + n_dairy - 2, 
         pvalue = 2*pt(-abs(t),df=df)) %>%
  arrange(abs(pvalue)) 

clr_tests_counts_phylum$qvalue <- qvalue::qvalue(clr_tests_counts_phylum$pvalue, pi0.method = "bootstrap", lambda = 0.5)$qvalues
write.csv(clr_tests_counts_phylum,"STable2.csv")


all_long_counts <- readRDS("all_long_counts.RDS")

clr_tests_counts <- all_long_counts %>%
  group_by(sample) %>% 
  mutate("clr" = clr(abundance_counts)) %>%
  group_by(clade_name, group) %>%
  summarise(mean = mean(clr), sd = sd(clr), n = n()) %>%
  filter(n > 2) %>% 
  pivot_wider(names_from=group, values_from=mean:n) %>%
  filter(n_dairy > 2, n_community > 2) %>%
  mutate(t = (mean_community - mean_dairy)/sqrt(sd_community^2/n_community + sd_dairy^2/n_dairy),
         df = n_community + n_dairy - 2, 
         pvalue = 2*pt(-abs(t),df=df)) %>%
  arrange(abs(pvalue)) 

clr_tests_counts$qvalue <- qvalue::qvalue(clr_tests_counts$pvalue)$qvalues
write.csv(clr_tests_counts,"STable1.csv")


# Figure S1 
pcoa_meta <- readRDS("pcoa_meta.RDS")
png("pcoa_metaphlan3.png", units="in", width=6, height=4, res=300)  
ggplot(pcoa_meta, aes(x = PC1, y = PC2, color = group)) + 
  geom_point(size = 5) + 
  ggtitle("PCoA Using Bray-Curtis Distances") + 
  theme_bw() + 
  scale_color_manual(values = c("dairy" = "darkorange", "community" = "darkblue"))
dev.off()


#### Metaphlan Comparisons with and without HMP cohort: Figure S2

combined_model_results <- readRDS("combined_model_results.RDS")
compare_betas <- ggplot(combined_model_results, aes(y = dairy_beta_model1, x = dairy_beta)) + 
  geom_point(size = 2, alpha = 0.8) + 
  theme_bw() + 
  xlab(expression(paste(beta[1], " dairy estimate (w/ HMP)"))) + 
  ylab(expression(paste(beta[1], " dairy estimate (w/o HMP)"))) + 
  geom_abline(intercept = 0, slope = 1) + 
  ggrepel::geom_text_repel(aes(label = species), size = 2) + 
  scale_x_continuous(limits = c(-7.5,6)) + 
  scale_y_continuous(limits = c(-7.5,6))

compare_pvalues<- ggplot(combined_model_results, aes(y = log10_p_model1, x = log10_p_model2)) + 
  geom_point(size = 2, alpha = 0.8) + 
  theme_bw() + 
  ylab("-log10 unadjusted p-value (w/o HMP cohort)") + 
  xlab("-log10 unadjusted p-value (w/ HMP cohort)") + 
  geom_abline(intercept = 0, slope = 1) + 
  ggrepel::geom_text_repel(aes(label = species), size = 2) + 
  scale_x_continuous(limits = c(0,16)) + 
  scale_y_continuous(limits = c(0,16))

library(ggpubr)
pdf("HMP_metaphlan_comparison.pdf", width = 12, height = 6)
ggarrange(compare_betas, compare_pvalues, 
          nrow = 1)
dev.off()

# Figure S3
heatmap_species_clr <- readRDS("heatmap_species_clr.RDS")
pheatmap(heatmap_species_clr,color = colorRampPalette(brewer.pal(9,"Blues")[3:9])(100), treeheight_row = 0, treeheight_col = 0, cluster_cols = F, cluster_rows = F, 
         cellwidth = 15, cellheight = 12, fontsize = 8, filename = "FigureS3.pdf",  na_col = "white")

