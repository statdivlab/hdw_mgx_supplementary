# gene CLR tests 
library(tidyverse)
library(pheatmap)
library(RColorBrewer)

# Code for Supplementary Figure S4
#card_heatmap_all <- readRDS("card_heatmap_all.RDS")
pheatmap(card_heatmap_all_plot, color = colorRampPalette(brewer.pal(9,"Blues")[3:9])(100), treeheight_row = 0, treeheight_col = 0, cluster_cols = F, cluster_rows = F, 
         cellwidth = 15, cellheight = 12, fontsize = 8, filename = "card_heatmap_log10_all.pdf",  na_col = "white")

# CLR t-test results for AMR and VF genes

clr <- function(x) {
  (x %>% log) - (x %>% log %>% mean)
}

#clr_genes_p1 <- readRDS("clr_genes_p1.RDS")
clr_genes_p1_long <- clr_genes_p1 %>% 
  rownames_to_column("sample") %>%
  pivot_longer(!sample, names_to = "genes", values_to = "gene_coverage") %>% 
  ungroup() %>% 
  group_by(sample) %>% 
  mutate(total_genecoverage = sum(gene_coverage)) %>% 
  mutate(relativeabundance = gene_coverage/total_genecoverage,
         "clr" = clr(relativeabundance)) %>% 
  ungroup() 
#clr_genes_p1_long <- readRDS("clr_genes_p1_long.RDS")
#genes_source_ref <- readRDS("genes_source_ref.RDS")
clr_results_genes_p1 <-  clr_genes_p1_long %>% 
  ungroup() %>% 
  group_by(sample) %>% 
  mutate(clr = clr(relativeabundance), 
         group = ifelse(sample %in% c("C02_S104","C02_S106",
                                      "C02_S111","C02_S112",
                                      "C02_S113","C02_S115"),"community","dairy")) %>% 
  ungroup() %>% 
  group_by(genes, group) %>%
  summarise(mean = mean(clr), sd = sd(clr), n = n()) %>%
  pivot_wider(names_from=group, values_from=mean:n) %>% 
  mutate(t = (mean_community - mean_dairy)/sqrt(sd_community^2/n_community + sd_dairy^2/n_dairy), 
         df = n_community + n_dairy - 2,
         pvalue = 2*pt(-abs(t),df=df), 
         mean_diff = mean_community - mean_dairy) %>%
  arrange(desc(abs(t))) %>% 
  left_join(genes_source_ref)

clr_results_genes_p1$qvalue <- qvalue::qvalue(clr_results_genes_p1$pvalue, pi0.method = "bootstrap", lambda = seq(0,0.95,0.05))$qvalue
# clr_results_genes_p1 <- readRDS("clr_results_genes_p1.RDS")


# Antibiotic classes CLR tests 
#clr_abxclass_p1 <- readRDS("clr_abxclass_p1.RDS")
clr_abxclass_p1_long <- clr_abxclass_p1 %>% 
  rownames_to_column("sample") %>%
  pivot_longer(!sample, names_to = "abxclass", values_to = "gene_coverage") %>% 
  ungroup() %>% 
  group_by(sample) %>% 
  mutate(total_genecoverage = sum(gene_coverage)) %>% 
  mutate(relativeabundance = gene_coverage/total_genecoverage,
         "clr" = clr(relativeabundance)) %>% 
  ungroup()
clr_results_abxclass_p1 <-  clr_abxclass_p1_long %>% 
  ungroup() %>% 
  group_by(sample) %>% 
  mutate(clr = clr(gene_coverage), 
         group = ifelse(sample %in% c("C02_S104","C02_S106",
                                      "C02_S111","C02_S112",
                                      "C02_S113","C02_S115"),"community","dairy")) %>% 
  filter(abxclass %in% c("aminoglycoside","fluoroquinolone","macrolide",
                         "tetracycline","cephamycin","glycopeptide","sulfonamide","cephalosporin")) %>% 
  ungroup() %>% 
  group_by(abxclass, group) %>%
  summarise(mean = mean(clr), sd = sd(clr), n = n()) %>%
  pivot_wider(names_from=group, values_from=mean:n) %>% 
  mutate(t = (mean_community - mean_dairy)/sqrt(sd_community^2/n_community + sd_dairy^2/n_dairy), 
         df = n_community + n_dairy - 2,
         pvalue = 2*pt(-abs(t),df=df), 
         mean_diff = mean_dairy - mean_community) %>%
  arrange(desc(abs(t))) 

clr_results_abxclass_p1$qvalue <- qvalue::qvalue(clr_results_abxclass_p1$pvalue, pi0.method = "bootstrap", lambda = seq(0,0.9,0.05))$qvalue
#saveRDS(clr_results_abxclass_p1, "clr_results_abxclass_p1.RDS")
