#### 
# Figures 1, 2 & 3 Update 
#################
#### Figure 1 
#################
library(Polychrome)
library(ggpubr)
phylum_plot <- readRDS("data/phylum_plot.RDS")
species_plot <- readRDS("data/species_plot.RDS")

set.seed(1302)
P48 = createPalette(48,  c("#ff0000", "#00ff00", "#0000ff"), range = c(30, 80), M = 10000)
P48 <- sortByHue(P48)
P48 <- as.vector(t(matrix(P48, ncol=4)))
names(P48) <- NULL

species_RA_plot <- ggplot(species_plot, aes(y=abundance, x=sample_names, fill=plot_species)) +
  geom_bar(stat="identity") + 
  theme(axis.text.x = element_text(angle=90)) + 
  scale_fill_manual(values = P48)+
  guides(fill=guide_legend(nrow=24, title = "Species", 
                           title.position = "top",title.hjust = 0.5)) +
  theme_bw() +
  theme(legend.position="bottom", 
        plot.title = element_text(size = 12, face = "bold"),
        legend.title=element_text(size=27), 
        legend.text=element_text(size=27)) + 
  xlab("") + 
  ylab("relative abundance") + 
  theme(axis.text.x=element_text(size=35), axis.text.y=element_text(size=35), axis.title=element_text(size=35))


custom_palette <- c("#661100","#88CCEE","#CC6677", "#DDCC77","#117733",
                    "#332288", "#AA4499", "#44AA99", "#882255")
phylum_RA_plot <- ggplot(phylum_plot, aes(y=abundance, x=sample_names, fill=phylum)) +
  geom_bar(stat="identity") + 
  theme(axis.text.x = element_text(angle=90)) + 
  scale_fill_manual(values = custom_palette) +
  guides(fill=guide_legend(nrow=2, title = "Phyla")) +
  theme_bw() + 
  theme(legend.position="bottom", 
        plot.title = element_text(size = 12, face = "bold"),
        legend.title=element_text(size=27), 
        legend.text=element_text(size=27)) + 
  xlab("") + 
  ylab("relative abundance") + 
  theme(axis.text.x=element_text(size=35), axis.text.y=element_text(size=35), axis.title=element_text(size=35))


pdf("figure1.pdf", width = 30, height = 25)
ggarrange(phylum_RA_plot, species_RA_plot, 
          labels = c("A","B"),
          font.label=list(color="black",size=40),
          nrow = 1)
dev.off()


###################
# Figure 2 
###################
Figure2 <- readRDS("data/Figure2_Dec2.RDS") %>% 
  mutate(figure_names = ifelse(group == "HMP", sample_id, sample_names))

Figure2_rev <- readRDS("data/hmp_contigs_data.RDS") %>% 
  rename(figure_names = sample_id) %>% left_join(Figure2) %>% 
  mutate(scaledreads = `total pairs after host removal`/1000000,
         dairy = ifelse(group == "dairy", "dairy", 
                        ifelse(group == "community", "community", NA)))


## Stats for text 
Figure2_rev %>% 
  group_by(group) %>% 
  summarise(mean_card = mean(n_card), 
            sd_card =  sd(n_card), 
            mean_vf = mean(n_vfdb), 
            sd_vf = sd(n_vfdb), 
            mean_sequencereads = mean(`total pairs after host removal`), 
            sd_reads = sd(`total pairs after host removal`))  

Figure2_rev %>% 
  group_by(group) %>% 
  summarise(min_vf = min(n_vfdb), 
            max_vf = max(n_vfdb),
            min_card = min(n_card), 
            max_card = max(n_card))

t.test(n_vfdb ~ dairy, data = Figure2_rev) # p = 0.08
t.test(n_card ~ dairy, data = Figure2_rev) # p = 0.09



card_fig <- ggplot(Figure2_rev, aes(x = scaledreads, y = n_card, color = group)) + 
  geom_point(size = 5, alpha = 0.85) + 
  xlim(c(0,250)) + ylim(c(0,60)) + 
  theme_bw() + 
  labs(x = expression(paste("Number of Reads Sequenced ", (10^6))), y = "Number of CARD genes") + 
  scale_color_manual(values = c("HMP" = "cornflowerblue","dairy" = "darkorange", "community" = "darkblue")) + 
  theme(legend.position="bottom") + 
  theme(axis.text.x=element_text(size=16), 
        axis.text.y=element_text(size=16), 
        axis.title=element_text(size=16), 
        legend.text=element_text(size=16), 
        legend.title = element_text(size=16)) + 
  ggrepel::geom_text_repel(aes(label = age), size = 5)

vfdb_fig <- ggplot(Figure2_rev, aes(x = scaledreads, y = n_vfdb, color = group)) + 
  geom_point(size = 5, alpha = 0.85) + 
  xlim(c(0,250)) + ylim(c(0,60)) + 
  theme_bw() + 
  labs(x = expression(paste("Number of Reads Sequenced ", (10^6))), y = "Number of VFDB genes") + 
  scale_color_manual(values = c("dairy" = "darkorange", "community" = "darkblue", "HMP" = "cornflowerblue")) + 
  theme(legend.position="bottom") + 
  theme(axis.text.x=element_text(size=16), 
        axis.text.y=element_text(size=16), 
        axis.title=element_text(size=16), 
        legend.text=element_text(size=16), 
        legend.title = element_text(size=16)) + 
  ggrepel::geom_text_repel(aes(label = age), size = 5)


png("Figure2.png", units = "in", width = 12, height = 6, res = 300)
ggpubr::ggarrange(vfdb_fig, card_fig, 
                            labels = c("A","B"),
                            ncol = 2,  nrow = 1,
                            common.legend=TRUE, 
                    legend = "bottom")
dev.off()


##############
## Figure 3
##############
# Figure 3 
library(pheatmap)
library(RColorBrewer)
old_card_heatmap <- readRDS("data/card_heatmap.RDS")
new_qvalues <- readRDS("data/2024_amr_pvalue_results.RDS") %>% 
  rename(genes = species) %>% 
  select(genes, amr_qvalue_dairy)
card_heatmap <- old_card_heatmap %>% as.data.frame() %>% 
  rownames_to_column("names") %>% 
  separate(names, c("genes","old_qval"), sep = " ") %>% 
  full_join(new_qvalues) %>% 
  mutate(abx_pvalue = paste0(genes," (q=",round(amr_qvalue_dairy,2),")")) %>%
  select(-c(genes, old_qval, amr_qvalue_dairy)) %>% 
  column_to_rownames("abx_pvalue")
  

card_heatmap <- card_heatmap[ , c("C1","C2","C3","C4","C5","C6","D1","D2","D3","D4","D5","D6","D7","D8","D9","D10")]
pheatmap(card_heatmap, color = colorRampPalette(brewer.pal(9,"Blues")[3:9])(100), treeheight_row = 0, treeheight_col = 0, cluster_cols = F, cluster_rows = F, 
         cellwidth = 15, cellheight = 12, fontsize = 8, filename = "Figure3.pdf",  na_col = "white")




