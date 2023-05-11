# Figure 1 
library(Polychrome)
library(ggpubr)
phylum_plot <- readRDS("phylum_plot.RDS")
species_plot <- readRDS("species_plot.RDS")

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
          nrow = 1)
dev.off()


