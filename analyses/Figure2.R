# Creation of Figure 2
Figure2 <- readRDS("Figure2_data.RDS")

card_fig <- ggplot(Figure2, aes(x = scaledreads, y = n_card, color = group)) + 
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

vfdb_fig <- ggplot(Figure2, aes(x = scaledreads, y = n_vfdb, color = group)) + 
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


png("figure2.png", units = "in", width = 12, height = 6, res = 300)
ggarrange(ggarrange(card_fig, vfdb_fig, ncol = 2,  common.legend=TRUE, 
                    legend = "bottom"), 
          labels = c("",""), nrow = 1)
dev.off()
