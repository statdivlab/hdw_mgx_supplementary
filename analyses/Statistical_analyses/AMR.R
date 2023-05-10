# let's look at the AMR gene search results with the CARD database first 
# Go to CARD and download the data 
# https://card.mcmaster.ca/download
# Download CARD Data and unzip 
# Use the aro_categories_index.tsv file for drug classes and DNA accessions 
install.packages("ggpattern")
library(tidyverse)
card <- read_tsv(file = '/Users/paulinetrinh/Documents/HDW/metagenomic/card-data/aro_categories_index.tsv')
dim(card)
colnames(card)
# checking for duplication accessions 
n_occur <- data.frame(table(card$`DNA Accession`))
n_occur[n_occur$Freq > 1,]
# well it looks like there are a lot of hits which makes sense why merging resulted in more observations than expected... 
# So what do I do about this? there are different antibiotic classes linking to each accession number 
# what is the appropriate level of granularity here 


files <- list.files(path = "/Users/paulinetrinh/Documents/HDW/metagenomic/contigs-amr-analysis/results", pattern = "*.tab", full.names = T)

tbl <- sapply(files, read_tsv, simplify=FALSE) %>% 
  bind_rows(.id = "ID")
colnames(tbl) # looks good! 

split_tbl <- mutate(tbl,'DNA Accession'= unlist(lapply(strsplit(ACCESSION, ':',fixed = TRUE), '[', 1)))

AMR <- merge(split_tbl, card, all = FALSE)

n_occur <- data.frame(table(AMR$`DNA Accession`))
n_occur[n_occur$Freq > 1,]

### Read in the full summary across all databases 
all_amr_results <- read_tsv(file ='/Users/paulinetrinh/Documents/HDW/metagenomic/contigs-amr-analysis/amr_results_summary.txt')
# Let's take out any columns that are all empty 
all_amr_results %>% 
  separate("#FILE", sep = "-", c("ID","intermediate")) %>% 
  separate("intermediate", sep = ".tab", c("database","dropme")) %>% 
  select(-c(dropme,NUM_FOUND)) -> all_amr_results

all_amr_results %>% 
  pivot_longer(!c(ID,database), names_to = "gene_name", values_to = "gene_coverage") %>% 
  filter(!gene_coverage == ".") %>% 
  separate("gene_coverage", sep = ";", c("gene1","gene2","gene3")) %>% 
  mutate(gene1 = as.numeric(as.character(gene1))*1, 
         gene2 = as.numeric(as.character(gene2))*1, 
         gene3 = as.numeric(as.character(gene3))*1, 
         thresholdmin =  ifelse(gene1 == 100, 1, 
                                ifelse(gene2 == 100, 1, 
                                       ifelse(gene3 == 100, 1, 0)))) -> amr_long

amr_long %>% 
  filter(thresholdmin == 1) -> complete_genes

complete_genes %>% 
  group_by() %>% 
  
complete_genes %>% 
  distinct(gene_name) -> distinct_genes

write_csv(distinct_genes,"distinct_genes.csv")

# So there are 1384 genes identified 



##### Okay let's read in the summary of how many of each ARGs were found in CARD and resfinder per sample? 
summary_amr <- read_tsv(file = '/Users/paulinetrinh/Documents/HDW/metagenomic/summary.txt')
card <- read_tsv(file = '/Users/paulinetrinh/Documents/HDW/metagenomic/card_summary.txt')
resfinder <- read_tsv(file = '/Users/paulinetrinh/Documents/HDW/metagenomic/resfinder_summary.txt')

summary_amr %>% 
    separate(FILE, c("person","source"), sep = "-") %>% 
    separate(source, c("database","tab"), sep = ".tab") %>%
    separate(person, c("group","ID")) %>% 
    select(-tab) %>% 
    mutate(xpos = rep(seq(1,16), each = 2)) -> AMR
View(AMR)

# Create data.frame with shading info
shading <- data.frame(min = seq(from = 0.5, to = max(as.numeric(as.factor(AMR$group))), by = 1),
                      max = seq(from = 1.5, to = max(as.numeric(as.factor(AMR$group))) + 0.5, by = 1),
                      col = c(0,1))
# create barplot of AMR frequency per sample
png(filename = "/Users/paulinetrinh/Documents/DissertationProposal/AMR_barplot.png", width = 8, height = 6, units = 'in', res = 300)
AMR %>%
  ggplot(aes(x=xpos, y=NUM_FOUND, fill=database)) +
  geom_bar(stat="identity", color = "black", position=position_dodge()) +
  geom_rect(aes(
    xmin = xpos - .5,
    xmax = xpos + .5,
    ymin = -Inf,
    ymax = Inf,
    fill = ID), alpha = .2) +
  geom_col(
    position = "dodge") +
  scale_x_continuous(
    breaks = AMR$xpos %>% unique(),
    labels = AMR$ID %>% unique(),
    expand = c(0, 0)
  )  + 
  scale_fill_manual(
    values = c(
      'lightcoral',  
      'cornflowerblue',  
      'white',  
      'white',  
      'white',  
      'white', 
      'white',
      'white',
      'grey',
      'grey',
      'grey',
      'grey',
      'grey',
      'grey',
      'grey',
      'grey',
      'grey',
      'grey'
    ),
    breaks = c(
      'card',
      'resfinder'
    )
  ) + 
  theme_bw() + 
  ylab("# of ARGs") + 
  xlab("Sample ID") + 
  theme(legend.position='bottom', legend.title= element_blank())  + 
  theme(axis.text.x = element_text(angle = 90)) 
dev.off()


AMR %>%
  filter(group == "C02") %>%
  ggplot(aes(x=ID, y=NUM_FOUND, fill=database)) +
  geom_bar(stat="identity", position=position_dodge(), width = 0.8) +
  scale_fill_manual(values = c("lightcoral","cornflowerblue")) +  
  theme_bw() + 
  ylab("# of ARGs") + 
  xlab("Sample ID") + 
  theme(legend.position='bottom', legend.title= element_blank())  + 
  theme(axis.text.x = element_text(angle = 90)) -> community

AMR %>%
  filter(group == "D03") %>%
  ggplot(aes(x=ID, y=NUM_FOUND, fill=database)) +
  geom_bar(stat="identity", position=position_dodge(), width = 0.8) +
  scale_fill_manual(values = c("lightcoral","cornflowerblue")) + 
  theme_bw() + 
  ylab("# of ARGs") + 
  xlab("Sample ID") + 
  theme(legend.position='bottom', legend.title= element_blank())  + 
  theme(axis.text.x = element_text(angle = 90)) -> dairy
write.csv(dairy_AMR, "dairy_resistancegenes.csv")
png("ARG_Figure.png", units="in", width=8, height=8, res=300)
ggpubr::ggarrange(community, dairy, ncol = 2, nrow = 1, common.legend = TRUE, legend = "bottom")
dev.off()


AMR %>% 
  filter(database=="card") %>% 
  mutate(worker = ifelse(group=="D03","dairy","community")) %>% 
  ggplot(aes(x = worker, y = NUM_FOUND, color = worker)) + 
  geom_jitter(size = 4, width = 0.12) + 
  scale_color_manual(values=c("darkblue", "orange")) +
  theme_bw() + 
  ylab("# of ARGs") + 
  xlab("") +
  theme(legend.position="none", legend.title= element_blank()) + 
  ggtitle("CARD") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) -> CARD

AMR %>% 
  filter(database=="resfinder") %>% 
  mutate(worker = ifelse(group=="D03","dairy","community")) %>% 
  ggplot(aes(x = worker, y = NUM_FOUND, color = worker)) + 
  geom_jitter(size = 4, width = 0.12) + 
  scale_color_manual(values=c("darkblue", "orange")) +
  theme_bw() + 
  ylab("# of ARGs") + 
  xlab("") +
  theme(legend.position="none", legend.title= element_blank()) + 
  ggtitle("ResFinder") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) -> resfinder

png("ARG_Figure_Presentation_rescale.png", units="in", width=8, height=6, res=300)
ggpubr::ggarrange(CARD, resfinder, ncol = 2, nrow = 1)
dev.off()
