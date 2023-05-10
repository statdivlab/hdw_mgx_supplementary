# Comparison of ARGs and VFDBs with HMP healthy human subjects cohort of 239 samples out of 300 participants
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

card_files <- list.files(path = "/Users/paulinetrinh/Documents/HDW/metagenomic/HMP_comparison/abricate_results", pattern = "*card.tab", full.names = T)
vfdb_files <- list.files(path = "/Users/paulinetrinh/Documents/HDW/metagenomic/HMP_comparison/abricate_results", pattern = "*vfdb.tab", full.names = T)

tbl <- sapply(card_files, read_tsv, simplify=FALSE) 
card_tbl <- do.call(rbind,tbl)

colnames(card_tbl) # looks good! 

split_tbl <- mutate(card_tbl,'DNA Accession'= unlist(lapply(strsplit(ACCESSION, ':',fixed = TRUE), '[', 1)))

# let's look at number of contigs 




#AMR <- merge(split_tbl, card, all = FALSE)

#n_occur <- data.frame(table(AMR$`DNA Accession`))
#n_occur[n_occur$Freq > 1,]

## VFDB files 
tbl <- sapply(card_files, read_tsv, simplify=FALSE) %>% 
  bind_rows(.id = "ID")
colnames(tbl) # looks good! 

split_tbl <- mutate(tbl,'DNA Accession'= unlist(lapply(strsplit(ACCESSION, ':',fixed = TRUE), '[', 1)))

AMR <- merge(split_tbl, card, all = FALSE)

n_occur <- data.frame(table(AMR$`DNA Accession`))
n_occur[n_occur$Freq > 1,]