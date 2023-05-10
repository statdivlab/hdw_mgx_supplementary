# Amy Willis
# Aug 17 2022
# Do a quick and dirty estimate of the total number of species (strains?) in these metagenomes

# Note that this is a kind of new approach so may be worth writing up into a mini paper

library(tidyverse)
library(magrittr)
library(breakaway)
dataframe <- readRDS("HDW_SCG_presence.RDS") %>% 
  as_tibble 



dataframe$present %>% unique

dataframe %<>%   
  filter(present == 1) %>%
  select(-present)
dataframe
dataframe %>% 
  select(1:3, 9:12)


dataframe %>% nrow
dataframe %>% distinct %>% nrow


### Calculate the number of times each species is observed wrt to each gene

dataframe %>%
  filter(sample == "C02_S104") %>%
  group_by("gene_name", "t_domain", "t_phylum", "t_class", "t_order", "t_family", "t_genus", 
             "t_species") %>%
  summarise(n())

frequency_counts <- dataframe %>%
  # distinct %>% # TODO confirm?!?!?!?!?! Not sure if this is the right thing to do. species vs strain diversity? Could also be a good use of breakaway_nof1
  mutate(taxon = paste(t_domain, t_phylum, t_class, t_order, t_family, t_genus, 
                       t_species)) %>%
  select(-c("t_domain", "t_phylum", "t_class", "t_order", "t_family", "t_genus", 
            "t_species")) %>%
  group_by(taxon, sample) %>%
  summarise(n = n()) %>%
  ungroup


frequency_counts %>%
  filter(sample == "C02_S104") %>%
  select(-taxon) %>%
  pull(n) %>%
  table %>%
  as_tibble %>% 
  rename(count = 1, f = 2) %>%
  mutate(count = as.numeric(count)) %>% 
  breakaway::breakaway(.)

bas <- frequency_counts %>%
  split(.$sample) %>%
  map(function(df) pull(.data=df, var=n)) %>%
  map(~table(.)) %>%
  map(~as_tibble(.)) %>%
  map(function(df) rename(.data=df, count = 1, f = 2)) %>%
  map(function(df) mutate(.data=df, count = as.numeric(count))) %>%
  map(~breakaway(.)) 

bas %>%
  alpha_estimates %>% 
  summary %>%
  inner_join(dataframe %>%
               select(sample, dairy) %>% 
               distinct, c("sample_names" = "sample")) %>%
  betta(formula = estimate ~ dairy, ses = error, data = .)
  


### Open questions: 
### What to do with these duplicates? Are they ok for species richness? Can we do strain richness? 

dataframe %>%
  mutate(taxon = paste(t_domain, t_phylum, t_class, t_order, t_family, t_genus, 
                       t_species)) %>%
  select(-c("t_domain", "t_phylum", "t_class", "t_order", "t_family", "t_genus", 
            "t_species")) %>%
  group_by(taxon, sample) %>%
  filter(sample == "C02_S104", taxon == "Bacteria Bacteroidota Bacteroidia Bacteroidales Bacteroidaceae Bacteroides None") %>%
  filter(gene_name == "Ribosomal_L22")


dataframe %>%
  mutate(taxon = paste(t_domain, t_phylum, t_class, t_order, t_family, t_genus, 
                       t_species)) %>%
  select(-c("t_domain", "t_phylum", "t_class", "t_order", "t_family", "t_genus", 
            "t_species")) %>%
  group_by(taxon, sample) %>%
  filter(sample == "C02_S104", taxon == "Bacteria Firmicutes Clostridia Oscillospirales Ruminococcaceae Faecalibacterium Faecalibacterium prausnitzii") %>%
  filter(gene_name == "Ribosomal_L22")
