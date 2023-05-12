# cephalosporins 

#####
tax_names <- read.csv("/Users/paulinetrinh/Documents/HDW/metagenomic/metacherchant_results/kraken2_seqs/C02_S106_genecall_186953.seqs", sep = "\t", header = F ) %>% 
  select(V3) %>% 
  unique() %>% 
  separate(V3, c("Taxa", "other"), sep = " [(]taxid ", extra = "merge") %>% 
  separate(other, c("taxid","other"), sep = "[)]") %>% 
  select(-other) %>% 
  as_tibble() 

node_info <- read.csv("/Users/paulinetrinh/Documents/HDW/metagenomic/metacherchant_results/kraken2_CIA_results/C02_S106_genecall_186953.classified", sep = "\t", header = F ) %>%
  filter(grepl(">", V1)) %>% 
  separate(V1, c("dropme","id","therest"),sep = " ", extra = "merge") %>% 
  select(-dropme) %>% 
  separate(therest, c("drop_these","taxid"), sep = "kraken:taxid[|]", extra = "merge") %>%  
  separate(id, c("dropthis","Node name"), sep = "Id") %>% 
  select(-c(drop_these,dropthis)) %>% 
  as_tibble() %>% 
  full_join(tax_names) %>% 
  filter(!Taxa == "root") %>% 
  filter(!Taxa %in% c("Bacteroidetes","Bacteroidales","Firmicutes"))

write.csv(node_info,"bandage_annotations/C02_S106_genecall_186953_annotation.csv", row.names = F, quote = F)
######


#####
tax_names <- read.csv("/Users/paulinetrinh/Documents/HDW/metagenomic/metacherchant_results/kraken2_seqs/D03_S26_genecall_164534.seqs", sep = "\t", header = F ) %>% 
  select(V3) %>% 
  unique() %>% 
  separate(V3, c("Taxa", "other"), sep = " [(]taxid ", extra = "merge") %>% 
  separate(other, c("taxid","other"), sep = "[)]") %>% 
  select(-other) %>% 
  as_tibble() 

node_info <- read.csv("/Users/paulinetrinh/Documents/HDW/metagenomic/metacherchant_results/kraken2_CIA_results/D03_S26_genecall_164534.classified", sep = "\t", header = F ) %>%
  filter(grepl(">", V1)) %>% 
  separate(V1, c("dropme","id","therest"),sep = " ", extra = "merge") %>% 
  select(-dropme) %>% 
  separate(therest, c("drop_these","taxid"), sep = "kraken:taxid[|]", extra = "merge") %>%  
  separate(id, c("dropthis","Node name"), sep = "Id") %>% 
  select(-c(drop_these,dropthis)) %>% 
  as_tibble() %>% 
  full_join(tax_names) %>% 
  filter(!Taxa == "root") %>% 
  filter(!Taxa %in% c("Bacteroidetes","Bacteroidales","Firmicutes"))

write.csv(node_info,"bandage_annotations/D03_S26_genecall_164534_annotation.csv", row.names = F, quote = F)
######


#####
tax_names <- read.csv("/Users/paulinetrinh/Documents/HDW/metagenomic/metacherchant_results/kraken2_seqs/D03_S32_genecall_4531.seqs", sep = "\t", header = F ) %>% 
  select(V3) %>% 
  unique() %>% 
  separate(V3, c("Taxa", "other"), sep = " [(]taxid ", extra = "merge") %>% 
  separate(other, c("taxid","other"), sep = "[)]") %>% 
  select(-other) %>% 
  as_tibble() 

node_info <- read.csv("/Users/paulinetrinh/Documents/HDW/metagenomic/metacherchant_results/kraken2_CIA_results/D03_S32_genecall_4531.classified", sep = "\t", header = F ) %>%
  filter(grepl(">", V1)) %>% 
  separate(V1, c("dropme","id","therest"),sep = " ", extra = "merge") %>% 
  select(-dropme) %>% 
  separate(therest, c("drop_these","taxid"), sep = "kraken:taxid[|]", extra = "merge") %>%  
  separate(id, c("dropthis","Node name"), sep = "Id") %>% 
  select(-c(drop_these,dropthis)) %>% 
  as_tibble() %>% 
  full_join(tax_names) %>% 
  filter(!Taxa == "root") %>% 
  filter(!Taxa %in% c("Bacteroidetes","Bacteroidales","Firmicutes"))

write.csv(node_info,"bandage_annotations/D03_S32_genecall_4531_annotation.csv", row.names = F, quote = F)
######
