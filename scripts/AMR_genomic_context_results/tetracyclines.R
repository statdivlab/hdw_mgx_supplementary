# tetracylines 
#####
tax_names <- read.csv("/Users/paulinetrinh/Documents/HDW/metagenomic/metacherchant_results/kraken2_seqs/D03_S25_genecall_1119.seqs", sep = "\t", header = F ) %>% 
  select(V3) %>% 
  unique() %>% 
  separate(V3, c("Taxa", "other"), sep = " [(]taxid ", extra = "merge") %>% 
  separate(other, c("taxid","other"), sep = "[)]") %>% 
  select(-other) %>% 
  as_tibble() 

node_info <- read.csv("/Users/paulinetrinh/Documents/HDW/metagenomic/metacherchant_results/kraken2_CIA_results/D03_S25_genecall_1119.classified", sep = "\t", header = F ) %>%
  filter(grepl(">", V1)) %>% 
  separate(V1, c("dropme","id","therest"),sep = " ", extra = "merge") %>% 
  select(-dropme) %>% 
  separate(therest, c("drop_these","taxid"), sep = "kraken:taxid[|]", extra = "merge") %>%  
  separate(id, c("dropthis","Node name"), sep = "Id") %>% 
  select(-c(drop_these,dropthis)) %>% 
  as_tibble() %>% 
  full_join(tax_names) %>% 
  filter(!Taxa == "root") %>% 
  filter(!Taxa %in% c("Bacteroidetes","Bacteroidales","Firmicutes","Bacteria","Terrabacteria group"))

write.csv(node_info,"bandage_annotations/D03_S25_genecall_1119_annotation.csv", row.names = F, quote = F)
######


#####
tax_names <- read.csv("/Users/paulinetrinh/Documents/HDW/metagenomic/metacherchant_results/kraken2_seqs/D03_S35_genecall_139914.seqs", sep = "\t", header = F ) %>% 
  select(V3) %>% 
  unique() %>% 
  separate(V3, c("Taxa", "other"), sep = " [(]taxid ", extra = "merge") %>% 
  separate(other, c("taxid","other"), sep = "[)]") %>% 
  select(-other) %>% 
  as_tibble() 

node_info <- read.csv("/Users/paulinetrinh/Documents/HDW/metagenomic/metacherchant_results/kraken2_CIA_results/D03_S35_genecall_139914.classified", sep = "\t", header = F ) %>%
  filter(grepl(">", V1)) %>% 
  separate(V1, c("dropme","id","therest"),sep = " ", extra = "merge") %>% 
  select(-dropme) %>% 
  separate(therest, c("drop_these","taxid"), sep = "kraken:taxid[|]", extra = "merge") %>%  
  separate(id, c("dropthis","Node name"), sep = "Id") %>% 
  select(-c(drop_these,dropthis)) %>% 
  as_tibble() %>% 
  full_join(tax_names) %>% 
  filter(!Taxa == "root") %>% 
  filter(!Taxa %in% c("Eubacteriales","Bacteroidetes","Bacteroidales","Firmicutes","Bacteria","Terrabacteria group","Actinobacteria"))

write.csv(node_info,"bandage_annotations/D03_S35_genecall_139914_annotation.csv", row.names = F, quote = F)
######


### tetM 
#####
tax_names <- read.csv("/Users/paulinetrinh/Documents/HDW/metagenomic/metacherchant_results/kraken2_seqs/C02_S113_genecall_79810.seqs", sep = "\t", header = F ) %>% 
  select(V3) %>% 
  unique() %>% 
  separate(V3, c("Taxa", "other"), sep = " [(]taxid ", extra = "merge") %>% 
  separate(other, c("taxid","other"), sep = "[)]") %>% 
  select(-other) %>% 
  as_tibble() 

node_info <- read.csv("/Users/paulinetrinh/Documents/HDW/metagenomic/metacherchant_results/kraken2_CIA_results/C02_S113_genecall_79810.classified", sep = "\t", header = F ) %>%
  filter(grepl(">", V1)) %>% 
  separate(V1, c("dropme","id","therest"),sep = " ", extra = "merge") %>% 
  select(-dropme) %>% 
  separate(therest, c("drop_these","taxid"), sep = "kraken:taxid[|]", extra = "merge") %>%  
  separate(id, c("dropthis","Node name"), sep = "Id") %>% 
  select(-c(drop_these,dropthis)) %>% 
  as_tibble() %>% 
  full_join(tax_names) %>% 
  filter(!Taxa == "root") %>% 
  filter(!Taxa %in% c("Eubacteriales","Bacteroidetes","Bacteroidales","Firmicutes","Bacteria","Terrabacteria group","Actinobacteria"))

write.csv(node_info,"bandage_annotations/C02_S113_genecall_79810_annotation.csv", row.names = F, quote = F)
######


#####
tax_names <- read.csv("/Users/paulinetrinh/Documents/HDW/metagenomic/metacherchant_results/kraken2_seqs/D03_S27_genecall_91374.seqs", sep = "\t", header = F ) %>% 
  select(V3) %>% 
  unique() %>% 
  separate(V3, c("Taxa", "other"), sep = " [(]taxid ", extra = "merge") %>% 
  separate(other, c("taxid","other"), sep = "[)]") %>% 
  select(-other) %>% 
  as_tibble() 

node_info <- read.csv("/Users/paulinetrinh/Documents/HDW/metagenomic/metacherchant_results/kraken2_CIA_results/D03_S27_genecall_91374.classified", sep = "\t", header = F ) %>%
  filter(grepl(">", V1)) %>% 
  separate(V1, c("dropme","id","therest"),sep = " ", extra = "merge") %>% 
  select(-dropme) %>% 
  separate(therest, c("drop_these","taxid"), sep = "kraken:taxid[|]", extra = "merge") %>%  
  separate(id, c("dropthis","Node name"), sep = "Id") %>% 
  select(-c(drop_these,dropthis)) %>% 
  as_tibble() %>% 
  full_join(tax_names) %>% 
  filter(!Taxa == "root") %>% 
  filter(!Taxa %in% c("Bacteria"))

write.csv(node_info,"bandage_annotations/D03_S27_genecall_91374_annotation.csv", row.names = F, quote = F)
######

#####
tax_names <- read.csv("/Users/paulinetrinh/Documents/HDW/metagenomic/metacherchant_results/kraken2_seqs/D03_S33_genecall_88214.seqs", sep = "\t", header = F ) %>% 
  select(V3) %>% 
  unique() %>% 
  separate(V3, c("Taxa", "other"), sep = " [(]taxid ", extra = "merge") %>% 
  separate(other, c("taxid","other"), sep = "[)]") %>% 
  select(-other) %>% 
  as_tibble() 

node_info <- read.csv("/Users/paulinetrinh/Documents/HDW/metagenomic/metacherchant_results/kraken2_CIA_results/D03_S33_genecall_88214.classified", sep = "\t", header = F ) %>%
  filter(grepl(">", V1)) %>% 
  separate(V1, c("dropme","id","therest"),sep = " ", extra = "merge") %>% 
  select(-dropme) %>% 
  separate(therest, c("drop_these","taxid"), sep = "kraken:taxid[|]", extra = "merge") %>%  
  separate(id, c("dropthis","Node name"), sep = "Id") %>% 
  select(-c(drop_these,dropthis)) %>% 
  as_tibble() %>% 
  full_join(tax_names) %>% 
  filter(!Taxa == "root") %>% 
  filter(!Taxa %in% c("Eubacteriales","Bacteroidetes","Bacteroidales","Firmicutes","Bacteria","Terrabacteria group","Actinobacteria"))

write.csv(node_info,"bandage_annotations/D03_S33_genecall_88214_annotation.csv", row.names = F, quote = F)
######


## tetO

#####
tax_names <- read.csv("/Users/paulinetrinh/Documents/HDW/metagenomic/metacherchant_results/kraken2_seqs/C02_S115_genecall_249987.seqs", sep = "\t", header = F ) %>% 
  select(V3) %>% 
  unique() %>% 
  separate(V3, c("Taxa", "other"), sep = " [(]taxid ", extra = "merge") %>% 
  separate(other, c("taxid","other"), sep = "[)]") %>% 
  select(-other) %>% 
  as_tibble() 

node_info <- read.csv("/Users/paulinetrinh/Documents/HDW/metagenomic/metacherchant_results/kraken2_CIA_results/C02_S115_genecall_249987.classified", sep = "\t", header = F ) %>%
  filter(grepl(">", V1)) %>% 
  separate(V1, c("dropme","id","therest"),sep = " ", extra = "merge") %>% 
  select(-dropme) %>% 
  separate(therest, c("drop_these","taxid"), sep = "kraken:taxid[|]", extra = "merge") %>%  
  separate(id, c("dropthis","Node name"), sep = "Id") %>% 
  select(-c(drop_these,dropthis)) %>% 
  as_tibble() %>% 
  full_join(tax_names) %>% 
  filter(!Taxa == "root") %>% 
  filter(!Taxa %in% c("Eubacteriales","Bacteroidetes","Bacteroidales","Firmicutes","Bacteria","Terrabacteria group","Actinobacteria"))

write.csv(node_info,"bandage_annotations/C02_S115_genecall_249987_annotation.csv", row.names = F, quote = F)
######

#####
tax_names <- read.csv("/Users/paulinetrinh/Documents/HDW/metagenomic/metacherchant_results/kraken2_seqs/D03_S27_genecall_325.seqs", sep = "\t", header = F ) %>% 
  select(V3) %>% 
  unique() %>% 
  separate(V3, c("Taxa", "other"), sep = " [(]taxid ", extra = "merge") %>% 
  separate(other, c("taxid","other"), sep = "[)]") %>% 
  select(-other) %>% 
  as_tibble() 

node_info <- read.csv("/Users/paulinetrinh/Documents/HDW/metagenomic/metacherchant_results/kraken2_CIA_results/D03_S27_genecall_325.classified", sep = "\t", header = F ) %>%
  filter(grepl(">", V1)) %>% 
  separate(V1, c("dropme","id","therest"),sep = " ", extra = "merge") %>% 
  select(-dropme) %>% 
  separate(therest, c("drop_these","taxid"), sep = "kraken:taxid[|]", extra = "merge") %>%  
  separate(id, c("dropthis","Node name"), sep = "Id") %>% 
  select(-c(drop_these,dropthis)) %>% 
  as_tibble() %>% 
  full_join(tax_names) %>% 
  filter(!Taxa == "root") %>% 
  filter(!Taxa %in% c("Bacteria"))

write.csv(node_info,"bandage_annotations/D03_S27_genecall_325_annotation.csv", row.names = F, quote = F)
######


#####
tax_names <- read.csv("/Users/paulinetrinh/Documents/HDW/metagenomic/metacherchant_results/kraken2_seqs/D03_S35_genecall_8048.seqs", sep = "\t", header = F ) %>% 
  select(V3) %>% 
  unique() %>% 
  separate(V3, c("Taxa", "other"), sep = " [(]taxid ", extra = "merge") %>% 
  separate(other, c("taxid","other"), sep = "[)]") %>% 
  select(-other) %>% 
  as_tibble() 

node_info <- read.csv("/Users/paulinetrinh/Documents/HDW/metagenomic/metacherchant_results/kraken2_CIA_results/D03_S35_genecall_8048.classified", sep = "\t", header = F ) %>%
  filter(grepl(">", V1)) %>% 
  separate(V1, c("dropme","id","therest"),sep = " ", extra = "merge") %>% 
  select(-dropme) %>% 
  separate(therest, c("drop_these","taxid"), sep = "kraken:taxid[|]", extra = "merge") %>%  
  separate(id, c("dropthis","Node name"), sep = "Id") %>% 
  select(-c(drop_these,dropthis)) %>% 
  as_tibble() %>% 
  full_join(tax_names) %>% 
  filter(!Taxa == "root") %>% 
  filter(!Taxa %in% c("Eubacteriales","Bacteroidetes","Bacteroidales","Firmicutes","Bacteria","Terrabacteria group","Actinobacteria"))


write.csv(node_info,"bandage_annotations/D03_S35_genecall_8048_annotation.csv", row.names = F, quote = F)
######

#####
tax_names <- read.csv("/Users/paulinetrinh/Documents/HDW/metagenomic/metacherchant_results/kraken2_seqs/D03_S37_genecall_53234.seqs", sep = "\t", header = F ) %>% 
  select(V3) %>% 
  unique() %>% 
  separate(V3, c("Taxa", "other"), sep = " [(]taxid ", extra = "merge") %>% 
  separate(other, c("taxid","other"), sep = "[)]") %>% 
  select(-other) %>% 
  as_tibble() 

node_info <- read.csv("/Users/paulinetrinh/Documents/HDW/metagenomic/metacherchant_results/kraken2_CIA_results/D03_S37_genecall_53234.classified", sep = "\t", header = F ) %>%
  filter(grepl(">", V1)) %>% 
  separate(V1, c("dropme","id","therest"),sep = " ", extra = "merge") %>% 
  select(-dropme) %>% 
  separate(therest, c("drop_these","taxid"), sep = "kraken:taxid[|]", extra = "merge") %>%  
  separate(id, c("dropthis","Node name"), sep = "Id") %>% 
  select(-c(drop_these,dropthis)) %>% 
  as_tibble() %>% 
  full_join(tax_names) %>% 
  filter(!Taxa == "root") %>% 
  filter(!Taxa %in% c("Eubacteriales","Bacteroidetes","Bacteroidales","Firmicutes","Bacteria","Terrabacteria group","Actinobacteria"))


write.csv(node_info,"bandage_annotations/D03_S37_genecall_53234_annotation.csv", row.names = F, quote = F)
######

######
# tetQ
######
#####
tax_names <- read.csv("/Users/paulinetrinh/Documents/HDW/metagenomic/metacherchant_results/kraken2_seqs/C02_S112_genecall_22182.seqs", sep = "\t", header = F ) %>% 
  select(V3) %>% 
  unique() %>% 
  separate(V3, c("Taxa", "other"), sep = " [(]taxid ", extra = "merge") %>% 
  separate(other, c("taxid","other"), sep = "[)]") %>% 
  select(-other) %>% 
  as_tibble() 

node_info <- read.csv("/Users/paulinetrinh/Documents/HDW/metagenomic/metacherchant_results/kraken2_CIA_results/C02_S112_genecall_22182.classified", sep = "\t", header = F ) %>%
  filter(grepl(">", V1)) %>% 
  separate(V1, c("dropme","id","therest"),sep = " ", extra = "merge") %>% 
  select(-dropme) %>% 
  separate(therest, c("drop_these","taxid"), sep = "kraken:taxid[|]", extra = "merge") %>%  
  separate(id, c("dropthis","Node name"), sep = "Id") %>% 
  select(-c(drop_these,dropthis)) %>% 
  as_tibble() %>% 
  full_join(tax_names) %>% 
  filter(!Taxa == "root") %>% 
  filter(!Taxa %in% c("Firmicutes","Bacteria","Bacteroidetes","Terrabacteria group","Actinobacteria"))


write.csv(node_info,"bandage_annotations/C02_S112_genecall_22182_annotation.csv", row.names = F, quote = F)
######


#####
tax_names <- read.csv("/Users/paulinetrinh/Documents/HDW/metagenomic/metacherchant_results/kraken2_seqs/C02_S115_genecall_6473.seqs", sep = "\t", header = F ) %>% 
  select(V3) %>% 
  unique() %>% 
  separate(V3, c("Taxa", "other"), sep = " [(]taxid ", extra = "merge") %>% 
  separate(other, c("taxid","other"), sep = "[)]") %>% 
  select(-other) %>% 
  as_tibble() 

node_info <- read.csv("/Users/paulinetrinh/Documents/HDW/metagenomic/metacherchant_results/kraken2_CIA_results/C02_S115_genecall_6473.classified", sep = "\t", header = F ) %>%
  filter(grepl(">", V1)) %>% 
  separate(V1, c("dropme","id","therest"),sep = " ", extra = "merge") %>% 
  select(-dropme) %>% 
  separate(therest, c("drop_these","taxid"), sep = "kraken:taxid[|]", extra = "merge") %>%  
  separate(id, c("dropthis","Node name"), sep = "Id") %>% 
  select(-c(drop_these,dropthis)) %>% 
  as_tibble() %>% 
  full_join(tax_names) %>% 
  filter(!Taxa == "root") %>% 
  filter(!Taxa %in% c("Firmicutes","Bacteria","Bacteroidetes","Terrabacteria group","Actinobacteria"))


write.csv(node_info,"bandage_annotations/C02_S115_genecall_6473_annotation.csv", row.names = F, quote = F)
######


#####
tax_names <- read.csv("/Users/paulinetrinh/Documents/HDW/metagenomic/metacherchant_results/kraken2_seqs/D03_S26_genecall_4527.seqs", sep = "\t", header = F ) %>% 
  select(V3) %>% 
  unique() %>% 
  separate(V3, c("Taxa", "other"), sep = " [(]taxid ", extra = "merge") %>% 
  separate(other, c("taxid","other"), sep = "[)]") %>% 
  select(-other) %>% 
  as_tibble() 

node_info <- read.csv("/Users/paulinetrinh/Documents/HDW/metagenomic/metacherchant_results/kraken2_CIA_results/D03_S26_genecall_4527.classified", sep = "\t", header = F ) %>%
  filter(grepl(">", V1)) %>% 
  separate(V1, c("dropme","id","therest"),sep = " ", extra = "merge") %>% 
  select(-dropme) %>% 
  separate(therest, c("drop_these","taxid"), sep = "kraken:taxid[|]", extra = "merge") %>%  
  separate(id, c("dropthis","Node name"), sep = "Id") %>% 
  select(-c(drop_these,dropthis)) %>% 
  as_tibble() %>% 
  full_join(tax_names) %>% 
  filter(!Taxa == "root") %>% 
  filter(!Taxa %in% c("Firmicutes","Bacteria","Bacteroidetes","Terrabacteria group","Actinobacteria"))


write.csv(node_info,"bandage_annotations/D03_S26_genecall_4527_annotation.csv", row.names = F, quote = F)
######

#####
tax_names <- read.csv("/Users/paulinetrinh/Documents/HDW/metagenomic/metacherchant_results/kraken2_seqs/D03_S37_genecall_60628.seqs", sep = "\t", header = F ) %>% 
  select(V3) %>% 
  unique() %>% 
  separate(V3, c("Taxa", "other"), sep = " [(]taxid ", extra = "merge") %>% 
  separate(other, c("taxid","other"), sep = "[)]") %>% 
  select(-other) %>% 
  as_tibble() 

node_info <- read.csv("/Users/paulinetrinh/Documents/HDW/metagenomic/metacherchant_results/kraken2_CIA_results/D03_S37_genecall_60628.classified", sep = "\t", header = F ) %>%
  filter(grepl(">", V1)) %>% 
  separate(V1, c("dropme","id","therest"),sep = " ", extra = "merge") %>% 
  select(-dropme) %>% 
  separate(therest, c("drop_these","taxid"), sep = "kraken:taxid[|]", extra = "merge") %>%  
  separate(id, c("dropthis","Node name"), sep = "Id") %>% 
  select(-c(drop_these,dropthis)) %>% 
  as_tibble() %>% 
  full_join(tax_names) %>% 
  filter(!Taxa == "root") %>% 
  filter(!Taxa %in% c("Firmicutes","Bacteria","Bacteroidetes","Terrabacteria group","Actinobacteria"))


write.csv(node_info,"bandage_annotations/D03_S37_genecall_60628_annotation.csv", row.names = F, quote = F)
######

### 
# tetW
###
#####
tax_names <- read.csv("/Users/paulinetrinh/Documents/HDW/metagenomic/metacherchant_results/kraken2_seqs/D03_S31_genecall_66817.seqs", sep = "\t", header = F ) %>% 
  select(V3) %>% 
  unique() %>% 
  separate(V3, c("Taxa", "other"), sep = " [(]taxid ", extra = "merge") %>% 
  separate(other, c("taxid","other"), sep = "[)]") %>% 
  select(-other) %>% 
  as_tibble() 

node_info <- read.csv("/Users/paulinetrinh/Documents/HDW/metagenomic/metacherchant_results/kraken2_CIA_results/D03_S31_genecall_66817.classified", sep = "\t", header = F ) %>%
  filter(grepl(">", V1)) %>% 
  separate(V1, c("dropme","id","therest"),sep = " ", extra = "merge") %>% 
  select(-dropme) %>% 
  separate(therest, c("drop_these","taxid"), sep = "kraken:taxid[|]", extra = "merge") %>%  
  separate(id, c("dropthis","Node name"), sep = "Id") %>% 
  select(-c(drop_these,dropthis)) %>% 
  as_tibble() %>% 
  full_join(tax_names) %>% 
  filter(!Taxa == "root") %>% 
  filter(!Taxa %in% c("Firmicutes","Bacteria","Bacteroidetes","Terrabacteria group","Actinobacteria"))


write.csv(node_info,"bandage_annotations/D03_S31_genecall_66817_annotation.csv", row.names = F, quote = F)
######

#####
tax_names <- read.csv("/Users/paulinetrinh/Documents/HDW/metagenomic/metacherchant_results/kraken2_seqs/D03_S32_genecall_21166.seqs", sep = "\t", header = F ) %>% 
  select(V3) %>% 
  unique() %>% 
  separate(V3, c("Taxa", "other"), sep = " [(]taxid ", extra = "merge") %>% 
  separate(other, c("taxid","other"), sep = "[)]") %>% 
  select(-other) %>% 
  as_tibble() 

node_info <- read.csv("/Users/paulinetrinh/Documents/HDW/metagenomic/metacherchant_results/kraken2_CIA_results/D03_S32_genecall_21166.classified", sep = "\t", header = F ) %>%
  filter(grepl(">", V1)) %>% 
  separate(V1, c("dropme","id","therest"),sep = " ", extra = "merge") %>% 
  select(-dropme) %>% 
  separate(therest, c("drop_these","taxid"), sep = "kraken:taxid[|]", extra = "merge") %>%  
  separate(id, c("dropthis","Node name"), sep = "Id") %>% 
  select(-c(drop_these,dropthis)) %>% 
  as_tibble() %>% 
  full_join(tax_names) %>% 
  filter(!Taxa == "root") %>% 
  filter(!Taxa %in% c("Firmicutes","Bacteria","Bacteroidetes","Terrabacteria group","Actinobacteria"))


write.csv(node_info,"bandage_annotations/D03_S32_genecall_21166_annotation.csv", row.names = F, quote = F)
######
#####
tax_names <- read.csv("/Users/paulinetrinh/Documents/HDW/metagenomic/metacherchant_results/kraken2_seqs/D03_S37_genecall_18323.seqs", sep = "\t", header = F ) %>% 
  select(V3) %>% 
  unique() %>% 
  separate(V3, c("Taxa", "other"), sep = " [(]taxid ", extra = "merge") %>% 
  separate(other, c("taxid","other"), sep = "[)]") %>% 
  select(-other) %>% 
  as_tibble() 

node_info <- read.csv("/Users/paulinetrinh/Documents/HDW/metagenomic/metacherchant_results/kraken2_CIA_results/D03_S37_genecall_18323.classified", sep = "\t", header = F ) %>%
  filter(grepl(">", V1)) %>% 
  separate(V1, c("dropme","id","therest"),sep = " ", extra = "merge") %>% 
  select(-dropme) %>% 
  separate(therest, c("drop_these","taxid"), sep = "kraken:taxid[|]", extra = "merge") %>%  
  separate(id, c("dropthis","Node name"), sep = "Id") %>% 
  select(-c(drop_these,dropthis)) %>% 
  as_tibble() %>% 
  full_join(tax_names) %>% 
  filter(!Taxa == "root") %>% 
  filter(!Taxa %in% c("Firmicutes","Bacteria","Bacteroidetes","Terrabacteria group","Actinobacteria"))


write.csv(node_info,"bandage_annotations/D03_S37_genecall_18323_annotation.csv", row.names = F, quote = F)
######


### 
# tet32 
#####
tax_names <- read.csv("/Users/paulinetrinh/Documents/HDW/metagenomic/metacherchant_results/kraken2_seqs/C02_S113_genecall_197758.seqs", sep = "\t", header = F ) %>% 
  select(V3) %>% 
  unique() %>% 
  separate(V3, c("Taxa", "other"), sep = " [(]taxid ", extra = "merge") %>% 
  separate(other, c("taxid","other"), sep = "[)]") %>% 
  select(-other) %>% 
  as_tibble() 

node_info <- read.csv("/Users/paulinetrinh/Documents/HDW/metagenomic/metacherchant_results/kraken2_CIA_results/C02_S113_genecall_197758.classified", sep = "\t", header = F ) %>%
  filter(grepl(">", V1)) %>% 
  separate(V1, c("dropme","id","therest"),sep = " ", extra = "merge") %>% 
  select(-dropme) %>% 
  separate(therest, c("drop_these","taxid"), sep = "kraken:taxid[|]", extra = "merge") %>%  
  separate(id, c("dropthis","Node name"), sep = "Id") %>% 
  select(-c(drop_these,dropthis)) %>% 
  as_tibble() %>% 
  full_join(tax_names) %>% 
  filter(!Taxa == "root") %>% 
  filter(!Taxa %in% c("Firmicutes","Bacteria","Bacteroidetes","Terrabacteria group","Actinobacteria"))


write.csv(node_info,"bandage_annotations/C02_S113_genecall_197758_annotation.csv", row.names = F, quote = F)
######

#####
tax_names <- read.csv("/Users/paulinetrinh/Documents/HDW/metagenomic/metacherchant_results/kraken2_seqs/C02_S115_genecall_147095.seqs", sep = "\t", header = F ) %>% 
  select(V3) %>% 
  unique() %>% 
  separate(V3, c("Taxa", "other"), sep = " [(]taxid ", extra = "merge") %>% 
  separate(other, c("taxid","other"), sep = "[)]") %>% 
  select(-other) %>% 
  as_tibble() 

node_info <- read.csv("/Users/paulinetrinh/Documents/HDW/metagenomic/metacherchant_results/kraken2_CIA_results/C02_S115_genecall_147095.classified", sep = "\t", header = F ) %>%
  filter(grepl(">", V1)) %>% 
  separate(V1, c("dropme","id","therest"),sep = " ", extra = "merge") %>% 
  select(-dropme) %>% 
  separate(therest, c("drop_these","taxid"), sep = "kraken:taxid[|]", extra = "merge") %>%  
  separate(id, c("dropthis","Node name"), sep = "Id") %>% 
  select(-c(drop_these,dropthis)) %>% 
  as_tibble() %>% 
  full_join(tax_names) %>% 
  filter(!Taxa == "root") %>% 
  filter(!Taxa %in% c("Firmicutes","Bacteria","Bacteroidetes","Terrabacteria group","Actinobacteria"))


write.csv(node_info,"bandage_annotations/C02_S115_genecall_147095_annotation.csv", row.names = F, quote = F)
######

##### 
# tetB 
#####
tax_names <- read.csv("/Users/paulinetrinh/Documents/HDW/metagenomic/metacherchant_results/kraken2_seqs/C02_S111_genecall_82947.seqs", sep = "\t", header = F ) %>% 
  select(V3) %>% 
  unique() %>% 
  separate(V3, c("Taxa", "other"), sep = " [(]taxid ", extra = "merge") %>% 
  separate(other, c("taxid","other"), sep = "[)]") %>% 
  select(-other) %>% 
  as_tibble() 

node_info <- read.csv("/Users/paulinetrinh/Documents/HDW/metagenomic/metacherchant_results/kraken2_CIA_results/C02_S111_genecall_82947.classified", sep = "\t", header = F ) %>%
  filter(grepl(">", V1)) %>% 
  separate(V1, c("dropme","id","therest"),sep = " ", extra = "merge") %>% 
  select(-dropme) %>% 
  separate(therest, c("drop_these","taxid"), sep = "kraken:taxid[|]", extra = "merge") %>%  
  separate(id, c("dropthis","Node name"), sep = "Id") %>% 
  select(-c(drop_these,dropthis)) %>% 
  as_tibble() %>% 
  full_join(tax_names) %>% 
  filter(!Taxa == "root") %>% 
  filter(!Taxa %in% c("Firmicutes","Bacteria","Bacteroidetes","Terrabacteria group","Actinobacteria"))


write.csv(node_info,"bandage_annotations/C02_S111_genecall_82947_annotation.csv", row.names = F, quote = F)
######


# tetG 
#####
tax_names <- read.csv("/Users/paulinetrinh/Documents/HDW/metagenomic/metacherchant_results/kraken2_seqs/D03_S27_genecall_78879.seqs", sep = "\t", header = F ) %>% 
  select(V3) %>% 
  unique() %>% 
  separate(V3, c("Taxa", "other"), sep = " [(]taxid ", extra = "merge") %>% 
  separate(other, c("taxid","other"), sep = "[)]") %>% 
  select(-other) %>% 
  as_tibble() 

node_info <- read.csv("/Users/paulinetrinh/Documents/HDW/metagenomic/metacherchant_results/kraken2_CIA_results/D03_S27_genecall_78879.classified", sep = "\t", header = F ) %>%
  filter(grepl(">", V1)) %>% 
  separate(V1, c("dropme","id","therest"),sep = " ", extra = "merge") %>% 
  select(-dropme) %>% 
  separate(therest, c("drop_these","taxid"), sep = "kraken:taxid[|]", extra = "merge") %>%  
  separate(id, c("dropthis","Node name"), sep = "Id") %>% 
  select(-c(drop_these,dropthis)) %>% 
  as_tibble() %>% 
  full_join(tax_names) %>% 
  filter(!Taxa == "root") %>% 
  filter(!Taxa %in% c("Firmicutes","Bacteria","Bacteroidetes","Terrabacteria group","Actinobacteria"))


write.csv(node_info,"bandage_annotations/D03_S27_genecall_78879_annotation.csv", row.names = F, quote = F)
######


##### 
# tet(40)
#####
tax_names <- read.csv("/Users/paulinetrinh/Documents/HDW/metagenomic/metacherchant_results/kraken2_seqs/C02_S111_genecall_148720.seqs", sep = "\t", header = F ) %>% 
  select(V3) %>% 
  unique() %>% 
  separate(V3, c("Taxa", "other"), sep = " [(]taxid ", extra = "merge") %>% 
  separate(other, c("taxid","other"), sep = "[)]") %>% 
  select(-other) %>% 
  as_tibble() 

node_info <- read.csv("/Users/paulinetrinh/Documents/HDW/metagenomic/metacherchant_results/kraken2_CIA_results/C02_S111_genecall_148720.classified", sep = "\t", header = F ) %>%
  filter(grepl(">", V1)) %>% 
  separate(V1, c("dropme","id","therest"),sep = " ", extra = "merge") %>% 
  select(-dropme) %>% 
  separate(therest, c("drop_these","taxid"), sep = "kraken:taxid[|]", extra = "merge") %>%  
  separate(id, c("dropthis","Node name"), sep = "Id") %>% 
  select(-c(drop_these,dropthis)) %>% 
  as_tibble() %>% 
  full_join(tax_names) %>% 
  filter(!Taxa == "root") %>% 
  filter(!Taxa %in% c("Firmicutes","Bacteria","Bacteroidetes","Terrabacteria group","Actinobacteria"))


write.csv(node_info,"bandage_annotations/C02_S111_genecall_148720_annotation.csv", row.names = F, quote = F)
######

#####
tax_names <- read.csv("/Users/paulinetrinh/Documents/HDW/metagenomic/metacherchant_results/kraken2_seqs/C02_S112_genecall_58368.seqs", sep = "\t", header = F ) %>% 
  select(V3) %>% 
  unique() %>% 
  separate(V3, c("Taxa", "other"), sep = " [(]taxid ", extra = "merge") %>% 
  separate(other, c("taxid","other"), sep = "[)]") %>% 
  select(-other) %>% 
  as_tibble() 

node_info <- read.csv("/Users/paulinetrinh/Documents/HDW/metagenomic/metacherchant_results/kraken2_CIA_results/C02_S112_genecall_58368.classified", sep = "\t", header = F ) %>%
  filter(grepl(">", V1)) %>% 
  separate(V1, c("dropme","id","therest"),sep = " ", extra = "merge") %>% 
  select(-dropme) %>% 
  separate(therest, c("drop_these","taxid"), sep = "kraken:taxid[|]", extra = "merge") %>%  
  separate(id, c("dropthis","Node name"), sep = "Id") %>% 
  select(-c(drop_these,dropthis)) %>% 
  as_tibble() %>% 
  full_join(tax_names) %>% 
  filter(!Taxa == "root") %>% 
  filter(!Taxa %in% c("Firmicutes","Bacteria","Bacteroidetes","Terrabacteria group","Actinobacteria","Eubacteriales"))


write.csv(node_info,"bandage_annotations/C02_S112_genecall_58368_annotation.csv", row.names = F, quote = F)
######


#####
tax_names <- read.csv("/Users/paulinetrinh/Documents/HDW/metagenomic/metacherchant_results/kraken2_seqs/C02_S113_genecall_45611.seqs", sep = "\t", header = F ) %>% 
  select(V3) %>% 
  unique() %>% 
  separate(V3, c("Taxa", "other"), sep = " [(]taxid ", extra = "merge") %>% 
  separate(other, c("taxid","other"), sep = "[)]") %>% 
  select(-other) %>% 
  as_tibble() 

node_info <- read.csv("/Users/paulinetrinh/Documents/HDW/metagenomic/metacherchant_results/kraken2_CIA_results/C02_S113_genecall_45611.classified", sep = "\t", header = F ) %>%
  filter(grepl(">", V1)) %>% 
  separate(V1, c("dropme","id","therest"),sep = " ", extra = "merge") %>% 
  select(-dropme) %>% 
  separate(therest, c("drop_these","taxid"), sep = "kraken:taxid[|]", extra = "merge") %>%  
  separate(id, c("dropthis","Node name"), sep = "Id") %>% 
  select(-c(drop_these,dropthis)) %>% 
  as_tibble() %>% 
  full_join(tax_names) %>% 
  filter(!Taxa == "root") %>% 
  filter(!Taxa %in% c("Firmicutes","Bacteria","Bacteroidetes","Terrabacteria group","Actinobacteria","Eubacteriales"))


write.csv(node_info,"bandage_annotations/C02_S113_genecall_45611_annotation.csv", row.names = F, quote = F)
######

#####
tax_names <- read.csv("/Users/paulinetrinh/Documents/HDW/metagenomic/metacherchant_results/kraken2_seqs/C02_S115_genecall_56302.seqs", sep = "\t", header = F ) %>% 
  select(V3) %>% 
  unique() %>% 
  separate(V3, c("Taxa", "other"), sep = " [(]taxid ", extra = "merge") %>% 
  separate(other, c("taxid","other"), sep = "[)]") %>% 
  select(-other) %>% 
  as_tibble() 

node_info <- read.csv("/Users/paulinetrinh/Documents/HDW/metagenomic/metacherchant_results/kraken2_CIA_results/C02_S115_genecall_56302.classified", sep = "\t", header = F ) %>%
  filter(grepl(">", V1)) %>% 
  separate(V1, c("dropme","id","therest"),sep = " ", extra = "merge") %>% 
  select(-dropme) %>% 
  separate(therest, c("drop_these","taxid"), sep = "kraken:taxid[|]", extra = "merge") %>%  
  separate(id, c("dropthis","Node name"), sep = "Id") %>% 
  select(-c(drop_these,dropthis)) %>% 
  as_tibble() %>% 
  full_join(tax_names) %>% 
  filter(!Taxa == "root") %>% 
  filter(!Taxa %in% c("Firmicutes","Bacteria","Bacteroidetes","Terrabacteria group","Actinobacteria","Eubacteriales"))


write.csv(node_info,"bandage_annotations/C02_S115_genecall_56302_annotation.csv", row.names = F, quote = F)
######

#####
tax_names <- read.csv("/Users/paulinetrinh/Documents/HDW/metagenomic/metacherchant_results/kraken2_seqs/D03_S25_genecall_53651.seqs", sep = "\t", header = F ) %>% 
  select(V3) %>% 
  unique() %>% 
  separate(V3, c("Taxa", "other"), sep = " [(]taxid ", extra = "merge") %>% 
  separate(other, c("taxid","other"), sep = "[)]") %>% 
  select(-other) %>% 
  as_tibble() 

node_info <- read.csv("/Users/paulinetrinh/Documents/HDW/metagenomic/metacherchant_results/kraken2_CIA_results/D03_S25_genecall_53651.classified", sep = "\t", header = F ) %>%
  filter(grepl(">", V1)) %>% 
  separate(V1, c("dropme","id","therest"),sep = " ", extra = "merge") %>% 
  select(-dropme) %>% 
  separate(therest, c("drop_these","taxid"), sep = "kraken:taxid[|]", extra = "merge") %>%  
  separate(id, c("dropthis","Node name"), sep = "Id") %>% 
  select(-c(drop_these,dropthis)) %>% 
  as_tibble() %>% 
  full_join(tax_names) %>% 
  filter(!Taxa == "root") %>% 
  filter(!Taxa %in% c("Firmicutes","Bacteria","Bacteroidetes","Terrabacteria group","Actinobacteria", "Eubacteriales"))


write.csv(node_info,"bandage_annotations/D03_S25_genecall_53651_annotation.csv", row.names = F, quote = F)
######

#####
tax_names <- read.csv("/Users/paulinetrinh/Documents/HDW/metagenomic/metacherchant_results/kraken2_seqs/D03_S27_genecall_324.seqs", sep = "\t", header = F ) %>% 
  select(V3) %>% 
  unique() %>% 
  separate(V3, c("Taxa", "other"), sep = " [(]taxid ", extra = "merge") %>% 
  separate(other, c("taxid","other"), sep = "[)]") %>% 
  select(-other) %>% 
  as_tibble() 

node_info <- read.csv("/Users/paulinetrinh/Documents/HDW/metagenomic/metacherchant_results/kraken2_CIA_results/D03_S27_genecall_324.classified", sep = "\t", header = F ) %>%
  filter(grepl(">", V1)) %>% 
  separate(V1, c("dropme","id","therest"),sep = " ", extra = "merge") %>% 
  select(-dropme) %>% 
  separate(therest, c("drop_these","taxid"), sep = "kraken:taxid[|]", extra = "merge") %>%  
  separate(id, c("dropthis","Node name"), sep = "Id") %>% 
  select(-c(drop_these,dropthis)) %>% 
  as_tibble() %>% 
  full_join(tax_names) %>% 
  filter(!Taxa == "root") %>% 
  filter(!Taxa %in% c("Firmicutes","Bacteria","Bacteroidetes","Terrabacteria group","Actinobacteria"))


write.csv(node_info,"bandage_annotations/D03_S27_genecall_324_annotation.csv", row.names = F, quote = F)
######

#####
tax_names <- read.csv("/Users/paulinetrinh/Documents/HDW/metagenomic/metacherchant_results/kraken2_seqs/D03_S32_genecall_75755.seqs", sep = "\t", header = F ) %>% 
  select(V3) %>% 
  unique() %>% 
  separate(V3, c("Taxa", "other"), sep = " [(]taxid ", extra = "merge") %>% 
  separate(other, c("taxid","other"), sep = "[)]") %>% 
  select(-other) %>% 
  as_tibble() 

node_info <- read.csv("/Users/paulinetrinh/Documents/HDW/metagenomic/metacherchant_results/kraken2_CIA_results/D03_S32_genecall_75755.classified", sep = "\t", header = F ) %>%
  filter(grepl(">", V1)) %>% 
  separate(V1, c("dropme","id","therest"),sep = " ", extra = "merge") %>% 
  select(-dropme) %>% 
  separate(therest, c("drop_these","taxid"), sep = "kraken:taxid[|]", extra = "merge") %>%  
  separate(id, c("dropthis","Node name"), sep = "Id") %>% 
  select(-c(drop_these,dropthis)) %>% 
  as_tibble() %>% 
  full_join(tax_names) %>% 
  filter(!Taxa == "root") %>% 
  filter(!Taxa %in% c("Firmicutes","Bacteria","Bacteroidetes","Terrabacteria group","Actinobacteria"))


write.csv(node_info,"bandage_annotations/D03_S32_genecall_75755_annotation.csv", row.names = F, quote = F)
######


#####
tax_names <- read.csv("/Users/paulinetrinh/Documents/HDW/metagenomic/metacherchant_results/kraken2_seqs/D03_S35_genecall_43778.seqs", sep = "\t", header = F ) %>% 
  select(V3) %>% 
  unique() %>% 
  separate(V3, c("Taxa", "other"), sep = " [(]taxid ", extra = "merge") %>% 
  separate(other, c("taxid","other"), sep = "[)]") %>% 
  select(-other) %>% 
  as_tibble() 

node_info <- read.csv("/Users/paulinetrinh/Documents/HDW/metagenomic/metacherchant_results/kraken2_CIA_results/D03_S35_genecall_43778.classified", sep = "\t", header = F ) %>%
  filter(grepl(">", V1)) %>% 
  separate(V1, c("dropme","id","therest"),sep = " ", extra = "merge") %>% 
  select(-dropme) %>% 
  separate(therest, c("drop_these","taxid"), sep = "kraken:taxid[|]", extra = "merge") %>%  
  separate(id, c("dropthis","Node name"), sep = "Id") %>% 
  select(-c(drop_these,dropthis)) %>% 
  as_tibble() %>% 
  full_join(tax_names) %>% 
  filter(!Taxa == "root") %>% 
  filter(!Taxa %in% c("Firmicutes","Bacteria","Bacteroidetes","Terrabacteria group","Actinobacteria"))


write.csv(node_info,"bandage_annotations/D03_S35_genecall_43778_annotation.csv", row.names = F, quote = F)
######

#####
tax_names <- read.csv("/Users/paulinetrinh/Documents/HDW/metagenomic/metacherchant_results/kraken2_seqs/D03_S37_genecall_18548.seqs", sep = "\t", header = F ) %>% 
  select(V3) %>% 
  unique() %>% 
  separate(V3, c("Taxa", "other"), sep = " [(]taxid ", extra = "merge") %>% 
  separate(other, c("taxid","other"), sep = "[)]") %>% 
  select(-other) %>% 
  as_tibble() 

node_info <- read.csv("/Users/paulinetrinh/Documents/HDW/metagenomic/metacherchant_results/kraken2_CIA_results/D03_S37_genecall_18548.classified", sep = "\t", header = F ) %>%
  filter(grepl(">", V1)) %>% 
  separate(V1, c("dropme","id","therest"),sep = " ", extra = "merge") %>% 
  select(-dropme) %>% 
  separate(therest, c("drop_these","taxid"), sep = "kraken:taxid[|]", extra = "merge") %>%  
  separate(id, c("dropthis","Node name"), sep = "Id") %>% 
  select(-c(drop_these,dropthis)) %>% 
  as_tibble() %>% 
  full_join(tax_names) %>% 
  filter(!Taxa == "root") %>% 
  filter(!Taxa %in% c("Firmicutes","Bacteria","Bacteroidetes","Terrabacteria group","Actinobacteria"))


write.csv(node_info,"bandage_annotations/D03_S37_genecall_18548_annotation.csv", row.names = F, quote = F)
######

### 
# emrK
###

#####
tax_names <- read.csv("/Users/paulinetrinh/Documents/HDW/metagenomic/metacherchant_results/kraken2_seqs/C02_S111_genecall_33194.seqs", sep = "\t", header = F ) %>% 
  select(V3) %>% 
  unique() %>% 
  separate(V3, c("Taxa", "other"), sep = " [(]taxid ", extra = "merge") %>% 
  separate(other, c("taxid","other"), sep = "[)]") %>% 
  select(-other) %>% 
  as_tibble() 

node_info <- read.csv("/Users/paulinetrinh/Documents/HDW/metagenomic/metacherchant_results/kraken2_CIA_results/C02_S111_genecall_33194.classified", sep = "\t", header = F ) %>%
  filter(grepl(">", V1)) %>% 
  separate(V1, c("dropme","id","therest"),sep = " ", extra = "merge") %>% 
  select(-dropme) %>% 
  separate(therest, c("drop_these","taxid"), sep = "kraken:taxid[|]", extra = "merge") %>%  
  separate(id, c("dropthis","Node name"), sep = "Id") %>% 
  select(-c(drop_these,dropthis)) %>% 
  as_tibble() %>% 
  full_join(tax_names) %>% 
  filter(!Taxa == "root") %>% 
  filter(!Taxa %in% c("Firmicutes","Bacteria","Bacteroidetes","Terrabacteria group","Actinobacteria"))


write.csv(node_info,"bandage_annotations/C02_S111_genecall_33194_annotation.csv", row.names = F, quote = F)
######


#####
tax_names <- read.csv("/Users/paulinetrinh/Documents/HDW/metagenomic/metacherchant_results/kraken2_seqs/C02_S111_genecall_33195.seqs", sep = "\t", header = F ) %>% 
  select(V3) %>% 
  unique() %>% 
  separate(V3, c("Taxa", "other"), sep = " [(]taxid ", extra = "merge") %>% 
  separate(other, c("taxid","other"), sep = "[)]") %>% 
  select(-other) %>% 
  as_tibble() 

node_info <- read.csv("/Users/paulinetrinh/Documents/HDW/metagenomic/metacherchant_results/kraken2_CIA_results/C02_S111_genecall_33195.classified", sep = "\t", header = F ) %>%
  filter(grepl(">", V1)) %>% 
  separate(V1, c("dropme","id","therest"),sep = " ", extra = "merge") %>% 
  select(-dropme) %>% 
  separate(therest, c("drop_these","taxid"), sep = "kraken:taxid[|]", extra = "merge") %>%  
  separate(id, c("dropthis","Node name"), sep = "Id") %>% 
  select(-c(drop_these,dropthis)) %>% 
  as_tibble() %>% 
  full_join(tax_names) %>% 
  filter(!Taxa == "root") %>% 
  filter(!Taxa %in% c("Firmicutes","Bacteria","Bacteroidetes","Terrabacteria group","Actinobacteria"))


write.csv(node_info,"bandage_annotations/C02_S111_genecall_33195_annotation.csv", row.names = F, quote = F)
######



