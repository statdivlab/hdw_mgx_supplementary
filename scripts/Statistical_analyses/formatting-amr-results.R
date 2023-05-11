library(dplyr)
library(readr)
df <- list.files(path="/Users/paulinetrinh/Documents/HDW/metagenomic/contigs-amr-analysis/results", full.names = TRUE) %>% 
  lapply(read_tsv)
df_bound <- do.call("rbind",df) %>% as_tibble()

df_bound %>% 
  rename("gene_callers_id" = SEQUENCE,
         "source" = DATABASE,
         "accession" = ACCESSION, 
         "function" = GENE) %>% 
  mutate(e_value = 0) -> df_bound

View(df_bound)

df_bound %>% 
  filter(`%COVERAGE` == 100 | `%IDENTITY` == 90) %>% 
  select(`#FILE`,gene_callers_id,START,END,`%COVERAGE`,`%IDENTITY`, RESISTANCE, source, "function", accession, PRODUCT, e_value)-> complete_genes
dim(complete_genes)

saveRDS(complete_genes,"complete_amr_genes.RDS")
complete_genes <- readRDS("complete_amr_genes.RDS")

names <- list()
names <- complete_genes$`#FILE` %>% unique()
card <- complete_genes %>% filter(source == "card")
megares <- complete_genes %>% filter(source == "megares") 
hdw_meta <- read_tsv("/Users/paulinetrinh/Documents/HDW/metagenomic/hdw_mgx_metadata.txt")

vfdb_genes <- complete_genes %>% filter(source == "vfdb") %>% 
  mutate(filename = `#FILE`) %>% 
  separate(`#FILE`,c("sample","delete"), "-") %>% 
  mutate(sample_names = ifelse(sample == "C02_S104", "HDW-C02-S104",
                               ifelse(sample == "C02_S106","HDW-C02-S106",
                                      ifelse(sample == "C02_S111","HDW-C02-S111",
                                             ifelse(sample == "C02_S112", "HDW-C02-S112", 
                                                    ifelse(sample == "C02_S113","HDW-C02-S113",
                                                           ifelse(sample == "C02_S115","HDW-C02-S115",
                                                                  ifelse(sample == "D03_S25","HDW-D03-S25",
                                                                         ifelse(sample == "D03_S26","HDW-D03-S26",
                                                                                ifelse(sample == "D03_S27","HDW-D03-S27",
                                                                                       ifelse(sample == "D03_S28","HDW-D03-S28",
                                                                                              ifelse(sample == "D03_S31","HDW-D03-S31",
                                                                                                     ifelse(sample == "D03_S32","HDW-D03-S32",
                                                                                                            ifelse(sample == "D03_S33","HDW-D03-S33",
                                                                                                                   ifelse(sample == "D03_S34","HDW-D03-S34",
                                                                                                                          ifelse(sample == "D03_S35","HDW-D03-S35",
                                                                                                                                 ifelse(sample == "D03_S37","HDW-D03-S37",NA)))))))))))))))))  %>% 
  select(filename,gene_callers_id,START,END,sample,sample_names)
write_tsv(vfdb_genes,"vfdb_formatted_metacherchant.txt", col_names = F)

metacherchant_list <- card_CIA %>% left_join(gene_coverage_data) %>% 
  filter(genecoverage > 15)

write.csv(metacherchant_list, "metacherchant_list.csv", row.names = F)

card_CIA <- complete_genes %>% filter(source == "card") %>% 
  filter(RESISTANCE %in% c("aminoglycoside","fluoroquinolone","macrolide",
                           "tetracycline","cephamycin","glycopeptide","sulfonamide","cephalosporin")) %>% 
  mutate(filename = `#FILE`) %>% 
  separate(`#FILE`,c("sample","delete"), "-") %>% 
  mutate(sample_names = ifelse(sample == "C02_S104", "HDW-C02-S104",
                               ifelse(sample == "C02_S106","HDW-C02-S106",
                                      ifelse(sample == "C02_S111","HDW-C02-S111",
                                             ifelse(sample == "C02_S112", "HDW-C02-S112", 
                                                    ifelse(sample == "C02_S113","HDW-C02-S113",
                                                           ifelse(sample == "C02_S115","HDW-C02-S115",
                                                                  ifelse(sample == "D03_S25","HDW-D03-S25",
                                                                         ifelse(sample == "D03_S26","HDW-D03-S26",
                                                                                ifelse(sample == "D03_S27","HDW-D03-S27",
                                                                                       ifelse(sample == "D03_S28","HDW-D03-S28",
                                                                                              ifelse(sample == "D03_S31","HDW-D03-S31",
                                                                                                     ifelse(sample == "D03_S32","HDW-D03-S32",
                                                                                                            ifelse(sample == "D03_S33","HDW-D03-S33",
                                                                                                                   ifelse(sample == "D03_S34","HDW-D03-S34",
                                                                                                                          ifelse(sample == "D03_S35","HDW-D03-S35",
                                                                                                                                 ifelse(sample == "D03_S37","HDW-D03-S37",NA)))))))))))))))))  %>% 
  select(filename,gene_callers_id,START,END,sample,sample_names,RESISTANCE)

card_CIA$sample_names <- as.factor(card_CIA$sample_names)

card_nonCIA <- complete_genes %>% filter(source == "card") %>% 
  filter(!RESISTANCE %in% c("aminoglycoside","fluoroquinolone","macrolide",
                           "tetracycline","cephamycin","glycopeptide","sulfonamide","cephalosporin")) %>% 
  mutate(filename = `#FILE`) %>% 
  separate(`#FILE`,c("sample","delete"), "-") %>% 
  mutate(sample_names = ifelse(sample == "C02_S104", "HDW-C02-S104",
                               ifelse(sample == "C02_S106","HDW-C02-S106",
                                      ifelse(sample == "C02_S111","HDW-C02-S111",
                                             ifelse(sample == "C02_S112", "HDW-C02-S112", 
                                                    ifelse(sample == "C02_S113","HDW-C02-S113",
                                                           ifelse(sample == "C02_S115","HDW-C02-S115",
                                                                  ifelse(sample == "D03_S25","HDW-D03-S25",
                                                                         ifelse(sample == "D03_S26","HDW-D03-S26",
                                                                                ifelse(sample == "D03_S27","HDW-D03-S27",
                                                                                       ifelse(sample == "D03_S28","HDW-D03-S28",
                                                                                              ifelse(sample == "D03_S31","HDW-D03-S31",
                                                                                                     ifelse(sample == "D03_S32","HDW-D03-S32",
                                                                                                            ifelse(sample == "D03_S33","HDW-D03-S33",
                                                                                                                   ifelse(sample == "D03_S34","HDW-D03-S34",
                                                                                                                          ifelse(sample == "D03_S35","HDW-D03-S35",
                                                                                                                                 ifelse(sample == "D03_S37","HDW-D03-S37",NA)))))))))))))))))  %>% 
  select(filename,gene_callers_id,START,END,sample,sample_names)


write.csv(card, "all_card_genes.csv")

write_tsv(card_CIA, "card_CIA_formatted_metacherchant.txt")
write_csv(card_CIA, "card_CIA_formatted_metacherchant.csv")

write.csv(card_nonCIA,"card_nonCIA_formatted_metacherchant.csv", row.names = F)
write_tsv(card_nonCIA,"card_nonCIA_formatted_metacherchant.txt")

card_summary <- card %>% 
  group_by(`#FILE`,RESISTANCE) %>% 
  summarise(n = n())

names_resistance <- card_summary$RESISTANCE %>% unique() %>% unlist()

classes <- c("aminoglycoside","cephamycin","fluoroquinolone","macrolide","tetracycline",
             "cephalosporin","rifamycin","aminocoumarin","fosfomycin","glycopeptide",
             "nitroimidazole","peptide","phenicol","lincosamide","mupirocin",
             "diaminopyrimidine","sulfonamide","nucleoside")

reclassify_card <- card %>% 
                    mutate(abx_class = ifelse(RESISTANCE %in% classes, RESISTANCE, "multi-drug"))
card_summary <- reclassify_card %>% 
  group_by(`#FILE`,abx_class) %>% 
  summarise(n = n())

saveRDS(reclassify_card, "reclassify_card.RDS")
saveRDS(card_summary, "card_summary.RDS")

ggplot(card_summary, aes(x = `#FILE`, y = n, color = abx_class)) + 
  geom_jitter(width = 0.3, height = 0.3) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


ggplot(card_summary, aes(x = `#FILE`, y = n, color = abx_class)) +
  geom_bar(aes(fill = abx_class), position = "dodge", stat="identity") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# "C02_S104-gene-calls.fa" "C02_S106-gene-calls.fa" "C02_S111-gene-calls.fa"
# [4] "C02_S112-gene-calls.fa" "C02_S113-gene-calls.fa" "C02_S115-gene-calls.fa"
# [7] "D03_S25-gene-calls.fa"  "D03_S26-gene-calls.fa"  "D03_S27-gene-calls.fa" 
# [10] "D03_S28-gene-calls.fa"  "D03_S31-gene-calls.fa"  "D03_S32-gene-calls.fa" 
# [13] "D03_S33-gene-calls.fa"  "D03_S34-gene-calls.fa"  "D03_S35-gene-calls.fa" 
# [16] "D03_S37-gene-calls.fa" 
my_data_list <- list()

for (i in 1:length(names)){
 
 names[[i]] -> name

 complete_genes %>% 
    filter(complete_genes$`#FILE` == name) -> my_data
    
 my_data_list[[i]] <- my_data %>% select(-(`#FILE`))
    
 write_tsv(my_data_list[[i]],paste0(name,".txt"))
  }






