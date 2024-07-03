## relative abundance of AMR genes run with regression instead of t-test 

rigr::regress
install.packages('rigr')
library(rigr)
library(tidyverse)

gene_coverage_data <- readRDS("data/gene_coverage_data.RDS")
gene_coverage_per_sample <- gene_coverage_data %>% 
  group_by(sample) %>%
  summarise(coverage_sum = sum(genecoverage))

complete_genes <- readRDS("data/complete_amr_genes.RDS")
dim(complete_genes)

complete_genes_cleaned <- complete_genes %>% 
  separate(`#FILE` ,c("sample","filename","trash"),sep = "-") %>% 
  select(-c(filename,trash))

pseudocount_decision <- complete_genes_cleaned %>% 
  mutate(DIFF = END-START) 
min(pseudocount_decision$DIFF) # smallest length is 130 sequences 
# 0.0000000001

card_hits <- complete_genes_cleaned %>% 
  filter(source == "card")

identify_card <- card_hits %>% full_join(gene_coverage_data, by = c("sample","gene_callers_id")) %>% 
  mutate(card_gene = ifelse(source == "card",1,0), 
         dairy = ifelse(sample %in% c("C02_S104","C02_S106","C02_S111","C02_S112","C02_S113","C02_S115"),0,1), 
         group = ifelse(sample %in% c("C02_S104","C02_S106","C02_S111","C02_S112","C02_S113","C02_S115"),"community","dairy")) 
 # filter(source == "card") #247 rows

reclassify_card <- readRDS("data/reclassify_card2.RDS")

abx_class_RA <- identify_card %>% 
  full_join(reclassify_card) %>%
  mutate(abx_class_edit = ifelse(is.na(card_gene), "non-abx", abx_class)) %>%
  group_by(sample,abx_class_edit) %>% 
  summarise(abundance_counts = sum(genecoverage)) %>% 
  full_join(gene_coverage_per_sample) %>% 
  rename(gene_callers_id = abx_class_edit) %>% 
  mutate(relativeabundance = abundance_counts/coverage_sum)

check_for_NAs <- abx_class_RA %>% 
  select(gene_callers_id, sample, abundance_counts) %>% 
  pivot_wider(names_from = gene_callers_id, values_from = abundance_counts)
check_for_NAs[is.na(check_for_NAs)] = 0.000000000000001
# lots of NAs so we're going to have to add a pseudocount that is very small since the 
# scale of gene coverage is so small 
clr <- function(x) {
  (x %>% log) - (x %>% log %>% mean)
}
abx_class_comparison_table <- check_for_NAs %>%
  pivot_longer(cols=`aminoglycoside`:`nucleoside`,
               names_to="abx_class", values_to="abundance") %>%
  group_by(sample) %>% 
  mutate(clr = clr(abundance)) %>% 
  select(-abundance) %>%
  pivot_wider(names_from = abx_class, values_from = clr)

hdw_clr_amr_meta <- readRDS("data/hdw_metadata_rademu.RDS") %>% 
  select(hdw_id, age) %>% 
  rename(sample = hdw_id) %>% 
  mutate(dairy = ifelse(sample %in% c("C02-S104","C02-S106","C02-S111","C02-S112","C02-S113","C02-S115"), 0,1))
hdw_clr_amr_meta$sample <- gsub("-", "_", hdw_clr_amr_meta$sample)

hdw_clr_amr_regressions <- hdw_clr_amr_meta %>% 
  full_join(abx_class_comparison_table)

fit_regress <- function(colnum){
  rigr::regress("mean", unlist(hdw_clr_amr_regressions[,colnum]) ~ dairy + age, data = hdw_clr_amr_regressions)
}

#regress("mean", aminoglycoside ~ dairy + age, data = hdw_clr_amr_regressions)
#library(rigr)

amr_comparison_results <- parallel::mclapply(4:22, fit_regress, mc.cores=6)

# testing robustness to pseudocount
dairy_pvalues2 <- lapply(amr_comparison_results, function(x) x$coefficients[20]) %>% unlist
cor(amr_dairy_pvalues, dairy_pvalues2)
# with 1e-14 0.99 pearson corr
# with 0.0001 0.98355 corr 

amr_dairy_pvalues <- lapply(amr_comparison_results, function(x) x$coefficients[20]) %>% unlist
amr_dairy_beta <- lapply(amr_comparison_results, function(x) x$coefficients[2]) %>% unlist
amr_dairy_se <- lapply(amr_comparison_results, function(x) x$coefficients[8]) %>% unlist
amr_dairy_tvalue <- lapply(amr_comparison_results, function(x) x$coefficients[17]) %>% unlist

amr_intercept_pvalues <- lapply(amr_comparison_results, function(x) x$coefficients[19]) %>% unlist
amr_intercept_beta <- lapply(amr_comparison_results, function(x) x$coefficients[1]) %>% unlist
amr_intercept_se <- lapply(amr_comparison_results, function(x) x$coefficients[7]) %>% unlist
amr_intercept_tvalue <- lapply(amr_comparison_results, function(x) x$coefficients[16]) %>% unlist

amr_age_pvalues <- lapply(amr_comparison_results, function(x) x$coefficients[21]) %>% unlist
amr_age_beta <- lapply(amr_comparison_results, function(x) x$coefficients[3]) %>% unlist
amr_age_se <- lapply(amr_comparison_results, function(x) x$coefficients[9]) %>% unlist
amr_age_tvalue <- lapply(amr_comparison_results, function(x) x$coefficients[18]) %>% unlist


amr_pvalue_results <- tibble("species" = colnames(hdw_clr_amr_regressions)[4:22], 
                             amr_intercept_beta,
                             amr_intercept_se,
                             amr_intercept_tvalue,
                             amr_intercept_pvalues, 
                             amr_dairy_beta,
                             amr_dairy_se,
                             amr_dairy_tvalue,
                             amr_dairy_pvalues, 
                             amr_age_beta, 
                             amr_age_se, 
                             amr_age_tvalue, 
                             amr_age_pvalues) %>% 
  arrange(amr_dairy_pvalues) %>% 
  filter(species %in% c("fluoroquinolone","aminoglycoside",
                        "macrolide","tetracycline",
                        "cephalosporin","cephamycin",
                        "glycopeptide","sulfonamide"))

amr_pvalue_results$amr_qvalue_dairy <- qvalue::qvalue_truncp(amr_pvalue_results$amr_dairy_pvalues, pi0.method = "bootstrap")$qvalues

write.csv(amr_pvalue_results,"data/2024_amr_abundance_pvalue_results.csv", row.names = F)

saveRDS(amr_pvalue_results, "data/2024_amr_pvalue_results.RDS")

