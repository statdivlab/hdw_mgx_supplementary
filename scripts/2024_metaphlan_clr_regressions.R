## CLR + regressions. Compare with radEmu results 

clr <- function(x) {
  (x %>% log) - (x %>% log %>% mean)
}
biobakery_counts <- readRDS("data/biobakery_counts.RDS")
all_wide_counts <- biobakery_counts %>%
  filter(grepl("s__", clade_name)) %>% 
  separate(filename, c("throw","tempID"), "profiled_HDW-") %>% 
  separate(tempID, c("sample","throw2"),".txt") %>% 
  select(clade_name,sample,estimated_number_of_reads_from_the_clade) %>% 
  pivot_wider(names_from = sample, values_from = estimated_number_of_reads_from_the_clade)

all_wide_counts[is.na(all_wide_counts)] = 1
all_long_counts <- all_wide_counts %>% 
  pivot_longer(cols =`C02-S104`: `D03-S37`, 
               names_to = "sample", values_to = "abundance_counts") %>% 
  mutate(group = ifelse(sample %in% c("C02-S104","C02-S106","C02-S111","C02-S112","C02-S113","C02-S115"),"community","dairy"))

#all_long_counts <- readRDS("all_long_counts.RDS")
hdw_format_species <- read_csv("data/hdw_format_species.csv")
#saveRDS(hdw_format_species,"hdw_format_species.RDS")

comparison_table <- all_long_counts %>% 
  group_by(sample) %>% 
  mutate(clr = clr(abundance_counts)) %>% 
  separate(clade_name, c("Leftover","Species"), sep = "s__", extra = "merge") %>%
  select(-Leftover) %>% 
  full_join(hdw_format_species, by = "Species") %>% 
  select(-c(Species, abundance_counts)) %>% 
  pivot_wider(names_from = species, values_from = clr)

hdw_clr_regressions <- readRDS("data/hdw_metadata_rademu.RDS") %>% 
  select(hdw_id, age) %>% 
  rename(sample = hdw_id) %>% 
  full_join(comparison_table)

fit_anova <- function(colnum){
  lm(unlist(hdw_clr_regressions[,colnum]) ~ group + age, data = hdw_clr_regressions)
}

metaphlan3_comparison_results <- parallel::mclapply(4:275, fit_anova, mc.cores=6)

intercept_pvalues <- lapply(metaphlan3_comparison_results, function(x) summary(x)$coefficients[10]) %>% unlist
intercept_beta <- lapply(metaphlan3_comparison_results, function(x) summary(x)$coefficients[1]) %>% unlist
intercept_se <- lapply(metaphlan3_comparison_results, function(x) summary(x)$coefficients[4]) %>% unlist
intercept_tvalue <- lapply(metaphlan3_comparison_results, function(x) summary(x)$coefficients[7]) %>% unlist

dairy_pvalues <- lapply(metaphlan3_comparison_results, function(x) summary(x)$coefficients[11]) %>% unlist
dairy_beta <- lapply(metaphlan3_comparison_results, function(x) summary(x)$coefficients[2]) %>% unlist
dairy_se <- lapply(metaphlan3_comparison_results, function(x) summary(x)$coefficients[5]) %>% unlist
dairy_tvalue <- lapply(metaphlan3_comparison_results, function(x) summary(x)$coefficients[8]) %>% unlist

age_beta <- lapply(metaphlan3_comparison_results, function(x) summary(x)$coefficients[3]) %>% unlist
age_se <- lapply(metaphlan3_comparison_results, function(x) summary(x)$coefficients[6]) %>% unlist
age_tvalue <- lapply(metaphlan3_comparison_results, function(x) summary(x)$coefficients[9]) %>% unlist
age_pvalues <- lapply(metaphlan3_comparison_results, function(x) summary(x)$coefficients[12]) %>% unlist

pvalue_results <- tibble("species" = colnames(hdw_clr_regressions)[4:275], 
                         intercept_beta,
                         intercept_se,
                         intercept_tvalue,
                         intercept_pvalues, 
                         dairy_beta,
                         dairy_se,
                         dairy_tvalue,
                         dairy_pvalues, 
                         age_beta, 
                         age_se, 
                         age_tvalue, 
                         age_pvalues) %>% 
  arrange(dairy_pvalues)

pvalue_results$qvalue_dairy <- qvalue::qvalue(pvalue_results$dairy_pvalues, pi0.method = "bootstrap")$qvalues
#saveRDS(pvalue_results, "data/2024_HDW_metaphlan_clr_tests_results_ageadjusted.RDS")
#write.csv(pvalue_results, "2024_HDW_metaphlan_clr_tests_results_ageadjusted.csv", row.names = F)


#### compare with radEmu results 






#################################
hmp47_long <-  readRDS("data/hmp47_long.RDS")
full_cohort_long_counts <- hdw_long_counts %>% 
  full_join(hmp47_long)

hmp_comparison_dairy <- full_cohort_long_counts %>% 
  select(-group) %>% 
  pivot_wider(names_from = sample, values_from = abundance_counts)

hmp_comparison_wide_counts <- hmp_comparison_dairy %>% remove_rownames %>% column_to_rownames(var="species")
hmp_comparison_wide_counts[is.na(hmp_comparison_wide_counts)] = 0
remove_rowzeros = hmp_comparison_wide_counts[rowSums(hmp_comparison_wide_counts[])>0,]

pseudocounts <- replace(remove_rowzeros, remove_rowzeros == 0, 1)

all_long_counts_hmp <- pseudocounts %>% as.data.frame() %>%
  rownames_to_column(var = "Species") %>% 
  pivot_longer(cols =`C02-S104`: `SRS048164`, 
               names_to = "sample", values_to = "abundance_counts")

comparison_table_hmp <- all_long_counts_hmp %>% 
  group_by(sample) %>% 
  mutate(clr = clr(abundance_counts)) %>% 
  select(-c(abundance_counts)) %>% 
  pivot_wider(names_from = Species, values_from = clr)


hmp_hdw_meta <- readRDS("data/hmp_hdw_rademu_meta.RDS") %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "sample") %>% 
  select(sample, age, BMI, group) %>% 
  mutate(bmi_cat = ifelse(BMI <= 18.5 , "<18.5", 
                          ifelse(BMI > 18.5 & BMI <=25, "18.5-24.9", 
                                 ifelse(BMI > 25 & BMI <= 29.9, "25-29.9",
                                        ifelse(29.9 < BMI, "30.0+", NA)))))

hmp_hdw_clr_regression <- comparison_table_hmp %>% 
  full_join(hmp_hdw_meta)
saveRDS(hmp_hdw_clr_regression, "data/hmp_hdw_clr_regression.RDS")

fit_anova_hmp <- function(colnum){
  lm(unlist(hmp_hdw_clr_regression[,colnum]) ~ group + age + bmi_cat, data = hmp_hdw_clr_regression)
}

metaphlan3_comparison_results_hmp <- parallel::mclapply(2:413, fit_anova_hmp, mc.cores=6)

hmp_intercept_pvalues <- lapply(metaphlan3_comparison_results_hmp, function(x) summary(x)$coefficients[10]) %>% unlist
hmp_intercept_beta <- lapply(metaphlan3_comparison_results_hmp, function(x) summary(x)$coefficients[1]) %>% unlist
hmp_intercept_se <- lapply(metaphlan3_comparison_results_hmp, function(x) summary(x)$coefficients[4]) %>% unlist
hmp_intercept_tvalue <- lapply(metaphlan3_comparison_results_hmp, function(x) summary(x)$coefficients[7]) %>% unlist

hmp_dairy_pvalues <- lapply(metaphlan3_comparison_results_hmp, function(x) summary(x)$coefficients[11]) %>% unlist
hmp_dairy_beta <- lapply(metaphlan3_comparison_results_hmp, function(x) summary(x)$coefficients[2]) %>% unlist
hmp_dairy_se <- lapply(metaphlan3_comparison_results_hmp, function(x) summary(x)$coefficients[5]) %>% unlist
hmp_dairy_tvalue <- lapply(metaphlan3_comparison_results_hmp, function(x) summary(x)$coefficients[8]) %>% unlist

hmp_age_beta <- lapply(metaphlan3_comparison_results_hmp, function(x) summary(x)$coefficients[3]) %>% unlist
hmp_age_se <- lapply(metaphlan3_comparison_results_hmp, function(x) summary(x)$coefficients[6]) %>% unlist
hmp_age_tvalue <- lapply(metaphlan3_comparison_results_hmp, function(x) summary(x)$coefficients[9]) %>% unlist
hmp_age_pvalues <- lapply(metaphlan3_comparison_results_hmp, function(x) summary(x)$coefficients[12]) %>% unlist

hmp_pvalue_results <- tibble("species" = colnames(hmp_hdw_clr_regression)[2:413], 
                             hmp_intercept_beta,
                             hmp_intercept_se,
                             hmp_intercept_tvalue,
                             hmp_intercept_pvalues, 
                             hmp_dairy_beta,
                             hmp_dairy_se,
                             hmp_dairy_tvalue,
                             hmp_dairy_pvalues, 
                             hmp_age_beta, 
                             hmp_age_se, 
                             hmp_age_tvalue, 
                             hmp_age_pvalues) %>% 
  arrange(hmp_dairy_pvalues)

pvalue_results$qvalue_dairy <- qvalue::qvalue(pvalue_results$dairy_pvalues, pi0.method = "bootstrap")$qvalues

