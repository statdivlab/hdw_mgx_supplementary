library(magrittr)
library(tidyverse)

### HDW results 

# read in estimation results
hdw_ests <- readRDS("results/hdw_fit.RDS")
hdw_coefs <- hdw_ests$coef[, 1:7]

# replace category names (from full taxonomy to species only)
hdw_wide_counts <- readRDS("data/hdw_wide_counts_metaphlan3.RDS")
hdw_wide_counts[is.na(hdw_wide_counts)] <- 0
hdw_format_species <- read_csv("data/hdw_format_species.csv")

hdw_new_counts <- hdw_wide_counts %>% 
  pivot_longer(cols =`C02-S104`: `D03-S37`, 
               names_to = "sample", values_to = "abundance_counts") %>% 
  mutate(group = ifelse(sample %in% c("C02-S104","C02-S106","C02-S111","C02-S112","C02-S113","C02-S115"),"community","dairy")) %>% 
  separate(clade_name, c("Leftover","Species"), sep = "s__", extra = "merge") %>%
  dplyr::select(-Leftover) %>% 
  full_join(hdw_format_species, by = "Species") %>%
  dplyr::select(-Species, -group) %>%
  pivot_wider(names_from = species, values_from = abundance_counts) 
rownames(hdw_new_counts) <- hdw_new_counts$sample
hdw_new_counts <- hdw_new_counts %>% dplyr::select(-sample)
hdw_coefs$category <- rep(colnames(hdw_new_counts), 2)

# read in score test results 
files <- list.files("score_results/hdw_only", full.names = T)
results <- vector(length(files), mode = "list")
counter <- 1
for (file in files) {
  res_so_far <- readRDS(file)
  res_so_far$cat <- parse_number((file %>% strsplit("/") %>% extract2(1)  %>% strsplit(".rds"))[[3]][1])
  results[[counter]] <- res_so_far
  counter <- counter + 1
}
results <- do.call(rbind, results) %>% as_tibble %>%
  arrange(cat) %>%
  mutate(covariate = "dairydairy",
         category_num = cat) %>%
  dplyr::select(-cat)

hdw_results <- left_join(hdw_coefs, results)
saveRDS(hdw_results, file = "results/hdw_results.rds")

# check how many tests converged
hdw_results %>%
  filter(covariate == "dairydairy") %>%
  group_by(converged) %>%
  count()
# 266 converged, 6 did not (due to computationally singular information matrix)

### HMP results 

# read in estimation results
hmp_ests <- readRDS("results/hmp_fit.RDS")
hmp_coefs <- hmp_ests$coef[, 1:7]

# read in score test results 
files <- list.files("score_results/hdw_and_hmp", full.names = T)
results <- vector(length(files), mode = "list")
counter <- 1
for (file in files) {
  res_so_far <- readRDS(file)
  res_so_far$cat <- parse_number((file %>% strsplit("/") %>% extract2(1)  %>% strsplit(".rds"))[[3]][1])
  results[[counter]] <- res_so_far
  counter <- counter + 1
}
results <- do.call(rbind, results) %>% as_tibble %>%
  arrange(cat) %>%
  mutate(covariate = "groupdairy",
         category_num = cat) %>%
  dplyr::select(-cat)

hmp_results <- left_join(hmp_coefs, results)
saveRDS(hmp_results, file = "results/hmp_results.rds")

# check how many tests converged
hmp_results %>%
  filter(covariate == "groupdairy") %>%
  group_by(convergence_dairy) %>%
  count()
hmp_results %>%
  filter(covariate == "groupdairy") %>%
  group_by(convergence_hmp) %>%
  count()

# all tests converged for the dairy group, 3 tests hit an iteration limit
# for the HMP group

# compare results 
full_results <- hdw_results %>%
  mutate(covariate = ifelse(covariate == "dairydairy", "groupdairy", "age")) %>%
  full_join(hmp_results, by = c("covariate", "category"))
ggplot(full_results, aes(x = estimate.x, y = estimate.y)) + 
  geom_point() + 
  labs(x = "HDW estimate", y = "HMP estimate")
ggplot(full_results, aes(x = pval, y = pval_dairy)) + 
  geom_point() + 
  labs(x = "HDW score test p-value", y = "HMP score test p-value")
ggsave("figures/pvalue_comparison.png")
cor(full_results$pval, full_results$pval_dairy, use = "complete.obs")

### HMP with BMI covariate results 

# read in estimation results
hmp_bmi_ests <- readRDS("results/hmp_fit_bmi.RDS")
hmp_bmi_coefs <- hmp_bmi_ests$coef[, 1:7]

# read in score test results 
files <- list.files("score_results/hdw_and_hmp_with_bmi/", full.names = T)
results <- vector(length(files), mode = "list")
counter <- 1
for (file in files) {
  res_so_far <- readRDS(file)
  res_so_far$cat <- parse_number((file %>% strsplit("/") %>% extract2(1)  %>% strsplit(".rds"))[[4]][1])
  results[[counter]] <- res_so_far
  counter <- counter + 1
}
results <- do.call(rbind, results) %>% as_tibble %>%
  arrange(cat) %>%
  mutate(covariate = "groupdairy",
         category_num = cat) %>%
  dplyr::select(-cat)

hmp_bmi_results <- left_join(hmp_bmi_coefs, results)
saveRDS(hmp_bmi_results, file = "results/hmp_bmi_results.rds")

# check how many tests converged
hmp_bmi_results %>%
  filter(covariate == "groupdairy") %>%
  group_by(convergence_dairy) %>%
  count()

# 6 tests failed for the dairy covariate, likely due to a computationally
# singular information matrix

# compare results 
combine_hmp_results <- hmp_results %>%
  full_join(hmp_bmi_results, by = c("covariate", "category"))
combine_hmp_results %>% 
  filter(covariate == "groupdairy") %>% 
  dplyr::select(contains("estimate")) %>%
  cor(use = "complete.obs")
combine_hmp_results %>% 
  filter(covariate == "groupdairy") %>% 
  dplyr::select(contains("pval_dairy")) %>%
  cor(use = "complete.obs")
combine_hmp_results %>% 
  filter(covariate == "groupdairy") %>% 
  ggplot(aes(x = estimate.x, y = estimate.y)) + 
  geom_point() + 
  labs(x = "HMP estimate", y = "HMP estimate (model with BMI)")
combine_hmp_results %>% 
  filter(covariate == "groupdairy") %>% 
  ggplot(aes(x = pval_dairy.x, y = pval_dairy.y)) + 
  geom_point() + 
  labs(x = "HMP score test p-val", y = "HMP score test p-val (model with BMI)")
cor(full_results$pval, full_results$pval_dairy, use = "complete.obs")
