library(magrittr)
library(tidyverse)

# hdw results 
hdw_ests <- readRDS("hdw_fit.RDS")
hdw_coefs <- hdw_ests$coef[, 1:7]

# replace category names 
hdw_wide_counts <- readRDS("hdw_wide_counts_metaphlan3.RDS")
hdw_wide_counts[is.na(hdw_wide_counts)] <- 0
hdw_format_species <- read_csv("hdw_format_species.csv")

hdw_new_counts <- hdw_wide_counts %>% 
  pivot_longer(cols =`C02-S104`: `D03-S37`, 
               names_to = "sample", values_to = "abundance_counts") %>% 
  mutate(group = ifelse(sample %in% c("C02-S104","C02-S106","C02-S111","C02-S112","C02-S113","C02-S115"),"community","dairy")) %>% 
  separate(clade_name, c("Leftover","Species"), sep = "s__", extra = "merge") %>%
  select(-Leftover) %>% 
  full_join(hdw_format_species, by = "Species") %>%
  select(-Species, -group) %>%
  pivot_wider(names_from = species, values_from = abundance_counts) 
rownames(hdw_new_counts) <- hdw_new_counts$sample
hdw_new_counts <- hdw_new_counts %>% select(-sample)
hdw_coefs$category <- rep(colnames(hdw_new_counts), 2)
  
# read in score test results 
files <- list.files("results", full.names = T)
results <- vector(length(files), mode = "list")
counter <- 1
for (file in files) {
  res_so_far <- readRDS(file)
  get_datas <- file %>% strsplit("_") %>% extract2(1)  %>% strsplit(".rds")
  res_so_far$cat <- parse_number((file %>% strsplit("/") %>% extract2(1)  %>% strsplit(".rds"))[[2]][1])
  results[[counter]] <- res_so_far
  counter <- counter + 1
}
results <- do.call(rbind, results) %>% as_tibble %>%
  arrange(cat) %>%
  mutate(covariate = "dairydairy",
         category_num = cat) %>%
  select(-cat)

hdw_results <- left_join(hdw_coefs, results)
saveRDS(hdw_results, file = "hdw_results.rds")

# check on CI's for 7 cases in which score test didn't converge due to a computationally 
# singular matrix
hdw_results %>% filter(covariate == "dairydairy") %>%
  filter(is.na(score))
# unforatunely these are categories with big effect sizes and CI's that are far from 
# including 0. Sarah will try to run these manually and see where the error is coming from

# hmp results 
hmp_ests <- readRDS("hmp_fit_sm.RDS")
hmp_coefs <- hmp_ests$coef[, 1:7]

# read in score test results 
files <- list.files("results_bigger", full.names = T)
results <- vector(length(files), mode = "list")
counter <- 1
for (file in files) {
  res_so_far <- readRDS(file)
  get_datas <- file %>% strsplit("_") %>% extract2(1)  %>% strsplit(".rds")
  res_so_far$cat <- parse_number((file %>% strsplit("/") %>% extract2(1)  %>% strsplit(".rds"))[[2]][1])
  results[[counter]] <- res_so_far
  counter <- counter + 1
}
score_results <- do.call(rbind, results) %>% as_tibble %>%
  mutate(category_num = cat) %>%
  select(-cat, -contains("pval")) %>% 
  pivot_longer(cols = 2:3, names_to = "covariate", values_to = "score") %>%
  tidyr::separate(covariate, into = c("remove", "covariate")) %>%
  select(-remove) %>%
  mutate(covariate = paste0("group", covariate))
pval_results <- do.call(rbind, results) %>% as_tibble %>%
  mutate(category_num = cat) %>%
  select(-cat, -contains("score")) %>% 
  pivot_longer(cols = 2:3, names_to = "covariate", values_to = "pval") %>%
  tidyr::separate(covariate, into = c("remove", "covariate")) %>%
  select(-remove, -time) %>%
  mutate(covariate = paste0("group", covariate))

hmp_results <- left_join(hmp_coefs, score_results, by = c("category_num", "covariate")) %>%
  left_join(pval_results, by = c("category_num", "covariate"))
saveRDS(hmp_results, "hmp_results.rds")

# combine results from hdw and hmp studies 
combined_results <- hdw_results %>%
  mutate(covariate = ifelse(covariate == "dairydairy", "groupdairy", "age")) %>%
  full_join(hmp_results, by = c("covariate", "category"))
saveRDS(combined_results, file = "combined_results.rds")
ggplot(combined_results, aes(x = estimate.x, y = estimate.y)) + 
  geom_point()
