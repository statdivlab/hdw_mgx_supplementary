library(tidyverse)
library(radEmu)

#### estimation for HDW data only 

# load data
hdw_otus_rademu <- readRDS("data/hdw_rademu_data.RDS")
hdw_meta_rademu <- readRDS("data/hdw_metadata_rademu.RDS") %>% 
  remove_rownames %>% 
  column_to_rownames(var="hdw_id")

# re-order metadata so that row names match up between datasets 
hdw_meta_reordered <- hdw_meta_rademu[rownames(hdw_otus_rademu), ]
all.equal(rownames(hdw_meta_reordered), rownames(hdw_otus_rademu))

# run estimation algorithm
est_start <- proc.time()
hdw_fit <- radEmu::emuFit(formula = ~ dairy + age, 
                          data = hdw_meta_reordered,
                          Y = hdw_otus_rademu,
                          run_score_tests = FALSE,
                          verbose = TRUE) 
est_end <- proc.time() - est_start
saveRDS(hdw_fit, "results/hdw_fit.RDS")
# this takes about 2 minutes to run locally

### estimation for HDW and HMP data

# load in data
hmp_otus_rademu <- readRDS("data/hmp_dairy_rademu_otus.RDS")
hmp_meta_rademu <- readRDS("data/hmp_dairy_radEmu_meta.RDS") 

# reorder metadata to match order of samples in otu table
hmp_meta_reordered <- hmp_meta_rademu[rownames(hmp_otus_rademu), ]
all.equal(rownames(hmp_meta_reordered), rownames(hmp_otus_rademu))

# make otu names the same between hmp and hdw datasets
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
hmp_in_hdw <- which(colnames(hmp_otus_rademu) %in% colnames(hdw_new_counts))

# filter HMP dataset to only include otus that are also in HDW dataset
hmp_otus_filtered <- hmp_otus_rademu[, hmp_in_hdw]
saveRDS(hmp_otus_filtered, "data/hmp_otus_filtered.RDS")

# run estimation algorithm
est_start <- proc.time()
hmp_fit <- radEmu::emuFit(formula = ~ group + age, 
                          data = hmp_meta_reordered,
                          Y = hmp_otus_filtered,
                          run_score_tests = FALSE,
                          verbose = TRUE) 
est_end <- proc.time() - est_start
saveRDS(hmp_fit, "results/hmp_fit.RDS")
# this takes about 3 minutes to run locally

### compare estimates 

combined_est <- hdw_fit$coef %>%
  mutate(covariate = ifelse(covariate == "dairydairy", "groupdairy", covariate)) %>%
  full_join(hmp_fit$coef, by = c("category_num", "covariate"))
ggplot(combined_est %>% filter(covariate == "groupdairy"), 
       aes(x = estimate.x, y = estimate.y)) + 
  geom_point() + 
  labs(x = "HDW estimate", y = "HMP estimate") + 
  geom_abline(slope = 1, color = "red")
ggsave("figures/estimate_comparison.png")
combined_est %>% 
  filter(covariate == "groupdairy") %>%
  dplyr::select(contains("estimate")) %>%
  cor(use = "complete.obs")

# look at discrepancies between plots 
upper_left_ind <- combined_est %>%
  filter(covariate == "groupdairy") %>%
  filter(estimate.x < -20 & estimate.y > 0) %>%
  pull(category_num)
as.vector(colSums(hdw_otus_rademu[, upper_left_ind] > 0))
as.vector(colSums(hmp_otus_rademu[, upper_left_ind] > 0))

lower_ind <- combined_est %>%
  filter(covariate == "groupdairy") %>%
  filter((estimate.x - estimate.y) > 10) %>%
  pull(category_num)
as.vector(colSums(hdw_otus_rademu[, lower_ind] > 0))
as.vector(colSums(hmp_otus_rademu[, lower_ind] > 0))

single_non_zero <- which(colSums(hdw_otus_rademu > 0) == 1)
combined_est %>% 
  filter(covariate == "groupdairy") %>%
  mutate(one_non_zero = ifelse(category_num %in% single_non_zero, TRUE, FALSE)) %>%
ggplot(aes(x = estimate.x, y = estimate.y, color = one_non_zero)) + 
  geom_abline(slope = 1, color = "black") + 
  geom_point() + 
  labs(x = "HDW estimate", y = "HMP estimate",
       color = "one >0 HDW sample") 
ggsave("figures/estimate_comparison_one_nonzero.png")

separated <- which(
  colSums(hdw_otus_rademu[hdw_meta_reordered$dairy == "dairy", ] > 0) == 0 |
    colSums(hdw_otus_rademu[hdw_meta_reordered$dairy != "dairy", ] > 0) == 0)
combined_est %>% 
  filter(covariate == "groupdairy") %>%
  mutate(separated = ifelse(category_num %in% separated, TRUE, FALSE)) %>%
  ggplot(aes(x = estimate.x, y = estimate.y, color = separated)) + 
  geom_abline(slope = 1, color = "black") + 
  geom_point() + 
  labs(x = "HDW estimate", y = "HMP estimate",
       color = "separation") 
ggsave("figures/estimate_comparison_separation.png")

combined_est %>% 
  filter(covariate == "groupdairy") %>%
  mutate(separated = ifelse(category_num %in% separated, TRUE, FALSE)) %>%
  filter(separated == TRUE) %>%
  dplyr::select(contains("estimate")) %>%
  cor(use = "complete.obs")
combined_est %>% 
  filter(covariate == "groupdairy") %>%
  mutate(separated = ifelse(category_num %in% separated, TRUE, FALSE)) %>%
  filter(separated == FALSE) %>%
  dplyr::select(contains("estimate")) %>%
  cor(use = "complete.obs")

# compare age estimates
ggplot(combined_est %>% filter(covariate == "age"), 
       aes(x = estimate.x, y = estimate.y)) + 
  geom_point() + 
  labs(x = "HDW age estimate", y = "HMP age estimate") + 
  geom_abline(slope = 1, color = "red")
ggsave("figures/age_estimate.png")

combined_est %>% 
  filter(covariate == "age") %>% 
  dplyr::select(contains("estimate")) %>%
  cor(use = "complete.obs")

left_ind <- combined_est %>%
  filter(covariate == "age") %>%
  filter(estimate.x < -4 & estimate.y > -1) %>%
  pull(category_num)
as.vector(colSums(hdw_otus_rademu[, left_ind] > 0))
as.vector(colSums(hmp_otus_rademu[, left_ind] > 0))

### estimation for HDW and HMP datasets with BMI added 

# load in data
hmp_meta_rademu_bmi <- readRDS("data/hmp_meta_bmi.RDS") 

# reorder metadata to match order of samples in otu table
hmp_meta_bmi_reordered <- hmp_meta_rademu_bmi[rownames(hmp_otus_rademu), ]
all.equal(rownames(hmp_meta_bmi_reordered), rownames(hmp_otus_rademu))

# run estimation algorithm
est_start <- proc.time()
hmp_fit_bmi <- radEmu::emuFit(formula = ~ group + age + bmi_cat, 
                          data = hmp_meta_bmi_reordered,
                          Y = hmp_otus_filtered,
                          run_score_tests = FALSE,
                          verbose = TRUE) 
est_end <- proc.time() - est_start
saveRDS(hmp_fit_bmi, "results/hmp_fit_bmi.RDS")
# this takes about 5.5 minutes to run locally

# compare estimates with and without bmi 
cor(hmp_fit$coef %>% filter(covariate == "groupdairy") %>% pull(estimate),
    hmp_fit_bmi$coef %>% filter(covariate == "groupdairy") %>% pull(estimate))
cor(hmp_fit$coef %>% filter(covariate == "groupHMP") %>% pull(estimate),
    hmp_fit_bmi$coef %>% filter(covariate == "groupHMP") %>% pull(estimate))
cor(hmp_fit$coef %>% filter(covariate == "age") %>% pull(estimate),
    hmp_fit_bmi$coef %>% filter(covariate == "age") %>% pull(estimate))
