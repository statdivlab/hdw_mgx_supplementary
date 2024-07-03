library(tidyverse)

hdw_otus_rademu <- readRDS("hdw_rademu_data.RDS")

hdw_meta_rademu <- readRDS("hdw_metadata_rademu.RDS") %>% 
  remove_rownames %>% 
  column_to_rownames(var="hdw_id")

est_start <- proc.time()
hdw_fit <- radEmu::emuFit(formula = ~ dairy + age, 
                          data = hdw_meta_rademu,
                          Y = hdw_otus_rademu,
                          run_score_tests = FALSE,
                          verbose = TRUE) 
est_end <- proc.time() - est_start
saveRDS(hdw_fit, "hdw_fit.RDS")

est_start <- proc.time()
hdw_fit_longer <- radEmu::emuFit(formula = ~ dairy + age, 
                          data = hdw_meta_rademu,
                          Y = hdw_otus_rademu,
                          run_score_tests = FALSE,
                          maxit = 5000,
                          verbose = TRUE) 
est_end <- proc.time() - est_start

other_est_start <- proc.time()
other_hdw_fit <- radEmu::emuFit(formula = ~ dairy + age, 
                          data = hdw_meta_rademu,
                          tolerance = 0.01,
                          Y = hdw_otus_rademu,
                          run_score_tests = FALSE) 
other_est_end <- proc.time() - other_est_start
# fit takes about 2 minutes to run 

# run larger dataset 
hmp_otus_rademu <- readRDS("hmp_dairy_rademu_otus.RDS")
hmp_meta_rademu <- readRDS("hmp_dairy_radEmu_meta.RDS") 

est_start <- proc.time()
hmp_fit <- radEmu::emuFit(formula = ~ group + age, 
                          data = hmp_meta_rademu,
                          Y = hmp_otus_rademu,
                          run_score_tests = FALSE) 
est_end <- proc.time() - est_start
saveRDS(hmp_fit, "hmp_fit.RDS")

# run just over the 272 taxa
est_start <- proc.time()
hmp_fit_sm <- radEmu::emuFit(formula = ~ group + age, 
                          data = hmp_meta_rademu,
                          Y = hmp_otus_rademu[, 1:272],
                          run_score_tests = FALSE) 
est_end <- proc.time() - est_start
saveRDS(hmp_fit_sm, "hmp_fit_sm.RDS")

# run larger dataset with constraint only over hdw taxa 
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
hmp_in_hdw <- which(colnames(hmp_otus_rademu) %in% colnames(hdw_new_counts))
# set constraint over first 272 taxa 
constraint_fn_hmp <- (function(x) {
  radEmu:::pseudohuber_center(x[1:272], d = 0.1) })
constraint_grad_fn_hmp <- (function(x) {
  grad <- rep(0, length(x))
  grad[1:272] <- radEmu:::dpseudohuber_center_dx(x[1:272], d = 0.1)
  return(grad)
})
est_start <- proc.time()
hmp_fit_hdw_constraint <- radEmu::emuFit(formula = ~ group + age, 
                          data = hmp_meta_rademu,
                          Y = hmp_otus_rademu,
                          run_score_tests = FALSE,
                          constraint_fn = constraint_fn_hmp,
                          constraint_grad_fn = constraint_grad_fn_hmp,
                          constraint_param = NA,
                          verbose = TRUE) 
est_end <- proc.time() - est_start
