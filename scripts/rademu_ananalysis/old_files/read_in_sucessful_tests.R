library(magrittr)
library(tidyverse)

# for first run of analysis of hdw data 
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
  arrange(cat)

failed_test <- which(!(1:272 %in% results$cat))

# attempt one of the 7 tests that failed after additional time on the cluster
# read in data
hdw_otus_rademu <- readRDS("hdw_rademu_data.RDS")
hdw_meta_rademu <- readRDS("hdw_metadata_rademu.RDS") %>% 
  remove_rownames %>% 
  column_to_rownames(var="hdw_id")

# read in estimation fit 
hdw_fit <- readRDS("hdw_fit.RDS")
covariate_to_test <- which("dairydairy" == hdw_fit$B %>% rownames)

# time and run robust score test 
score_start <- proc.time() 
robust_score <- emuFit(formula = ~ dairy + age,
                       data = hdw_meta_rademu,
                       fitted_model = hdw_fit,
                       refit = FALSE,
                       test_kj = data.frame(k = covariate_to_test, 
                                            j = 148), 
                       Y = as.matrix(hdw_otus_rademu),
                       rho_init = 1,
                       tau = 5,
                       constraint_tol = 1e-3, 
                       verbose = TRUE)
score_end <- proc.time() - score_start

# for first run of analysis of hmp data 
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
results <- do.call(rbind, results) %>% as_tibble %>%
  arrange(cat)

passed_test <- results$cat
