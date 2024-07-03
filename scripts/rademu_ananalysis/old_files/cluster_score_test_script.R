library(radEmu)
library(tidyverse)

# get batch as command line argument 
args <- commandArgs(trailingOnly = FALSE)
if (length(args) == 0) {
  batch <- 1
} else {
  arg <- args[length(args)]
  batch <- abs(readr::parse_number(arg))
}

print(batch)

# read in data
hdw_otus_rademu <- readRDS("hdw/hdw_rademu_data.RDS")
hdw_meta_rademu <- readRDS("hdw/hdw_metadata_rademu.RDS") %>% 
  remove_rownames %>% 
  column_to_rownames(var="hdw_id")

# read in estimation fit 
hdw_fit <- readRDS("hdw/hdw_fit.RDS")
covariate_to_test <- which("dairydairy" == hdw_fit$B %>% rownames)

# time and run robust score test 
score_start <- proc.time() 
robust_score <- emuFit(formula = ~ dairy + age,
                       data = hdw_meta_rademu,
                       fitted_model = hdw_fit,
                       refit = FALSE,
                       test_kj = data.frame(k = covariate_to_test, 
                                            j = batch), 
                       Y = as.matrix(hdw_otus_rademu),
                       rho_init = 1,
                       tau = 5,
                       constraint_tol = 1e-3)
score_end <- proc.time() - score_start
ind <- which(!(is.na(robust_score$coef$score_stat)))
return_df <- data.frame(time = score_end[3],
                        score = robust_score$coef$score_stat[ind],
                        pval = robust_score$coef$pval[ind])
saveRDS(return_df, paste0("hdw/results/res", batch, ".rds"))
                                                                   
                                                                   
                                                                   