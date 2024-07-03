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
redo_ind <- c(33, 45, 46, 96, 107, 108, 109, 116, 136, 148, 149, 152, 155, 156, 158,
              160, 162, 165, 170, 174, 179, 181, 190, 191, 215, 226, 233)
cat <- redo_ind[batch]

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
                                            j = cat), 
                       Y = as.matrix(hdw_otus_rademu),
                       rho_init = 1,
                       tau = 2, # decrease based on not converging the first time
                       inner_maxit = 50, # increase based on not converging the first time
                       constraint_tol = 1e-3)
score_end <- proc.time() - score_start
ind <- which(!(is.na(robust_score$coef$score_stat)))
return_df <- data.frame(time = score_end[3],
                        score = robust_score$coef$score_stat[ind],
                        pval = robust_score$coef$pval[ind])
saveRDS(return_df, paste0("hdw/results/res", cat, ".rds"))