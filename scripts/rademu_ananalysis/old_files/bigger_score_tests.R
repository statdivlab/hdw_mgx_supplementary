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
hmp_otus_rademu <- readRDS("hdw/hmp_dairy_rademu_otus.RDS")
hmp_meta_rademu <- readRDS("hdw/hmp_dairy_radEmu_meta.RDS") 

# read in estimation fit 
hmp_fit <- readRDS("hdw/hmp_fit.RDS")
covariate_to_test1 <- which("groupdairy" == hmp_fit$B %>% rownames)
covariate_to_test2 <- which("groupHMP" == hmp_fit$B %>% rownames)

# time and run robust score test 
score_start <- proc.time() 
robust_score <- emuFit(formula = ~ group + age,
                       data = hmp_meta_rademu,
                       fitted_model = hmp_fit,
                       refit = FALSE,
                       test_kj = data.frame(k = c(covariate_to_test1, covariate_to_test2), 
                                            j = batch), 
                       Y = as.matrix(hmp_otus_rademu),
                       rho_init = 1,
                       tau = 1.5,
                       inner_maxit = 500,
                       constraint_tol = 1e-3)
score_end <- proc.time() - score_start
ind1 <- which(robust_score$coef$covariate == "groupdairy" & 
                robust_score$coef$category_num == batch)
ind2 <- which(robust_score$coef$covariate == "groupHMP" & 
                robust_score$coef$category_num == batch)
return_df <- data.frame(time = score_end[3],
                        score_dairy = robust_score$coef$score_stat[ind1],
                        pval_dairy = robust_score$coef$pval[ind1],
                        score_hmp = robust_score$coef$score_stat[ind2],
                        pval_hmp = robust_score$coef$pval[ind2])
saveRDS(return_df, paste0("hdw/results_bigger/res", batch, ".rds"))