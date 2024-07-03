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
hmp_otus <- readRDS("hdw/data/hmp_otus_filtered.RDS")
hmp_meta_rademu_bmi <- readRDS("hdw/data/hmp_meta_bmi.RDS") 

# reorder metadata to match order of samples in otu table
hmp_meta_bmi_reordered <- hmp_meta_rademu_bmi[rownames(hmp_otus), ]
all.equal(rownames(hmp_meta_bmi_reordered), rownames(hmp_otus))

# read in estimation fit 
hmp_fit_bmi <- readRDS("hdw/data/hmp_fit_bmi.RDS")
covariate_to_test <- which("groupdairy" == hmp_fit_bmi$B %>% rownames)

# time and run robust score test 
score_start <- proc.time() 
robust_score <- emuFit(formula = ~ group + age + bmi_cat,
                       data = hmp_meta_bmi_reordered,
                       fitted_model = hmp_fit_bmi,
                       refit = FALSE,
                       test_kj = data.frame(k = covariate_to_test, 
                                            j = batch), 
                       Y = as.matrix(hmp_otus),
                       rho_init = 1,
                       tau = 2,
                       inner_maxit = 50,
                       constraint_tol = 1e-3)
score_end <- proc.time() - score_start
ind <- which(robust_score$coef$covariate == "groupdairy" & 
                robust_score$coef$category_num == batch)
return_df <- data.frame(time = score_end[3],
                        score_dairy = robust_score$coef$score_stat[ind],
                        pval_dairy = robust_score$coef$pval[ind],
                        convergence_dairy = robust_score$score_test_hyperparams$converged)
saveRDS(return_df, paste0("hdw/hmp_bmi_results/res", batch, ".rds"))