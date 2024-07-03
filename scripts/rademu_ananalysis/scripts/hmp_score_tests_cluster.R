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
hmp_meta <- readRDS("hdw/data/hmp_dairy_radEmu_meta.RDS") 

# reorder metadata to match order of samples in otu table
hmp_meta_reordered <- hmp_meta[rownames(hmp_otus), ]
all.equal(rownames(hmp_meta_reordered), rownames(hmp_otus))

# read in estimation fit 
hmp_fit <- readRDS("hdw/data/hmp_fit.RDS")
covariate_to_test1 <- which("groupdairy" == hmp_fit$B %>% rownames)
covariate_to_test2 <- which("groupHMP" == hmp_fit$B %>% rownames)

# time and run robust score test 
score_start <- proc.time() 
robust_score <- emuFit(formula = ~ group + age,
                       data = hmp_meta_reordered,
                       fitted_model = hmp_fit,
                       refit = FALSE,
                       test_kj = data.frame(k = c(covariate_to_test1, covariate_to_test2), 
                                            j = batch), 
                       Y = as.matrix(hmp_otus),
                       rho_init = 1,
                       tau = 2,
                       inner_maxit = 50,
                       constraint_tol = 1e-3)
score_end <- proc.time() - score_start
ind1 <- which(robust_score$coef$covariate == "groupdairy" & 
                robust_score$coef$category_num == batch)
ind2 <- which(robust_score$coef$covariate == "groupHMP" & 
                robust_score$coef$category_num == batch)
return_df <- data.frame(time = score_end[3],
                        score_dairy = robust_score$coef$score_stat[ind1],
                        pval_dairy = robust_score$coef$pval[ind1],
                        convergence_dairy = robust_score$score_test_hyperparams$converged[1],
                        score_hmp = robust_score$coef$score_stat[ind2],
                        pval_hmp = robust_score$coef$pval[ind2],
                        convergence_hmp = robust_score$score_test_hyperparams$converged[2])
saveRDS(return_df, paste0("hdw/hmp_results/res", batch, ".rds"))