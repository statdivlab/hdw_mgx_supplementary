# happi analyses for CARD and VFDB identified genes 
library(tidyverse)
library(happi)
library(qvalue)

# CARD happi analysis
happi_card <- readRDS("happi_card.RDS")
x_matrix_AMR <- model.matrix(~dairy, data = happi_card)
run_happi_AMR <- function(colnum) {
  happi_results <-  happi(outcome=unlist(happi_card[,colnum]), 
                          covariate=x_matrix_AMR, 
                          quality_var=happi_card$reads_scaled,
                          method="splines", 
                          firth=TRUE, 
                          spline_df=3,
                          max_iterations=1000, 
                          change_threshold=0.1, 
                          epsilon=0.05, 
                          nstarts = 1, 
                          seed = 81)
  npLRT_results <- npLRT(happi_results, 
                         P = 1000, 
                         epsilon = 0.05, 
                         method="splines", 
                         firth=TRUE, 
                         spline_df=3,
                         max_iterations=1000, 
                         change_threshold=0.1)
  return(npLRT_results)
}

run_happi_AMR_beta <- function(colnum) {
  happi_results <-  happi(outcome=unlist(happi_card[,colnum]), 
                          covariate=x_matrix_AMR, 
                          quality_var=happi_card$reads_scaled,
                          method="splines", 
                          firth=TRUE, 
                          spline_df=3,
                          max_iterations=1000, 
                          change_threshold=0.1, 
                          epsilon=0.05, 
                          nstarts = 1, 
                          seed = 81)
  
  return(happi_results)
}

#run_happi_AMR_adjusted(3)
set.seed(8)
happi_amr_results <- parallel::mclapply(2:86, run_happi_AMR, mc.cores=6)
set.seed(8)
happi_amr_results_betas <- parallel::mclapply(2:86, run_happi_AMR_beta, mc.cores=6)

happi_amr_results_betas <- readRDS("happi_amr_results_betas.RDS")
happi_amr_results <- readRDS("happi_amr_results.RDS")
pvalue_amr <- happi_amr_results %>% unlist
beta_amr <- lapply(happi_amr_results_betas, function(x) tail(x$beta[!is.na(x$beta[,1]),],1)) %>% do.call("rbind",.)

hyp_results <- tibble("gene" = colnames(happi_card)[2:86], 
                         pvalue_amr, 
                         beta_amr) %>%
  arrange(pvalue_amr) %>% 
  mutate("qvalue" = qvalue::qvalue(pvalue_amr)$qvalue) 
write.csv(hyp_results,"hyp_results_qvalue_amr.csv")



## VFDB happi analysis 

happi_vfdb <- readRDS("happi_vfdb.RDS")

x_matrix_vfdb <- model.matrix(~dairy, data = happi_vfdb)
run_happi_vfdb_betas <- function(colnum) {
  happi(outcome=unlist(happi_vfdb[,colnum]), 
        covariate=x_matrix_vfdb, 
        quality_var=happi_vfdb$reads_scaled,
        method="splines", 
        firth=TRUE, 
        spline_df=3,
        max_iterations=1000, 
        change_threshold=0.1, 
        epsilon=0.05, 
        nstarts = 1, 
        seed = 81)
}

run_happi_vfdb <- function(colnum) {
  happi_results <-  happi(outcome=unlist(happi_vfdb[,colnum]), 
                          covariate=x_matrix_vfdb, 
                          quality_var=happi_vfdb$reads_scaled,
                          method="splines", 
                          firth=TRUE, 
                          spline_df=3,
                          max_iterations=1000, 
                          change_threshold=0.1, 
                          epsilon=0.05, 
                          nstarts = 1, 
                          seed = 81)
  npLRT_results <- npLRT(happi_results, 
                         P = 1000, 
                         epsilon = 0.05, 
                         method="splines", 
                         firth=TRUE, 
                         spline_df=3,
                         max_iterations=1000, 
                         change_threshold=0.1)
  return(npLRT_results)
  
}


set.seed(101)
happi_vfdb_results <- parallel::mclapply(2:38, run_happi_vfdb, mc.cores=6)
#saveRDS(happi_vfdb_results,"happi_vfdb_results.RDS")
#happi_vfdb_results <- readRDS("happi_vfdb_results.RDS")
set.seed(101)
happi_vfdb_betas <- parallel::mclapply(2:38, run_happi_vfdb_betas, mc.cores=6)
#saveRDS(happi_vfdb_betas,"happi_vfdb_betas.RDS")
#happi_vfdb_betas <- readRDS("happi_vfdb_betas.RDS")
pvalue_vf <- happi_vfdb_results %>% unlist
beta_vf <- lapply(happi_vfdb_betas, function(x) tail(x$beta[!is.na(x$beta[,1]),],1)) %>% do.call("rbind",.)
hyp_results_vfdb<- tibble("gene" = colnames(happi_vfdb)[2:38], 
                              pvalue_vf, 
                              beta_vf) %>%
  arrange(pvalue_vf) %>% 
  mutate("qvalue" = qvalue::qvalue(pvalue_vf)$qvalue) 
write.csv(hyp_results_vfdb,"hyp_results_vfdb_qvalue.csv")
