library(tidyverse)
library(radEmu)

#### estimation for HDW data at the phylum level 

# load data
hdw_phylum <- readRDS("data/hdw_rademu_phylum.RDS")
hdw_meta_rademu <- readRDS("data/hdw_metadata_rademu.RDS") %>% 
  remove_rownames %>% 
  column_to_rownames(var="hdw_id")

# re-order metadata so that row names match up between datasets 
hdw_meta_reordered <- hdw_meta_rademu[rownames(hdw_phylum), ]
all.equal(rownames(hdw_meta_reordered), rownames(hdw_phylum))

# run estimation algorithm and score tests
est_start <- proc.time()
hdw_phylum_fit <- radEmu::emuFit(formula = ~ dairy + age, 
                          data = hdw_meta_reordered,
                          Y = hdw_phylum,
                          run_score_tests = TRUE,
                          verbose = TRUE) 
est_end <- proc.time() - est_start
saveRDS(hdw_phylum_fit, "results/hdw_phylum_fit.RDS")
saveRDS(hdw_phylum_fit$coef, "results/hdw_phylum_results.RDS")
