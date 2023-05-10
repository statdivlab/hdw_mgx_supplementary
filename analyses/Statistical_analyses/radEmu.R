# radEmu script 
# read in counts table taken from metaphlan
counts_table <- readRDS("wide_counts_table.RDS") %>% 
  t() %>% 
  as.matrix()

dim(counts_table)

meta <- read_tsv("/Users/paulinetrinh/Documents/GitHub/statdivlab/HDW/mgx/hdw_mgx_metadata.txt") %>% 
  separate(id, sep="HDW-", c("trash","id")) %>% 
  select(-trash) %>% 
  column_to_rownames("id") %>% 
  mutate(dairy = ifelse(group == "dairy",2,1)) %>% 
  as_tibble() 
meta$group <- factor(meta$group,
                     levels = c("community","dairy"))

design_mat <- cbind(rep(1, nrow(meta)), meta$group)
dim(counts_table)
# now, fit the radEmu model 
hdw_mod_fit <- radEmu:::emuFit(X = design_mat,
                  Y = counts_table,
                  tolerance = .1,
                  method = "FL",
                  verbose = TRUE, 
                  reweight = TRUE)
saveRDS(hdw_mod_fit, "hdw_mod_fit.RDS")
hdw_mod_fit <- readRDS("hdw_mod_fit.RDS")


# next, bootstrap to build confidence intervals 
hdw_fit_cis <- radEmu:::emuCI(emuMod = hdw_mod_fit,
                 nboot = 1000,
                 parallel = TRUE,
                 ncore = 8)
saveRDS(hdw_fit_cis,"hdw_fit_cis.RDS")
hdw_fit_cis <- readRDS("hdw_fit_cis.RDS")
  
  
# plot confidence intervals 
rademu_plot <- hdw_fit_cis %>% 
  filter(row == 2) %>% 
  mutate(species = colnames(counts_table)) %>% 
  ggplot(aes(x = outcome_index, y = estimate)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                width = 0.25) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  coord_cartesian(ylim = c(-0.0005,0.0005))

##

significant_outcomes_smaller <- hdw_fit_cis %>% 
  filter(estimate < 0 & upper < 0) 
dim(significant_outcomes_smaller)

significant_outcomes_larger <- hdw_fit_cis %>% 
  filter(estimate > 0 & lower > 0) 
dim(significant_outcomes_larger)

## plots 
# taxa with smaller ratios in dairy workers
rademu_plot_smaller <- significant_outcomes_smaller %>% 
  filter(row == 2) %>% 
  ggplot(aes(x = outcome_index, y = estimate)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                width = 0.25) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  coord_cartesian(ylim = c(-0.0005,0.0005))

## larger ratios in dairy workers 
rademu_plot_larger <- significant_outcomes_larger %>% 
  filter(row == 2) %>% 
  ggplot(aes(x = outcome_index, y = estimate)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                width = 0.25) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  coord_cartesian(ylim = c(-0.0005,0.0005))


