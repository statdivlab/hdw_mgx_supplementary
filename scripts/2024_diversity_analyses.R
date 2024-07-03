# Diversity analyses 
## Bray-Curtis comparison 
##### Create a phyloseq obj 
############################
hdw_wide_counts <- readRDS("data/hdw_wide_counts_metaphlan3.RDS")
hdw_metadata_clean <- readRDS("data/HMP_comparison/hdw_metadata_rademu.RDS") %>% 
  mutate(sample_id = hdw_id, 
         group = dairy, 
         hmp_study = 0)

hdw_wide_counts_clean <- hdw_wide_counts %>% 
  separate(clade_name, c("Leftover","Species"), sep = "s__", extra = "merge") %>% 
  select(-Leftover) %>% 
  remove_rownames %>% 
  column_to_rownames(var = "Species")

hdw_meta_phy <- hdw_metadata_clean %>% 
  select(hdw_id, group, age) %>% 
  remove_rownames %>% 
  column_to_rownames(var = "hdw_id") %>% 
  mutate(group2 = as.factor(group)) %>% 
  select(-group)

hdw_tax_table <- hdw_wide_counts %>% 
  select(clade_name) %>% 
  separate(clade_name, c("Leftover","my_rownames"), sep = "s__", extra = "merge") %>%
  select(-Leftover) 

tax_levels <- strsplit(hdw_wide_counts$clade_name, "\\|")
hdw_tax_table$Kingdom <- sapply(tax_levels, function(x) x[1])
hdw_tax_table$Phylum <- sapply(tax_levels, function(x) x[2])
hdw_tax_table$Class <- sapply(tax_levels, function(x) x[3])
hdw_tax_table$Order <- sapply(tax_levels, function(x) x[4])
hdw_tax_table$Family <- sapply(tax_levels, function(x) x[5])
hdw_tax_table$Genus <- sapply(tax_levels, function(x) x[6])
hdw_tax_table$Species <- sapply(tax_levels, function(x) x[7])

hdw_tax_table2 <- hdw_tax_table %>% 
  remove_rownames %>% 
  column_to_rownames(var = "my_rownames") %>% 
  as.matrix()

my_tax <- tax_table(hdw_tax_table2)
my_OTU <- otu_table(hdw_wide_counts_clean, taxa_are_rows = TRUE)
my_meta <- sample_data(hdw_meta_phy)

my_hdw_mgx_physeq <- phyloseq(my_OTU, my_tax, my_meta)
saveRDS(my_hdw_mgx_physeq, "data/my_hdw_mgx_physeq.RDS")
hdw_species <- my_hdw_mgx_physeq %>%
  phyloseq::tax_glom(taxrank="Species")

hdw_phylum <- my_hdw_mgx_physeq %>%
  phyloseq::tax_glom(taxrank="Phylum")

#richness_hdw <- metaphlan_obj %>% breakaway

meta <- my_hdw_mgx_physeq %>%
  sample_data %>%
  as_tibble %>%
  mutate("sample_names" = my_hdw_mgx_physeq %>% sample_names )


set.seed(20200318)
divnet_species <- hdw_species %>%
  divnet(ncores = 4, tuning = "careful")
saveRDS(divnet_species, "data/divnet_species.RDS")

# X <- diag(nrow(sample_data(my_hdw_mgx_physeq)))
# colnames(X) <- rownames(X) <-  as.character(1:ncol(X))
# divnet_species_bc <-  divnet(W = my_hdw_mgx_physeq,
#                              X = X,
#                              variance = "none")
# saveRDS(divnet_species_bc, "divnet_species_bc.RDS")
set.seed(454)
bc_test <- testBetaDiversity(divnet_species, 
                             h0 = "bray-curtis",
                             n_boot = 10000,
                             sample_specimen_matrix = divnet_species$X,
                             groups = sample_data(my_hdw_mgx_physeq)$group2)
bc_test$p_value
# 0.399
saveRDS(bc_test, "data/bc_test.RDS")

library(magrittr)
estimates <- divnet_species$shannon %>% summary %$% estimate
ses <- sqrt(divnet_species$`shannon-variance`)
#X <- breakaway::make_design_matrix(divnet_species, "group2")
X <- model.matrix(~group2 + age, data = meta)

betta(estimates, ses, X)$table
# Estimates Standard Errors p-values
# (Intercept)  2.9318201514     0.129153913    0.000
# group2dairy -0.1904043460     0.163368213    0.244
# age          0.0008200099     0.002953343    0.781

