### metaphlan with radEmu

cmd_hmp2012 <- curatedMetagenomicData::sampleMetadata %>% 
  filter(study_name == "HMP_2012") |> 
  filter(body_site == "stool") |>
  select(where(~ !all(is.na(.x)))) |>
  returnSamples("relative_abundance", rownames = "short", counts = TRUE)

phyobj <- mia::makePhyloseqFromTreeSummarizedExperiment(cmd_hmp2012, assay.type = "relative_abundance", abund_values = "relative_abundance")

hmp2012_otu <- phyobj %>% phyloseq::otu_table()

hmp2012_tax <- phyobj %>% phyloseq::tax_table() %>% as.matrix() %>% as.data.frame()
#colSums(hmp2012_otu) 

hmp_2012_long <- hmp2012_otu %>% 
  as.matrix %>% as.data.frame() %>% 
  rownames_to_column("species") %>% 
  pivot_longer(cols=`SRS021948`:`SRS048164`,
               names_to="sample_id", values_to="abundance_counts")

hmp47_ids <- read.table("data/hmp47_ids.txt")
colnames(hmp47_ids) <- "sample_id"

# subset to the 85 people we care about
hmp47_long <- hmp_2012_long %>% right_join(hmp47_ids) %>% 
  mutate(sample = sample_id, 
         group = "HMP") %>% 
  select(-sample_id)
hmp47_long[is.na(hmp47_long)] = 0

#saveRDS(hmp47_long, "data/hmp47_long.RDS")

######## 
files <- list.files(path="data/results_w_read_stats", full.names = TRUE)
biobakery_counts <- data.table::rbindlist(sapply(files, fread,simplify = FALSE), idcol = 'filename')
biobakery_counts %>%
  filter(grepl("s__", clade_name)) %>%
  select(-clade_taxid, -clade_name, -coverage, -estimated_number_of_reads_from_the_clade) %>%
  group_by(filename) %>% 
  summarise(sum = sum(relative_abundance)) ## great sums to 100 
biobakery_counts <- readRDS("data/biobakery_counts.RDS")
all_wide_counts <- biobakery_counts %>%
  filter(grepl("s__", clade_name)) %>% 
  separate(filename, c("throw","tempID"), "profiled_HDW-") %>% 
  separate(tempID, c("sample","throw2"),".txt") %>% 
  select(clade_name,sample,estimated_number_of_reads_from_the_clade) %>% 
  pivot_wider(names_from = sample, values_from = estimated_number_of_reads_from_the_clade)

#saveRDS(all_wide_counts, "hdw_wide_counts_metaphlan3.RDS")
#all_wide_counts <- readRDS("hdw_wide_counts_metaphlan3.RDS")
all_wide_counts[is.na(all_wide_counts)] = 0
all_wide_counts <- all_wide_counts %>% remove_rownames %>% column_to_rownames(var="clade_name")
transpose_hdw_metaphlan_counts <- t(all_wide_counts)
#saveRDS(transpose_hdw_metaphlan_counts, "hdw_rademu_data.RDS") 

all_wide_counts_phylum <- biobakery_counts %>%
  filter(clade_name %in% c("k__Eukaryota|p__Eukaryota_unclassified",
                           "k__Bacteria|p__Verrucomicrobia",
                           "k__Bacteria|p__Spirochaetes", 
                           "k__Bacteria|p__Proteobacteria", 
                           "k__Bacteria|p__Firmicutes", 
                           "k__Bacteria|p__Bacteroidetes", 
                           "k__Archaea|p__Euryarchaeota", 
                           "k__Bacteria|p__Actinobacteria", 
                           "k__Bacteria|p__Synergistetes")) %>% 
  separate(filename, c("throw","tempID"), "profiled_HDW-") %>% 
  separate(tempID, c("sample","throw2"),".txt") %>% 
  select(clade_name,sample,estimated_number_of_reads_from_the_clade) %>% 
  pivot_wider(names_from = sample, values_from = estimated_number_of_reads_from_the_clade)
all_wide_counts_phylum[is.na(all_wide_counts_phylum)] = 0
all_wide_counts_phylum <- all_wide_counts_phylum %>% remove_rownames %>% column_to_rownames(var="clade_name")
transpose_hdw_metaphlan_counts_phylum <- t(all_wide_counts_phylum)
saveRDS(transpose_hdw_metaphlan_counts_phylum, "hdw_rademu_phylum.RDS")

hdw_otus_rademu <- readRDS("~/Documents/HDW/metagenomic/hdw_rademu_data.RDS")

hdw_meta_rademu <- readRDS("HMP_comparison/hdw_metadata_rademu.RDS") %>% 
  select(-hdw) %>%
  remove_rownames %>% 
  column_to_rownames(var="hdw_id")

hdw_fit <- radEmu::emuFit(formula = ~ dairy + age, 
                 data = hdw_meta_rademu,
                 Y = hdw_otus_rademu,
                 run_score_tests = FALSE) 
# fit takes about 2 minutes to run 

hdw_fit$B %>% rownames
covariate_to_test <- which("dairydairy" == hdw_fit$B %>% rownames)

robust_score <- emuFit(formula = ~ dairy + age,
                       data = hdw_meta_rademu,
                       fitted_model = hdw_fit,
                       refit = FALSE,
                       test_kj = data.frame(k = covariate_to_test, 
                                            j = 1), 
                       Y = as.matrix(hdw_otus_rademu),
                       rho_init = 1,
                       tau = 5,
                       constraint_tol = 1e-3)

ncores <- parallel::detectCores() - 1
ncores


emuTest <- function(category) {
  score_res <- emuFit(formula = ~ dairy + age,
                      data = hdw_meta_rademu,
                      fitted_model = hdw_fit,
                      refit = FALSE,
                      test_kj = data.frame(k = covariate_to_test, 
                                           j = category), 
                      Y = as.matrix(hdw_otus_rademu), 
                      rho_init = 1,
                      tau = 5,
                      constraint_tol = 1e-3)
  return(score_res)
}
library(parallel)
if (.Platform$OS.type != "windows") {
  # run if we are on a Mac or Linux machine
  score_res <- mclapply(1,
                        emuTest,
                        mc.cores = ncores)
} else {
  # don't run if we are on a Windows machine
  score_res <- NULL
}  

########################
# merge with hmp dataset 
########################
hdw_wide_counts <- readRDS("/Users/paulinetrinh/Documents/HDW/metagenomic/HMP_comparison/hdw_wide_counts_metaphlan3.RDS")
hdw_wide_counts[is.na(hdw_wide_counts)] = 0
hdw_format_species <- read_csv("/Users/paulinetrinh/Documents/HDW/metagenomic/HMP_comparison/hdw_format_species.csv")

hdw_long_counts <- hdw_wide_counts %>% 
  pivot_longer(cols =`C02-S104`: `D03-S37`, 
               names_to = "sample", values_to = "abundance_counts") %>% 
  mutate(group = ifelse(sample %in% c("C02-S104","C02-S106","C02-S111","C02-S112","C02-S113","C02-S115"),"community","dairy")) %>% 
  separate(clade_name, c("Leftover","Species"), sep = "s__", extra = "merge") %>%
  select(-Leftover) %>% 
  full_join(hdw_format_species, by = "Species") %>% 
  select(-c(Species))

full_cohort_long_counts <- hdw_long_counts %>% 
  full_join(hmp47_long)

hmp_comparison_dairy <- full_cohort_long_counts %>% 
  select(-group) %>% 
  pivot_wider(names_from = sample, values_from = abundance_counts)

hmp_comparison_wide_counts <- hmp_comparison_dairy %>% remove_rownames %>% column_to_rownames(var="species")
hmp_comparison_wide_counts[is.na(hmp_comparison_wide_counts)] = 0
remove_rowzeros = hmp_comparison_wide_counts[rowSums(hmp_comparison_wide_counts[])>0,]

transpose_hmp_dairy_metaphlan_counts <- t(remove_rowzeros)
saveRDS(transpose_hmp_dairy_metaphlan_counts, "2024_hmp_dairy_metaphlan_counts_otu_rademu.RDS")
transpose_hmp_dairy_metaphlan_counts <- readRDS("2024_hmp_dairy_metaphlan_counts_otu_rademu.RDS")
###########################################

my_data <- curatedMetagenomicData::sampleMetadata
hmp_done <- my_data %>% 
  filter(study_name == "HMP_2012")  %>% 
  filter(body_site == "stool")

hmp47_meta <- hmp_done %>% right_join(hmp47_ids) %>%
  mutate(group = "HMP", 
         hmp_study = 1)
hdw_metadata_clean <- readRDS("/Users/paulinetrinh/Documents/HDW/metagenomic/HMP_comparison/hdw_metadata_rademu.RDS") %>% 
  mutate(sample_id = hdw_id, 
         group = dairy, 
         hmp_study = 0)

hmp_hdw_rademu_meta <- hmp47_meta %>% 
  mutate(sex = gender) %>% 
  select(-gender) %>%
  full_join(hdw_metadata_clean) %>% 
  mutate(occupation = ifelse(group == "dairy", "dairy worker", 
                             ifelse(group == "community", "field worker", "unknown"))) %>% 
  remove_rownames %>% 
  column_to_rownames(var="sample_id")
saveRDS(hmp_hdw_rademu_meta,"hmp_hdw_rademu_meta.RDS")
hmp_hdw_rademu_meta <- readRDS("hmp_hdw_rademu_meta.RDS")

############################################
hmp_dairy_rademu_otus <- transpose_hmp_dairy_metaphlan_counts
saveRDS(hmp_dairy_rademu_otus, "hmp_dairy_rademu_otus.RDS")

hmp_dairy_rademu_meta <- hmp_hdw_rademu_meta %>%
  mutate( bmi_cat = ifelse(BMI <= 18.5 , "<18.5", 
                           ifelse(BMI > 18.5 & BMI <=25, "18.5-24.9", 
                                  ifelse(BMI > 25 & BMI <= 29.9, "25-29.9",
                                         ifelse(29.9 < BMI, "30.0+", NA))))) %>% 
  select(hmp_study, group, age, bmi_cat, BMI) 
saveRDS(hmp_dairy_rademu_meta,"hmp_dairy_rademu_meta.RDS")

hmp_dairy_rademu_otus <- readRDS("hmp_dairy_rademu_otus.RDS")
hmp_dairy_rademu_meta <- readRDS("hmp_dairy_rademu_meta.RDS")

hmp_hdw_fit <- radEmu::emuFit(formula = ~ group + age, 
                          data = hmp_dairy_rademu_meta,
                          tolerance = 0.01,
                          Y = hmp_dairy_rademu_otus,
                          run_score_tests = FALSE) 
# fit takes a few minutes to run 

covariate_to_test <- rbind(which("groupdairy" == hmp_hdw_fit$B %>% rownames),
                           which("groupHMP" == hmp_hdw_fit$B %>% rownames))

emuTest_hmp <- function(category) {
  
  score_res <- emuFit(formula = ~ group + age,
                      data = hmp_dairy_rademu_meta,
                      fitted_model = hmp_hdw_fit,
                      refit = FALSE,
                      test_kj = data.frame(k = 2, 
                                           j = 1), 
                      Y = as.matrix(hmp_dairy_rademu_otus), 
                      rho_init = 1,
                      tau = 5,
                      constraint_tol = 1e-3)
  
  return(score_res)
}
library(parallel)
if (.Platform$OS.type != "windows") {
  # run if we are on a Mac or Linux machine
  score_res <- mclapply(1:412,
                        emuTest_hmp,
                        mc.cores = ncores)
} else {
  # don't run if we are on a Windows machine
  score_res <- NULL
}  


#################################
## Read in radEmu results files #
#################################

hdw_results <- readRDS("/Users/paulinetrinh/Documents/HDW/metagenomic/rademu_results/hdw_results_fixed.rds")
hdw_results2 <- hdw_results %>% mutate(fold_change = exp(estimate)) %>% 
  filter(covariate == "dairydairy") %>%
  select(category, estimate)
hdw_format_species <- read_csv("/Users/paulinetrinh/Documents/HDW/metagenomic/HMP_comparison/hdw_format_species.csv")

clr_results <- readRDS("clr_tests_counts.RDS") %>% 
  separate(clade_name, c("Leftover","Species"), sep = "s__", extra = "merge") %>% 
  full_join(hdw_format_species) %>% 
  rename(category = species)

compare_results <- hdw_results2 %>% 
  full_join(clr_results)

cor_coefs <- cor.test(compare_results$estimate, compare_results$mean_diff)

ggplot(data = compare_results, aes(x = estimate, y = mean_diff)) + 
  geom_point()

combined_results <- readRDS("/Users/paulinetrinh/Documents/HDW/metagenomic/rademu_results/combined_results.rds")

combined_results2 <- combined_results %>% 
  filter(covariate == "groupdairy") %>%
  select(category, estimate.x, estimate.y)

compare_results <- combined_results2 %>% 
  full_join(clr_results)

ggplot(data = compare_results, aes(x = estimate.x, y = dairy_beta)) + 
  geom_point()

comparison_table <- readRDS("comparison_table.RDS")


hmp_results <- readRDS("/Users/paulinetrinh/Downloads/paulines_rademu_analysis_share/hmp_results.rds")




