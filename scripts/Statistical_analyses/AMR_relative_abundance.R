# AMR relative abundances April 28, 2022 
library(tidyverse)
df_list <- list.files(path="/Users/paulinetrinh/Documents/HDW/metagenomic/manual-refine/contigs-profiles/genecoverages", full.names = TRUE) %>% 
  lapply(read_tsv)

reformatted_data <- list()
for (i in 1:length(df_list)){
  colnames(df_list[[i]])[2] -> sample_name
  df_list[[i]] %>% 
    mutate(sample = sample_name) %>% 
    rename(genecoverage = paste0(sample_name), 
           gene_callers_id = key) ->   reformatted_data[[i]]
}

gene_coverage_data <- do.call("rbind",reformatted_data) %>% as_tibble()

unique_names <- unique(gene_coverage_data$sample) %>% 
  as.data.frame() %>% 
  rename(ID = ".") %>% 
  mutate(dairy = ifelse(ID %in% c("C02_S104","C02_S106","C02_S111","C02_S112","C02_S113","C02_S115"),0,1), 
         group = ifelse(ID %in% c("C02_S104","C02_S106","C02_S111","C02_S112","C02_S113","C02_S115"),"community","dairy"))
  

gene_counts <- gene_coverage_data %>% 
  pivot_wider(names_from = sample, values_from = genecoverage) 

gene_counts[is.na(gene_counts)] = 0

my_counts <- otu_table(gene_counts[,-1], taxa_are_rows = TRUE)
# Pull out the gene ID to use as "taxa" names
round_df <- function(df, digits) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  
  df[,nums] <- round(df[,nums], digits = digits)
  
  (df)
}
counts_rounded <- round(my_counts, 0)
taxa_names(counts_rounded) <- gene_counts %>% pull(gene_callers_id)


taxa_names(my_counts) <- gene_counts %>% pull(gene_callers_id)

my_samp <- sample_data(unique_names[,-1])
# Pull out the sample ID to use as sample names
# Note that I have to use the original data, 
# as I removed this column from my_samp
sample_names(my_samp) <- unique_names %>% pull(ID)

my_gene_physeq <- phyloseq(my_counts, my_samp)
my_gene_physeq <- phyloseq(counts_rounded, my_samp)

set.seed(1)
my_output <- corncob::differentialTest(formula = ~ dairy,
                              formula_null = ~ 1,
                              phi.formula = ~ dairy,
                              phi.formula_null = ~ dairy,
                              data = my_gene_physeq, 
                              test = "Wald", boot = FALSE,
                              fdr_cutoff = 0.05)
saveRDS(my_output,"corncob_gene_DA.RDS")

set.seed(8)
my_output_nofdr <- corncob::differentialTest(formula = ~ dairy,
                                       formula_null = ~ 1,
                                       phi.formula = ~ dairy,
                                       phi.formula_null = ~ dairy,
                                       data = my_gene_physeq, 
                                       test = "Wald", boot = FALSE)
#saveRDS(my_output_nofdr,"my_output_nofdr.RDS")

saveRDS(my_output_nofdr,"my_output_nofdr_genecounts.RDS")
my_output$significant_taxa

## differentially abundant genes
fdr_genes <- my_output_nofdr$p_fdr %>% 
  as.data.frame() %>% 
  rename(fdr_pvalue = ".") %>% 
  filter(fdr_pvalue < 0.05)

unique_gene_callers <-   card_hits %>% 
  group_by(gene_callers_id,source,`function`,accession) %>% 
  summarise(sum = sum())

### need to think about ways to cluster this
gene_RA <- gene_coverage_data %>% 
  filter(genecoverage > 0) %>% 
  group_by(sample) %>% 
  mutate(genecoverage_sum = sum(genecoverage)) %>% 
  mutate(relativeabundance = genecoverage/genecoverage_sum,
         "clr" = clr(relativeabundance), 
         log10_abundance = log10(relativeabundance))

download.file("https://raw.githubusercontent.com/statdivlab/stamps2019/master/labs/data/count_data.csv", "count_data.csv")
count_data <- read_csv("./count_data.csv")
download.file("https://raw.githubusercontent.com/statdivlab/stamps2019/master/labs/data/sample_data.csv", "sample_data.csv")
sample_data <- read_csv("./sample_data.csv")

gene_RA %>% 
  group_by(sample) %>% 
  summarise(sum_RA = sum(relativeabundance)) # check that it sums to 1 



clr <- function(x) {
  (x %>% log) - (x %>% log %>% mean)
}

complete_genes <- readRDS("/Users/paulinetrinh/Documents/HDW/metagenomic/complete_amr_genes.RDS")
dim(complete_genes)

complete_genes %>% 
  separate(`#FILE` ,c("sample","filename"),sep = "-") %>% 
  select(-filename) -> complete_genes_cleaned
dim(complete_genes_cleaned)

unique(complete_genes_cleaned$source)

complete_genes_cleaned %>% 
  filter(source == "vfdb") -> vfdb_hits 

complete_genes_cleaned %>% 
  filter(source == "card") -> card_hits

complete_genes_cleaned %>% 
  filter(source == "argannot") -> arg_hits

card_hits %>% left_join(gene_RA, by = c("sample","gene_callers_id")) -> card_RA

identify_card <- card_hits %>% full_join(gene_RA, by = c("sample","gene_callers_id")) %>% 
  mutate(card_gene = ifelse(source == "card",1,0), 
         dairy = ifelse(sample %in% c("C02_S104","C02_S106","C02_S111","C02_S112","C02_S113","C02_S115"),0,1), 
         group = ifelse(sample %in% c("C02_S104","C02_S106","C02_S111","C02_S112","C02_S113","C02_S115"),"community","dairy")) %>%
  filter(source == "card")
# 251  14

### Goal: CLR comparisons of antibiotic classes 
## If we sum the abundance of genes in a certain abx_class
## then CLR transform them 
## then compared the differences by group 
abx_class_RA <- identify_card %>% 
  full_join(reclassify_card) %>%
  group_by(sample,abx_class) %>% 
  summarise(relativeabundance = sum(relativeabundance)) %>% 
  rename(gene_callers_id = abx_class) %>% 
  mutate(source = "card",
         card_gene = 1)

all_genes_RA <- card_hits %>% full_join(gene_RA, by = c("sample","gene_callers_id")) %>% 
  mutate(card_gene = ifelse(source == "card",1,0)) %>% 
  filter(is.na(card_gene))

all_genes_RA2 <-  all_genes_RA %>% # 3236914      12
  full_join(abx_class_RA) %>% 
  mutate(dairy = ifelse(sample %in% c("C02_S104","C02_S106","C02_S111","C02_S112","C02_S113","C02_S115"),0,1), 
         group = ifelse(sample %in% c("C02_S104","C02_S106","C02_S111","C02_S112","C02_S113","C02_S115"),"community","dairy")) %>% 
  mutate("abx_clr" = clr(relativeabundance)) 
#3237165      14
#3237027      14
all_genes_RA2 %>% group_by(sample) %>% 
  summarise(sum_RA = sum(relativeabundance)) # check that it sums to 1 

ttests_genes <- all_genes_RA2 %>% 
group_by(`function`, group) %>%
  summarise(mean = mean(abx_clr), sd = sd(abx_clr), n = n()) %>%
  filter(n > 1) %>% 
  pivot_wider(names_from=group, values_from=mean:n) %>%
  filter(n_dairy > 1, n_community > 1) %>%
  mutate(t = (mean_community - mean_dairy)/sqrt(sd_community^2/n_community + sd_dairy^2/n_dairy), 
         df = n_community + n_dairy - 2) %>%
  arrange(desc(abs(t)))

abx_class_ttests <- ttests_genes %>% 
  filter(gene_callers_id %in% c("aminoglycoside","fluoroquinolone","macrolide",
                                "tetracycline","cephamycin","glycopeptide","sulfonamide","cephalosporin")) %>% 
  mutate(pvalue = pt(t,df,lower.tail = TRUE))
  

pt(q, df, lower.tail = TRUE)

sample <- c("C02_S104",
            "C02_S106",
            "C02_S111",
            "C02_S112",
            "C02_S113",
            "C02_S115",
            "D03_S25",
            "D03_S26",
            "D03_S27",
            "D03_S28",
            "D03_S31",
            "D03_S32",
            "D03_S33",
            "D03_S34",
            "D03_S35",
            "D03_S37")
sequencedreads <- c(26513204,
                    30480580,
                    29942131,
                    25262104,
                    28375384,
                    22302102,
                    23993557,
                    19736275,
                    23637209,
                    19562229,
                    23989230,
                    22233110,
                    22129067,
                    18687775,
                    20513598,
                    33675630)

sequence_reads_data <- data.frame(sample,sequencedreads)
community_t <- genome_diversity_plot %>% 
  filter(group == "community") %>% 
  select(estimated_genomes)

dairy_t <- genome_diversity_plot %>% 
  filter(group == "dairy") %>%
  select(estimated_genomes)
t.test(community_t, dairy_t, paired = F)

genome_numbers <- read_tsv("/Users/paulinetrinh/Documents/HDW/metagenomic/genome_diversity_results_compiled.txt")

genome_diversity_plot <- genome_numbers %>%
  left_join(sequence_reads_data, by = "sample") %>%
  mutate(group = ifelse(sample %in% c("C02_S104","C02_S106","C02_S111","C02_S112","C02_S113","C02_S115"),"community","dairy"))
  
amr_group_plot <- identify_card %>% 
  filter(source == "card") %>% 
  group_by(sample,group) %>% 
  summarise(n = n()) %>% 
  left_join(sequence_reads_data, by = "sample")

happi_card <- identify_card %>% 
  filter(source == "card") %>% 
  select(sample,`function`) %>% 
  group_by(sample,`function`) %>% 
  summarise(n = n()) %>% 
  mutate(presence = 1) %>%
  select(-n) %>% 
  pivot_wider(names_from = `function`, values_from = presence) %>% 
  left_join(sequence_reads_data, by = "sample") %>% 
  mutate(reads_scaled = sequencedreads/(1000000)) %>% 
  select(-sequencedreads) %>% 
  mutate(dairy = ifelse(sample %in% c("C02_S104","C02_S106",
                     "C02_S111","C02_S112",
                     "C02_S113","C02_S115"), 0,1))

happi_card[is.na(happi_card)] = 0

library(happi)

x_matrix_AMR <- model.matrix(~dairy, data = happi_card)
run_happi_AMR <- function(colnum) {
  happi(outcome=unlist(happi_card[,colnum]), 
        covariate=x_matrix_AMR, 
        quality_var=happi_card$reads_scaled,
        method="splines", 
        firth=T, 
        spline_df=4,
        max_iterations=1000, 
        change_threshold=0.01, 
        epsilon=0)
}ge

happi_amr_results <- parallel::mclapply(2:87, run_happi_AMR, mc.cores=8)

#saveRDS(happi_amr_results, "happi_amr_results.RDS")
happi_amr_results <- readRDS("happi_amr_results.RDS")
pvalue_amr <- lapply(happi_amr_results, function(x) tail(x$loglik$pvalue[!is.na(x$loglik$pvalue)], 1)) %>% unlist
beta_amr <- lapply(happi_amr_results, function(x) tail(x$beta[!is.na(x$beta[,1]),],1)) %>% do.call("rbind",.)
set.seed(123)
x <- rnorm(50, mean = c(rep(0, 25), rep(3, 25)))
p <- 2*pnorm(sort(-abs(x)))

round(p, 3)
round(p.adjust(p), 3)
round(p.adjust(sort(pvalue_amr), "BH", n = 86),5)


amr_hyp_results <- tibble("gene" = colnames(happi_card)[2:87], 
                                  pvalue_amr,
                                  beta_amr[,1],
                                  beta_amr[,2]) %>% 
  arrange(pvalue_amr) %>% 
  mutate(qvalue = qvalue(amr_hyp_results$pvalue_amr, pi0.method = "bootstrap", lambda = seq(0,0.9,0.05))$qvalue)

write.csv(amr_hyp_results,"amr_hyp_results_Aug8.csv")


my_pvalues <- amr_hyp_results$pvalue_amr 
p.adjust(my_pvalues, "BY")
hist(amr_hyp_results$pvalue_amr)
qvalues_happi <- qvalue(amr_hyp_results$pvalue_amr, pi0.method = "smoother", lambda = seq(0,0.90,0.05))
hist(qvalues_happi)

plot(qvalues_happi)

hist(qvalues_happi$pvalue, nclass = 20)



library(qvalue)
data(hedenfalk)
pvalues <- hedenfalk$p
qobj <- qvalue(p = pvalues)
summary(qobj)
hist(qobj)
plot(qobj)
# change tuning parameters to previous defaults

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("sva")

plot(qvalues_happi)
# weird p-values should use bootstrap method? 
amr_hyp_results = amr_hyp_results[order(amr_hyp_results$pvalue_amr),]
amr_hyp_results$BH = 
  p.adjust(amr_hyp_results$pvalue_amr, 
           method = "BH")
n <- length(amr_hyp_results$pvalue_amr)
p <- amr_hyp_results$pvalue_amr
# Given a vector of p-values p of length n 
i = n:1  # The reverse rank order
o <- order(p, decreasing = TRUE)
ro <- order(o)
pmin(1, cummin(n/i * p[o]))[ro]
## okay so I can definitely recreate the output here 

# So basically 

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("qvalue")
library(qvalue)

pvalues <- amr_hyp_results$pvalue_amr
qobj <- qvalue(p = pvalues,lambda = 0)

amr_hyp_results$BY = 
  p.adjust(amr_hyp_results$pvalue_amr, 
           method = "BH")



saveRDS(amr_hyp_results,"amr_hyp_results.RDS")
write.csv(amr_hyp_results,"amr_hyp_results.csv")

generichness <- read_csv("/Users/paulinetrinh/Documents/HDW/geneshot/estimated_richness.csv") %>% 
  separate(specimen, c("study","group_short","ID"),"-")
genome_diversity_generichness <- genome_diversity_plot %>% 
  separate(sample, c("group_short","ID"),"_") %>% 
  full_join(generichness, by = c("group_short","ID"))

png("group_genome_generichness_diversity_scatter.png", units="in", width=6, height=4, res=300)  
ggplot(genome_diversity_generichness,aes(x = estimate, y = estimated_genomes, color = group)) + 
  geom_point(size = 5, alpha = 0.8, position = position_jitter(w = 0.15, h = 0)) + 
  theme_bw() + 
  ylab("Number of Genomes Estimated Using SCGs") + 
  xlab("Estimated Gene Richness") + 
  scale_color_manual(values = c("dairy" = "darkorange", "community" = "darkblue")) + 
  theme(legend.position="right")
dev.off()


png("group_genome_diversity_scatter.png", units="in", width=6, height=4, res=300)  
ggplot(genome_diversity_plot,aes(x = sequencedreads, y = estimated_genomes, color = group)) + 
  geom_point(size = 5, alpha = 0.8, position = position_jitter(w = 0.15, h = 0)) + 
  theme_bw() + 
  labs(x = expression(paste("Number of Paired Reads Sequenced ", (10^6))), y = "Number of Genomes Estimated Using SCGs") + 
  scale_color_manual(values = c("dairy" = "darkorange", "community" = "darkblue")) + 
  theme(legend.position="right")
dev.off()
genome_diversity_plot <- genome_diversity_plot %>% 
  mutate(reads_scaled = sequencedreads/1000000)

genomediversity_fig <- ggplot(genome_diversity_plot,aes(x = reads_scaled, y = estimated_genomes, color = group)) + 
  geom_point(size = 6, alpha = 0.8, position = position_jitter(w = 0.15, h = 0)) + 
  theme_bw() + 
  labs(x = expression(paste("No. of PE Reads Sequenced ", (10^6))), y = "No. of Genomes Estimated Using SCGs") + 
  scale_color_manual(values = c("dairy" = "darkorange", "community" = "darkblue")) + 
  theme(legend.position="right") + 
  xlim(0,35) 

amr_group_plot <- amr_group_plot %>% 
  mutate(reads_scaled = sequencedreads/1000000)
png("group_card_scatter.png", units="in", width=6, height=4, res=300)  
ggplot(amr_group_plot,aes(x = reads_scaled, y = n, color = group)) + 
  geom_point(size = 5, alpha = 0.8, position = position_jitter(w = 0.15, h = 0)) + 
  theme_bw() + 
  labs(x = expression(paste("Number of Paired Reads Sequenced ", (10^6))), y = "Number of CARD genes") + 
  scale_color_manual(values = c("dairy" = "darkorange", "community" = "darkblue")) + 
  theme(legend.position="right")
dev.off()

card_fig <- ggplot(amr_group_plot,aes(x = reads_scaled, y = n, color = group)) + 
  geom_point(size = 6, alpha = 0.8, position = position_jitter(w = 0.15, h = 0)) + 
  theme_bw() + 
  labs(x = expression(paste("No. of PE Reads Sequenced ", (10^6))), y = "No. of CARD genes") + 
  scale_color_manual(values = c("dairy" = "darkorange", "community" = "darkblue")) + 
  theme(legend.position="right") + 
  xlim(0,35)


duplicate_entry <- identify_card %>% 
  slice(213:214) 
# duplicate_entry$genecoverage[1] + duplicate_entry$genecoverage[2]
merge_entries <- c("D03_S32",NA, "card", "CfxA2", "AF118110.1:1037-71",0,165.0104,
                   3466522.24782629,4.760114e-05,1,1,"dairy")
names(merge_entries) <- colnames(duplicate_entry)

### Make supplementary table 2 
STable2 <- identify_card %>% 
  filter(source == "card") %>% 
  slice(-c(213:214)) %>% 
  rbind(merge_entries) %>% 
  mutate(gene_coverage = round(as.numeric(genecoverage), digits = 2)) %>% 
  select(sample,"function",accession, gene_coverage) %>% 
  pivot_wider(names_from = sample, values_from = gene_coverage)
write.csv(STable2,"STable2.csv")


identify_vfdb <- vfdb_hits %>% full_join(gene_RA, by = c("sample","gene_callers_id")) %>% 
  mutate(vfdb_gene = ifelse(source == "vfdb",1,0), 
         dairy = ifelse(sample %in% c("C02_S104","C02_S106","C02_S111","C02_S112","C02_S113","C02_S115"),0,1), 
         group = ifelse(sample %in% c("C02_S104","C02_S106","C02_S111","C02_S112","C02_S113","C02_S115"),"community","dairy")) %>% 
  filter(source == "vfdb")
#######

clr_AMR_sample <- gene_RA %>% 
  left_join(reclassify_card) %>% 
  select(-clr) %>% 
  mutate(card_genes = ifelse(is.na(source), gene_callers_id, source)) %>% 
  group_by(sample,card_genes) %>% 
  summarise(RA_of_card = sum(relativeabundance))

clr_tests_card_all <- clr_AMR_sample %>% 
  ungroup() %>% 
  group_by(sample) %>% 
  mutate(clr = clr(RA_of_card), 
         group = ifelse(sample %in% c("C02_S104","C02_S106",
                                      "C02_S111","C02_S112",
                                      "C02_S113","C02_S115"),"community","dairy")) %>% 
  #  filter(abx_class_groupings %in% c("aminoglycoside","fluoroquinolone","macrolide",
  #                                    "tetracycline","cephamycin","glycopeptide","sulfonamide","cephalosporin")) %>% 
  ungroup() %>% 
  group_by(card_genes,group) %>%
  summarise(mean = mean(clr), sd = sd(clr), n = n()) %>%
  filter(n > 2) %>% 
  pivot_wider(names_from=group, values_from=mean:n) %>% 
  mutate(t = (mean_community - mean_dairy)/sqrt(sd_community^2/n_community + sd_dairy^2/n_dairy), 
         df = n_community + n_dairy - 2,
         pvalue = 2*pt(-abs(t),df=df)) %>%
  arrange(desc(abs(t))) %>% 
  mutate("qvalue" = p.adjust(pvalue, "BH"))
saveRDS(clr_tests_card_all,"clr_tests_card_all.RDS")

clr_tests_card_all %>% filter(card_genes == "card") -> card_RA_differences


clr_AMR_sample %>% filter(card_genes == "card")  %>% 
  mutate(group = ifelse(sample %in% c("C02_S104","C02_S106",
                                      "C02_S111","C02_S112",
                                      "C02_S113","C02_S115"),"community","dairy")) -> only_card


## mean relative abundances by sample 
identify_vfdb %>% group_by(sample) %>% 
  summarise(RA_of_vfdb = sum(relativeabundance), n = n()) %>% 
  mutate(group = ifelse(sample %in% c("C02_S104","C02_S106",
                                      "C02_S111","C02_S112",
                                      "C02_S113","C02_S115"),"community","dairy")) %>% 
  ungroup() %>% 
  group_by(group) %>% 
  summarise(mean = mean(RA_of_vfdb), sd = sd(RA_of_vfdb)) -> summarise_means_vfdb

data_for_vf_clrs <- vfdb_hits %>% full_join(gene_RA, by = c("sample","gene_callers_id")) %>% 
  mutate(vfdb_gene = ifelse(source == "vfdb",1,0), 
         dairy = ifelse(sample %in% c("C02_S104","C02_S106","C02_S111","C02_S112","C02_S113","C02_S115"),0,1), 
         group = ifelse(sample %in% c("C02_S104","C02_S106","C02_S111","C02_S112","C02_S113","C02_S115"),"community","dairy")) %>% 
  mutate(vfdb_genes = ifelse(is.na(source), gene_callers_id, source)) %>% 
  group_by(sample,vfdb_genes) %>% 
  summarise(RA_of_vfdb = sum(relativeabundance))
  saveRDS(data_for_vf_clrs,"data_for_vf_clrs.RDS")
  
  clr_tests_vfdb_all <- data_for_vf_clrs %>% 
    ungroup() %>% 
    group_by(sample) %>% 
    mutate(clr = clr(RA_of_vfdb), 
           group = ifelse(sample %in% c("C02_S104","C02_S106",
                                        "C02_S111","C02_S112",
                                        "C02_S113","C02_S115"),"community","dairy")) %>% 
    #  filter(abx_class_groupings %in% c("aminoglycoside","fluoroquinolone","macrolide",
    #                                    "tetracycline","cephamycin","glycopeptide","sulfonamide","cephalosporin")) %>% 
    ungroup() %>% 
    group_by(vfdb_genes,group) %>%
    summarise(mean = mean(clr), sd = sd(clr), n = n()) %>%
    pivot_wider(names_from=group, values_from=mean:n) %>% 
    mutate(t = (mean_community - mean_dairy)/sqrt(sd_community^2/n_community + sd_dairy^2/n_dairy), 
           df = n_community + n_dairy - 2,
           pvalue = 2*pt(-abs(t),df=df)) %>%
    arrange(desc(abs(t))) %>% 
    mutate("qvalue" = p.adjust(pvalue, "BH"))
saveRDS(clr_tests_vfdb_all,"clr_tests_vfdb_all.RDS")
clr_tests_vfdb_all %>% 
  filter(vfdb_genes == "vfdb") -> clr_vfdb_results


## let's use a nonparametric test 
ranked_test <- data_for_vf_clrs %>% 
  ungroup() %>% 
  group_by(sample) %>% 
  mutate(clr = clr(RA_of_vfdb), 
         group = ifelse(sample %in% c("C02_S104","C02_S106",
                                      "C02_S111","C02_S112",
                                      "C02_S113","C02_S115"),"community","dairy")) %>% 
  filter(vfdb_genes=="vfdb")
  #  filter(abx_class_groupings %in% c("aminoglycoside","fluoroquinolone","macrolide",
  #                                    "tetracycline","cephamycin","glycopeptide","sulfonamide","cephalosporin")) %>% 
  
ranked_test_results_vfdb <- wilcox.test(clr~ group,
                   data = ranked_test,
                   exact = FALSE)
# p = 0.1489
vfdb_plot_relativeabundance <- ranked_test %>% 
  left_join(sequence_reads_data, by = "sample") %>% 
  mutate(reads_scaled = sequencedreads/1000000)

png("group_vfdb_scatter_relativeabundance.png", units="in", width=6, height=4, res=300)  
ggplot(vfdb_plot_relativeabundance,aes(x = reads_scaled, y = RA_of_vfdb, color = group)) + 
  geom_point(size = 5, alpha = 0.8, position = position_jitter(w = 0.15, h = 0)) + 
  theme_bw() + 
  labs(x = expression(paste("Number of Paired Reads Sequenced ", (10^6))), y = "Relative Abundance of VFDB genes") + 
  scale_color_manual(values = c("dairy" = "darkorange", "community" = "darkblue")) + 
  theme(legend.position="right")
dev.off()


  
vfdb_plot <- identify_vfdb %>% 
  group_by(sample,group) %>% 
  summarise(n = n()) %>% 
  left_join(sequence_reads_data, by = "sample")

vfdb_plot <- vfdb_plot %>% 
  mutate(reads_scaled = sequencedreads/1000000)
library(mmtable2)
library(gt)
library(tidyverse)

data_wrangle <- vfdb_plot %>% 
  mutate(sample_names = ifelse(sample == "C02_S104", "C1",
                               ifelse(sample == "C02_S106","C2",
                                      ifelse(sample == "C02_S111","C3",
                                             ifelse(sample == "C02_S112", "C4", 
                                                    ifelse(sample == "C02_S113","C5",
                                                           ifelse(sample == "C02_S115","C6",
                                                                  ifelse(sample == "D03_S25","D1",
                                                                         ifelse(sample == "D03_S26","D2",
                                                                                ifelse(sample == "D03_S27","D3",
                                                                                       ifelse(sample == "D03_S28","D4",
                                                                                              ifelse(sample == "D03_S31","D5",
                                                                                                     ifelse(sample == "D03_S32","D6",
                                                                                                            ifelse(sample == "D03_S33","D7",
                                                                                                                   ifelse(sample == "D03_S34","D8",
                                                                                                                          ifelse(sample == "D03_S35","D9",
                                                                                                                                 ifelse(sample == "D03_S37","D10",NA))))))))))))))))) %>% 

  pivot_longer(
    cols = c(sequencedreads,n), 
    names_to = "my_y_values", 
    values_to = "my_data"
  ) %>% 
  mutate(my_variables = ifelse(my_y_values == "sequencedreads","No. of PE Reads Sequenced","No. of VFDB genes")) %>% 
  arrange(my_variables) %>% 
  ungroup()
data_wrangle$my_variables  <- factor(data_wrangle$my_variables, 
                                             levels=c("No. of VFDB genes","No. of PE Reads Sequenced"))

library(mmtable2)
library(gt)
main_table <- data_wrangle %>%
  mmtable(cells = my_data, table_name = "Number of VFDB genes") +
  # Specify Headers
  header_top(sample_names) +
  header_left(my_variables) + 
  # Specify formatting
  header_format(sample_names, list(cell_text(transform = "uppercase"))) +
  table_format(
    locations = list(
      cells_body(rows = c(2, 3))
    ),
    style = list(
      cell_borders(sides = "top", color = "grey")
    )
  )

png("group_vfdb_scatter.png", units="in", width=6, height=4, res=300)  
ggplot(vfdb_plot,aes(x = reads_scaled, y = n, color = group)) + 
  geom_point(size = 5, alpha = 0.8, position = position_jitter(w = 0.15, h = 0)) + 
  theme_bw() + 
  labs(x = expression(paste("Number of Paired Reads Sequenced ", (10^6))), y = "Number of VFDB genes") + 
  scale_color_manual(values = c("dairy" = "darkorange", "community" = "darkblue")) + 
  theme(legend.position="right")
dev.off()

vfdb_fig <- ggplot(vfdb_plot,aes(x = reads_scaled, y = n, color = group)) + 
  geom_point(size = 5, alpha = 0.8, position = position_jitter(w = 0.15, h = 0)) + 
  theme_bw() + 
  labs(x = expression(paste("Number of Paired Reads Sequenced ", (10^6))), y = "Number of VFDB genes") + 
  scale_color_manual(values = c("dairy" = "darkorange", "community" = "darkblue")) + 
  theme(legend.position="right")


happi_vfdb <- identify_vfdb %>% 
  filter(source == "vfdb") %>% 
  select(sample,`function`) %>% 
  group_by(sample,`function`) %>% 
  summarise(n = n()) %>% 
  mutate(presence = 1) %>%
  select(-n) %>% 
  pivot_wider(names_from = `function`, values_from = presence) %>% 
  left_join(sequence_reads_data, by = "sample") %>% 
  mutate(reads_scaled = sequencedreads/(10000000)) %>% 
  select(-sequencedreads) %>% 
  mutate(dairy = ifelse(sample %in% c("C02_S104","C02_S106",
                                      "C02_S111","C02_S112",
                                      "C02_S113","C02_S115"), 0,1))

happi_vfdb[is.na(happi_vfdb)] = 0

x_matrix_vfdb <- model.matrix(~dairy, data = happi_vfdb)
run_happi_vfdb <- function(colnum) {
  happi(outcome=unlist(happi_vfdb[,colnum]), 
        covariate=x_matrix_vfdb, 
        quality_var=happi_vfdb$reads_scaled,
        method="splines", 
        firth=T, 
        spline_df=4,
        max_iterations=1000, 
        change_threshold=0.01, 
        epsilon=0)
}

happi_vfdb_results <- parallel::mclapply(2:46, run_happi_vfdb, mc.cores=8)

saveRDS(happi_vfdb_results, "happi_vfdb_results.RDS")
happi_vfdb_results <- readRDS("happi_vfdb_results.RDS")
pvalue_vfdb <- lapply(happi_vfdb_results, function(x) tail(x$loglik$pvalue[!is.na(x$loglik$pvalue)], 1)) %>% unlist
beta_vfdb <- lapply(happi_vfdb_results, function(x) tail(x$beta[!is.na(x$beta[,1]),],1)) %>% do.call("rbind",.)

vfdb_hyp_results <- tibble("gene" = colnames(happi_vfdb)[2:46], 
                          pvalue_vfdb,
                          beta_vfdb[,1],
                          beta_vfdb[,2]) %>%  
  arrange(pvalue_vfdb) %>% 
  mutate("qvalue" = qvalue(vfdb_hyp_results$pvalue_vfdb, pi0.method = "bootstrap")$qvalue) %>% 
  arrange(pvalue_vfdb)

saveRDS(vfdb_hyp_results,"vfdb_hyp_results_Aug8.RDS")
write.csv(vfdb_hyp_results,"vfdb_hyp_results_Aug8.csv")

vfdb_hyp_results <- readRDS("vfdb_hyp_results.RDS")
vfdb_qobj <- qvalue(vfdb_hyp_results$pvalue_vfdb,pi0.method = "smoother", lambda = seq(0,0.95,0.05))

qvalues_happi <- qvalue(amr_hyp_results$pvalue_amr, pi0.method = "bootstrap", lambda = seq(0,0.9,0.05))

#######
library(ggpubr)

vfdb_table_data <- data_wrangle %>% 
  select(sample_names,my_data,my_variables) %>% 
  pivot_wider(names_from = sample_names, values_from = my_data) 
vfdb_table_data$my_variables  <- factor(vfdb_table_data$my_variables, 
                                     levels=c("No. of VFDB genes","No. of PE Reads Sequenced"))

vfdb_table_data <- vfdb_table_data %>% arrange(my_variables) %>% 
  rename(" "= my_variables)

vfdb_table_data.p <- ggtexttable(vfdb_table_data, rows = NULL, 
                        theme = ttheme("classic") )  + theme(plot.margin=unit(c(1,1,-0.5,1),"cm"))

                        
main.title <- "C. Number of Virulence Factor Genes Identified"
withtitle <- vfdb_table_data.p %>% 
  tab_add_title(text = main.title, face = "bold", padding = unit(1, "line"))

pdf("genes_counts_figure.pdf", width = 9, height = 9)
ggarrange(ggarrange(card_fig, genomediversity_fig, ncol = 2, labels = c("A", "B"), common.legend=TRUE, 
                    legend = "bottom"), withtitle, 
          labels = c("",""), nrow = 2)
dev.off()

### Make supplementary table 3 
duplicate_entry_vfdb <- identify_vfdb %>% 
  slice(55,60) 
duplicate_entry_vfdb$genecoverage[1] + duplicate_entry_vfdb$genecoverage[2]
merge_entries_vfdb <- c("C02_S113 ",NA, "vfdb", "gmhA/lpcA", "NP_439337",0,17.80684,
                        5149179,3.45819e-06,1,0,"community")
names(merge_entries_vfdb) <- colnames(duplicate_entry_vfdb)

STable3 <- identify_vfdb %>% 
  slice(-c(55,60)) %>% 
  rbind(merge_entries_vfdb) %>% 
  mutate(gene_coverage = round(as.numeric(genecoverage), digits = 2)) %>% 
  select(sample,"function",accession, gene_coverage) %>% 
  pivot_wider(names_from = sample, values_from = gene_coverage)
write.csv(STable3,"STable3.csv")

###### Differential abundance testing between antibiotic classes June 21, 2022
# I have gene_RA which provides me with all the relative abundances by gene_call 
# need to aggregate gene abundances by antibiotic classes 
# then do CLR transformation on those 
# then perform differential abundance testing through t-tests 
clr_AMR_Class <- gene_RA %>% 
  left_join(reclassify_card) %>% 
  select(-clr) %>% 
  mutate(abx_class_groupings = ifelse(is.na(abx_class), gene_callers_id, abx_class)) %>% 
  group_by(sample,abx_class_groupings) %>% 
  summarise(RA_of_classes = sum(relativeabundance))


clr_AMR_Class %>% 
  ungroup() %>% 
  group_by(sample) %>% 
  summarise(sum_RA = sum(RA_of_classes)) # yay still sums to 1 whoooo 

clr_tests_amr <- clr_AMR_Class %>% 
  ungroup() %>% 
  group_by(sample) %>% 
  mutate(clr = clr(RA_of_classes), 
        group = ifelse(sample %in% c("C02_S104","C02_S106",
                                             "C02_S111","C02_S112",
                                             "C02_S113","C02_S115"),"community","dairy")) %>% 
  filter(abx_class_groupings %in% c("aminoglycoside","fluoroquinolone","macrolide",
                                 "tetracycline","cephamycin","glycopeptide","sulfonamide","cephalosporin")) %>% 
  ungroup() %>% 
  group_by(abx_class_groupings, group) %>%
  summarise(mean = mean(clr), sd = sd(clr), n = n()) %>%
  filter(n > 2) %>% 
  pivot_wider(names_from=group, values_from=mean:n) %>% 
  mutate(t = (mean_community - mean_dairy)/sqrt(sd_community^2/n_community + sd_dairy^2/n_dairy), 
         df = n_community + n_dairy - 2,
         pvalue = 2*pt(-abs(t),df=df)) %>%
  arrange(desc(abs(t))) %>% 
  rename(abx_class = abx_class_groupings)
saveRDS(clr_tests_amr,"clr_tests_amr.RDS")

gene_RA %>% 
  left_join(reclassify_card) %>% 
  select(-clr) %>% 
  mutate(abx_class_groupings = ifelse(is.na(abx_class), gene_callers_id, abx_class)) -> gene_clrs
gene_clrs %>% 
  ungroup() %>% 
  group_by(sample) %>% 
  summarise(sum_RA = sum(relativeabundance)) # great sums to 1 
#View(gene_clrs)
clr_tests_args <- gene_clrs %>% 
  ungroup() %>% 
  group_by(sample) %>% 
  mutate(clr = clr(relativeabundance), 
         group = ifelse(sample %in% c("C02_S104","C02_S106",
                                      "C02_S111","C02_S112",
                                      "C02_S113","C02_S115"),"community","dairy")) %>% 
#  filter(abx_class_groupings %in% c("aminoglycoside","fluoroquinolone","macrolide",
#                                    "tetracycline","cephamycin","glycopeptide","sulfonamide","cephalosporin")) %>% 
  ungroup() %>% 
  group_by(`function`, group) %>%
  summarise(mean = mean(clr), sd = sd(clr), n = n()) %>%
  filter(n > 2) %>% 
  pivot_wider(names_from=group, values_from=mean:n) %>% 
  mutate(t = (mean_community - mean_dairy)/sqrt(sd_community^2/n_community + sd_dairy^2/n_dairy), 
         df = n_community + n_dairy - 2,
         pvalue = 2*pt(-abs(t),df=df)) %>%
  arrange(desc(abs(t))) %>% 
  mutate("qvalue" = p.adjust(pvalue, "BH"))  

saveRDS(gene_clrs, "gene_clrs.RDS")
saveRDS(clr_tests_args, "clr_tests_args.RDS")

write.csv(clr_tests_args,"clr_tests_args.csv")
#### 
### abundance by sample 
clr_AMR_sample <- gene_RA %>% 
  left_join(reclassify_card) %>% 
  select(-clr) %>% 
  mutate(card_genes = ifelse(is.na(source), gene_callers_id, source)) %>% 
  group_by(sample,card_genes) %>% 
  summarise(RA_of_card = sum(relativeabundance))

clr_tests_card_all <- clr_AMR_sample %>% 
  ungroup() %>% 
  group_by(sample) %>% 
  mutate(clr = clr(RA_of_card), 
         group = ifelse(sample %in% c("C02_S104","C02_S106",
                                      "C02_S111","C02_S112",
                                      "C02_S113","C02_S115"),"community","dairy")) %>% 
  #  filter(abx_class_groupings %in% c("aminoglycoside","fluoroquinolone","macrolide",
  #                                    "tetracycline","cephamycin","glycopeptide","sulfonamide","cephalosporin")) %>% 
  ungroup() %>% 
  group_by(card_genes,group) %>%
  summarise(mean = mean(clr), sd = sd(clr), n = n()) %>%
  filter(n > 2) %>% 
  pivot_wider(names_from=group, values_from=mean:n) %>% 
  mutate(t = (mean_community - mean_dairy)/sqrt(sd_community^2/n_community + sd_dairy^2/n_dairy), 
         df = n_community + n_dairy - 2,
         pvalue = 2*pt(-abs(t),df=df)) %>%
  arrange(desc(abs(t))) %>% 
  mutate("qvalue" = p.adjust(pvalue, "BH"))
saveRDS(clr_tests_card_all,"clr_tests_card_all.RDS")

clr_tests_card_all %>% filter(card_genes == "card") -> card_RA_differences


clr_AMR_sample %>% filter(card_genes == "card")  %>% 
  mutate(group = ifelse(sample %in% c("C02_S104","C02_S106",
                                      "C02_S111","C02_S112",
                                      "C02_S113","C02_S115"),"community","dairy")) -> only_card


saveRDS(only_card,"card_abundances.RDS")
only_card %>% 
  group_by(group) %>% 
  summarise(mean = mean(RA_of_card), sd =  sd(RA_of_card), n = n())

amr_group_plot_relativeabundances <- only_card %>% 
  left_join(sequence_reads_data, by = "sample")

amr_group_plot_relativeabundances <- amr_group_plot_relativeabundances %>% 
  mutate(reads_scaled = sequencedreads/1000000)
png("group_card_scatter_relativeabundances.png", units="in", width=6, height=4, res=300)  
ggplot(amr_group_plot_relativeabundances,aes(x = reads_scaled, y = RA_of_card, color = group)) + 
  geom_point(size = 6, alpha = 0.8, position = position_jitter(w = 0.15, h = 0)) + 
  theme_bw() + 
  labs(x = expression(paste("No. of PE Reads Sequenced ", (10^6))), y = "Relative Abundance of CARD genes") + 
  scale_color_manual(values = c("dairy" = "darkorange", "community" = "darkblue")) + 
  theme(legend.position="right") + 
  xlim(0,35)
dev.off()






## making abundance heatmap 
reclassify_card <- readRDS("reclassify_card.RDS")

reclassify_card <- reclassify_card %>% 
  separate(`#FILE` ,c("sample","filename"),sep = "-") %>% 
  select(-filename)

colmn <- paste0("resistance", 1:12)

multidrug <- identify_card %>% 
  full_join(reclassify_card) %>%
  group_by(sample,abx_class) %>% 
  filter(abx_class == "multi-drug") %>% 
  tidyr::separate(data = .,
    col = RESISTANCE,
    sep = ";",
    into = colmn,
    remove = FALSE)

unique_multidrug_genes <- unique(multidrug$`function`)  

write.csv(multidrug,"multidruggenes.csv")
write.csv(unique_multidrug_genes,"unique_multidrug_genes.csv")

intermediate_heatmap <- identify_card %>% 
  full_join(reclassify_card) %>%
  mutate(sample_names = ifelse(sample == "C02_S104", "C1",
                               ifelse(sample == "C02_S106","C2",
                                      ifelse(sample == "C02_S111","C3",
                                             ifelse(sample == "C02_S112", "C4", 
                                                    ifelse(sample == "C02_S113","C5",
                                                           ifelse(sample == "C02_S115","C6",
                                                                  ifelse(sample == "D03_S25","D1",
                                                                         ifelse(sample == "D03_S26","D2",
                                                                                ifelse(sample == "D03_S27","D3",
                                                                                       ifelse(sample == "D03_S28","D4",
                                                                                              ifelse(sample == "D03_S31","D5",
                                                                                                     ifelse(sample == "D03_S32","D6",
                                                                                                            ifelse(sample == "D03_S33","D7",
                                                                                                                   ifelse(sample == "D03_S34","D8",
                                                                                                                          ifelse(sample == "D03_S35","D9",
                                                                                                                                 ifelse(sample == "D03_S37","D10",NA))))))))))))))))) 
intermediate_heatmap$sample_names <- as.character(intermediate_heatmap$sample_names)
intermediate_heatmap$sample_names  <- factor(intermediate_heatmap$sample_names , 
                                             levels=c("C1", "C2",  "C3",  "C4",  "C5" , "C6",  "D1" , "D2",  "D3",  "D4",  "D5" , "D6",  "D7",  "D8" ,
                                                                                          "D9", "D10" ))
card_heatmap <- intermediate_heatmap %>% 
  group_by(sample_names,abx_class) %>% 
  filter(abx_class %in% c("aminoglycoside","fluoroquinolone","macrolide",
                          "tetracycline","cephamycin","glycopeptide","sulfonamide","cephalosporin")) %>% 
  summarise(AbxClass_relativeabundance = sum(relativeabundance))  %>% 
  mutate(log10_abundance = log10(AbxClass_relativeabundance)) %>% 
  select(-AbxClass_relativeabundance) %>% 
  left_join(clr_tests_amr, by = "abx_class") %>% 
  mutate(abx_pvalue = paste0(abx_class," (p=",pvalue,")")) %>% 
  arrange(sample_names,pvalue) %>% 
  select(sample_names,log10_abundance,abx_pvalue) %>%
  pivot_wider(names_from = sample_names, values_from = log10_abundance) %>% 
  column_to_rownames(var = "abx_pvalue") %>%
  as.matrix()

clr_tests_amr$pvalue <- round(clr_tests_amr$pvalue,2)

card_heatmap_all <- intermediate_heatmap %>% 
  group_by(sample_names,abx_class) %>% 
  summarise(AbxClass_relativeabundance = sum(relativeabundance))  %>% 
  mutate(log10_abundance = log10(AbxClass_relativeabundance)) %>% 
  select(-AbxClass_relativeabundance) %>% 
  select(sample_names,log10_abundance,abx_class) %>%
  pivot_wider(names_from = sample_names, values_from = log10_abundance) %>% 
  column_to_rownames(var = "abx_class") %>%
  as.matrix()


# card_heatmap_clr <- identify_card %>% 
#   full_join(reclassify_card) %>% 
#   filter(abx_class %in% c("aminoglycoside","fluoroquinolone","macrolide",
#                           "tetracycline","cephamycin","glycopeptide","sulfonamide","cephalosporin")) %>% 
#   pivot_wider(names_from = sample, values_from = log10_abundance) %>% 
#   column_to_rownames(var = "abx_class") %>%
#   as.matrix()


RA_abx_class <- identify_card %>% 
  full_join(reclassify_card) %>%
  group_by(sample,abx_class) %>% 
  filter(abx_class %in% c("aminoglycoside","fluoroquinolone","macrolide","multi-drug",
                          "tetracycline","cephamycin","glycopeptide","sulfonamide","cephalosporin")) %>% 
  summarise(AbxClass_relativeabundance = sum(relativeabundance)) %>% 
  mutate(dairy = ifelse(sample %in% c("C02_S104","C02_S106",
                                      "C02_S111","C02_S112",
                                      "C02_S113","C02_S115"), 0,1)) %>% 
  group_by(dairy,abx_class) %>% 
  summarise(mean = mean(AbxClass_relativeabundance), sd = sd(AbxClass_relativeabundance), 
            n = n()) %>% 
  arrange(abx_class)


medians <- identify_card %>% 
  full_join(reclassify_card) %>%
  group_by(sample,abx_class) %>% 
  filter(abx_class %in% c("aminoglycoside","fluoroquinolone","macrolide","multi-drug",
                          "tetracycline","cephamycin","glycopeptide","sulfonamide","cephalosporin")) %>% 
  summarise(AbxClass_relativeabundance = sum(relativeabundance)) %>% 
  mutate(dairy = ifelse(sample %in% c("C02_S104","C02_S106",
                                      "C02_S111","C02_S112",
                                      "C02_S113","C02_S115"), 0,1)) %>% 
  group_by(dairy,abx_class) %>% 
  summarise(median = median(AbxClass_relativeabundance), iqr = IQR(AbxClass_relativeabundance), 
            n = n()) %>% 
  arrange(abx_class)

look_RA <- identify_card %>% 
  full_join(reclassify_card) %>%
  group_by(sample,abx_class) %>% 
  #filter(abx_class %in% c("aminoglycoside","fluoroquinolone","macrolide","multi-drug",
#                          "tetracycline","cephamycin","glycopeptide","sulfonamide","cephalosporin"))
filter(abx_class == "cephamycin") %>% 
  summarise(AbxClass_relativeabundance = sum(relativeabundance)) %>% 
  mutate(dairy = ifelse(sample %in% c("C02_S104","C02_S106",
                                      "C02_S111","C02_S112",
                                      "C02_S113","C02_S115"), "community","dairy"))

ggplot(look_RA, aes(y = AbxClass_relativeabundance, x = as.factor(dairy))) +
  geom_point(aes(color = as.factor(dairy), fill = as.factor(dairy), size = 5)) +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  scale_fill_manual(values = c("#00AFBB", "#E7B800")) + 
  ylab("Cephamycin Relative Abundances")
  



clr_amr_tests <- identify_card %>% 
  full_join(reclassify_card) %>%
  group_by(sample,abx_class) %>% 
  filter(abx_class %in% c("aminoglycoside","fluoroquinolone","macrolide",
                          "tetracycline","cephamycin","glycopeptide","sulfonamide","cephalosporin")) %>% 
  group_by(gene_callers_id, group) %>%
  summarise(mean = mean(clr), sd = sd(clr), n = n()) %>%
  filter(n > 2) %>% 
  pivot_wider(names_from=group, values_from=mean:n) %>%
  filter(n_dairy > 2, n_community > 2) %>%
  mutate(t = (mean_community - mean_dairy)/sqrt(sd_community^2/n_community + sd_dairy^2/n_dairy)) %>%
  arrange(desc(abs(t)))
  


library(pheatmap)
pheatmap(card_heatmap, treeheight_row = 0, treeheight_col = 0, cluster_cols = F, cluster_rows = F, 
         cellwidth = 15, cellheight = 12, fontsize = 8, filename = "card_heatmap_log10_simplified.pdf")

pheatmap(card_heatmap_all, treeheight_row = 0, treeheight_col = 0, cluster_cols = F, cluster_rows = F, 
         cellwidth = 15, cellheight = 12, fontsize = 8, filename = "card_heatmap_log10_all.pdf")

clr <- function(x) {
  (x %>% log) - (x %>% log %>% mean)
}


clr_amr_tests <- identify_card %>% 
  mutate("clr" = clr(relativeabundance)) %>%
  group_by(gene_callers_id, group) %>%
  summarise(mean = mean(clr), sd = sd(clr), n = n()) %>%
  filter(n > 2) %>% 
  pivot_wider(names_from=group, values_from=mean:n) %>%
  filter(n_dairy > 2, n_community > 2) %>%
  mutate(t = (mean_community - mean_dairy)/sqrt(sd_community^2/n_community + sd_dairy^2/n_dairy)) %>%
  arrange(desc(abs(t))) ## hmm infinity.... are there too many genes? 
# ^^ This is with each gene at the most granular level without any correction to FDR 

# okay so I have all my relative abundances by AMR group. What if I summed the AMR group to get an idea of proportion 
identify_card %>% 
  filter(card_gene == 1) %>% 
  group_by(sample) %>% 
  summarise(sum_RA = sum(relativeabundance)) 
# sample      sum_RA
# <chr>        <dbl>
#   1 C02_S104 0.000186 
# 2 C02_S106 0.0000633
# 3 C02_S111 0.000583 
# 4 C02_S112 0.000298 
# 5 C02_S113 0.000143 
# 6 C02_S115 0.000202 
# 7 D03_S25  0.000102 
# 8 D03_S26  0.000248 
# 9 D03_S27  0.000198 
# 10 D03_S28  0.000124 
# 11 D03_S31  0.0000833
# 12 D03_S32  0.000351 
# 13 D03_S33  0.0000525
# 14 D03_S34  0.0000787
# 15 D03_S35  0.000224 
# 16 D03_S37  0.000359 
identify_card %>% 
  filter(card_gene == 1) %>% 
  group_by(sample,group) %>% 
  summarise(sum_RA = sum(relativeabundance)) %>% 
  group_by(group) %>% 
  summarise(mean = mean(sum_RA))
# group         mean
# <chr>        <dbl>
# 1 community 0.000246
# 2 dairy     0.000182
# Okay so from this proportion standpoint there's more abundance of AMR genes in the community than in the dairy 
# 0.02% vs. 0.01% (but these are so tiny tiny)



# how much of the gene hits do AMR genes take up 
card_RA %>% 
  group_by(sample) %>% 
  summarise(sum_RA = sum(relativeabundance)) %>% 
  mutate("clr" = clr(sum_RA)) %>% 
  mutate(dairy = ifelse(sample %in% c("C02_S104","C02_S106","C02_S111","C02_S112","C02_S113","C02_S115"),0,1), 
         group = ifelse(sample %in% c("C02_S104","C02_S106","C02_S111","C02_S112","C02_S113","C02_S115"),"community","dairy")) %>% 
  group_by(group)  %>% 
  summarise(mean = mean(clr), sd = sd(clr),n = n()) %>% 
  
  ungroup()


card_RA %>%
  group_by(sample) %>% 
  mutate("clr" = clr(relativeabundance)) %>%
  summarise(mean = mean(clr), sd = sd(clr), n = n()) %>% 
  mutate(dairy = ifelse(sample %in% c("C02_S104","C02_S106","C02_S111","C02_S112","C02_S113","C02_S115"),0,1)) %>% 
  as_tibble()-> mean_clr
test_RA <- t.test(mean ~ dairy, data = mean_clr, paired = FALSE)
test_RA <- glm(mean ~ dairy, data = mean_clr)
summary(test_RA)
# There doesn't appear to be a significant difference in the relative abundance of AMR genes identified 
# using CARD, between dairy vs. community controls 
# Hey that's good! 

# let's read in the motus data 
# let's read in the data 

library(tidyverse)
meta <- read_tsv("/Users/paulinetrinh/Documents/GitHub/statdivlab/HDW/mgx/hdw_mgx_metadata.txt")
motus <- read_tsv("/Users/paulinetrinh/Documents/HDW/metagenomic/HDW_metaphlan/mgx/motus/HDW-motus-mod.txt")
# let's remove species that don't show up at all 
motus[apply(motus[,-1], 1, function(x) !all(x==0)),] -> motus2
dim(motus2)

colSums(motus[,-1])
# the sum isn't perfect but it's pretty close to 1 

View(motus2)

motus_long <- motus2 %>%
  pivot_longer(cols=`HDW-C02-S104`:`HDW-D03-S37`,
               names_to="sample", values_to="abundance") %>%
  filter(abundance > 0) %>%
  group_by(sample) %>%
  right_join(meta, by=c("sample" = "id")) %>%
  select(consensus_taxonomy:group) 
# what does the last line -1 mean? unidentified reads? 

motus_long %>% 
  group_by(sample) %>% 
  summarise(sum_RA = sum(abundance))  # okay so these all sum to 1 by sample

clr <- function(x) {
  (x %>% log) - (x %>% log %>% mean)
}


clr_motus_ttest <- motus_long %>%
  mutate("clr" = clr(abundance)) %>%
  group_by(consensus_taxonomy, group) %>%
  summarise(mean = mean(clr), sd = sd(clr), n = n()) %>%
  filter(n > 2) %>% 
  pivot_wider(names_from=group, values_from=mean:n) %>%
  filter(n_dairy > 2, n_community > 2) %>%
  mutate(t = (mean_community - mean_dairy)/sqrt(sd_community^2/n_community + sd_dairy^2/n_dairy)) %>%
  arrange(desc(abs(t)))



