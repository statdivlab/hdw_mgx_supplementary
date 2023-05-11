# script running CLR tests at the species and phylum level 
library(data.table)

clr <- function(x) {
  (x %>% log) - (x %>% log %>% mean)
}

all_long_counts_phylum <- readRDS("all_long_counts_phylum.RDS")

clr_tests_counts_phylum <- all_long_counts_phylum %>%
  group_by(sample) %>% 
  mutate("clr" = clr(abundance_counts)) %>%
  group_by(clade_name, group) %>%
  summarise(mean = mean(clr), sd = sd(clr), n = n()) %>%
  filter(n > 2) %>% 
  pivot_wider(names_from=group, values_from=mean:n) %>%
  filter(n_dairy > 2, n_community > 2) %>%
  mutate(t = (mean_community - mean_dairy)/sqrt(sd_community^2/n_community + sd_dairy^2/n_dairy),
         df = n_community + n_dairy - 2, 
         pvalue = 2*pt(-abs(t),df=df)) %>%
  arrange(abs(pvalue)) 

clr_tests_counts_phylum$qvalue <- qvalue::qvalue(clr_tests_counts_phylum$pvalue, pi0.method = "bootstrap", lambda = 0.5)$qvalues
write.csv(clr_tests_counts_phylum,"STable2.csv")


all_long_counts <- readRDS("all_long_counts.RDS")

clr_tests_counts <- all_long_counts %>%
  group_by(sample) %>% 
  mutate("clr" = clr(abundance_counts)) %>%
  group_by(clade_name, group) %>%
  summarise(mean = mean(clr), sd = sd(clr), n = n()) %>%
  filter(n > 2) %>% 
  pivot_wider(names_from=group, values_from=mean:n) %>%
  filter(n_dairy > 2, n_community > 2) %>%
  mutate(t = (mean_community - mean_dairy)/sqrt(sd_community^2/n_community + sd_dairy^2/n_dairy),
         df = n_community + n_dairy - 2, 
         pvalue = 2*pt(-abs(t),df=df)) %>%
  arrange(abs(pvalue)) 

clr_tests_counts$qvalue <- qvalue::qvalue(clr_tests_counts$pvalue)$qvalues
write.csv(clr_tests_counts,"STable1.csv")

