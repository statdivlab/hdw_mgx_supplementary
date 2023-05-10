# Format manifest file to pull appropriate samples from HMP 
library(tidyverse)
HMP <- read_tsv("/Users/paulinetrinh/Documents/HDW/metagenomic/hmp_manifest_metadata_31e9c3cd44.tsv") 
# 14765     9
HMP2 <- read_tsv("/Users/paulinetrinh/Downloads/hmp_manifest_metadata_32194fc63d.tsv")
HMP_unique <- distinct(HMP) %>% filter(sample_body_site == "feces") %>% filter(visit_number == 1)
# 248 x 9 

HMP_links <- read_tsv("/Users/paulinetrinh/Documents/HDW/metagenomic/hmp_manifest_1a1d0ad5cf.tsv")
#14765     5
HMP_manifest <- HMP_unique %>% left_join(HMP_links, by = "sample_id") %>% 
  filter(str_detect(urls, 'hmasm2')) %>% 
  select(file_id, md5, size, urls, sample_id)

dim(HMP_manifest)
#239x5

write.table(HMP_manifest, "HMP_manifest.tsv", row.names = FALSE, quote=FALSE, sep='\t')
