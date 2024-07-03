library(tidyverse)

hmp_otus_rademu <- readRDS("hmp_dairy_rademu_otus.RDS")
hmp_meta_rademu <- readRDS("hmp_dairy_radEmu_meta.RDS") 

otus_community <- hmp_otus_rademu[which(hmp_meta_rademu$group == "community"), ]
otus_dairy <- hmp_otus_rademu[which(hmp_meta_rademu$group == "dairy"), ]
otus_hmp <- hmp_otus_rademu[which(hmp_meta_rademu$group == "HMP"), ]

sum(colSums(otus_community) == 0)
sum(colSums(otus_dairy) == 0)
sum(colSums(otus_hmp) == 0)
