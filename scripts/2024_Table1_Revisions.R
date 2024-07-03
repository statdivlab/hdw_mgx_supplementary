######################################
# April 15, 2024 Addressing Revisions 
######################################
library(curatedMetagenomicData)
library(tidyverse)
# Create table 1 of differences in demographics 
hdw_metadata <- read.delim("data/hdw_mgx_metadata.txt", sep="\t") %>%
  mutate(sample_id = id_number)

############################################
# read in diet data
############################################
diet_data <- read.csv("data/diet_data_20240424.csv") %>% 
  filter(redcap_event_name %in% c("baseline_arm_1","baseline_arm_2")) %>% 
  filter(identifier %in% c("HDW-C02-S104","HDW-C02-S106", "HDW-C02-S111",
                           "HDW-C02-S112","HDW-C02-S113", "HDW-C02-S115", 
                           "HDW-D03-S25",  "HDW-D03-S26",  "HDW-D03-S27",  
                           "HDW-D03-S28",  "HDW-D03-S31",  "HDW-D03-S32",  
                           "HDW-D03-S33",  "HDW-D03-S34", "HDW-D03-S35", 
                           "HDW-D03-S37")) %>% 
  mutate(id_number = identifier)

############################################
############################################
hdw_metadata_clean <- hdw_metadata %>% 
  separate(spiro_height, c('feet', 'inches'), "'", convert = TRUE) %>% 
  mutate(height_cm = (12*feet + inches)*2.54, 
         height_in = (12*feet + inches), 
         bmi = spiro_weight / (height_in)^2 *703, 
         BMI = round(bmi, digits = 0), 
         sex = "male", 
         dairy = ifelse(sample_id %in% c("HDW-C02-S104","HDW-C02-S106", "HDW-C02-S111",
                                         "HDW-C02-S112","HDW-C02-S113", "HDW-C02-S115"),"community", 
                        ifelse(sample_id %in% c("HDW-D03-S25",  "HDW-D03-S26",  "HDW-D03-S27",  
                                                "HDW-D03-S28",  "HDW-D03-S31",  "HDW-D03-S32",  
                                                "HDW-D03-S33",  "HDW-D03-S34", "HDW-D03-S35", "HDW-D03-S37"),"dairy","HMP")), 
         smoke_ever_cleaned = ifelse(dairy == "dairy", smoke_ever, smoke_ever_v3), 
         smoke_ever_table = ifelse(smoke_ever_cleaned == 1, "Yes", 
                                   ifelse(smoke_ever_cleaned == 2, "No", "Declined to answer")), 
         smoke_now_cleaned = ifelse(dairy == "dairy", smoke_now, smoke_now_v3), 
         smoke_now_table = ifelse(smoke_now_cleaned == 1, "Yes", 
                                   ifelse(smoke_now_cleaned == 2, "No", "Declined to answer")), 
         new_id = sample_id) %>% 
  separate(new_id, into = c("hdw","hdw_id"), sep = "HDW-") %>% 
  select(-hdw) %>% 
  full_join(diet_data, by = "id_number") %>% 
  mutate(alcohol_usage = ifelse(dairy == "dairy", food_alcohol, food_alcohol_v3), 
         alcohol = ifelse(alcohol_usage == 0, "Never", 
                          ifelse(alcohol_usage == 1, "Weekly 1-7 times/week", 
                                 ifelse(alcohol_usage == 2, "Weekly 1-7 times/week", 
                                        ifelse(alcohol_usage == 3, "Rarely <1/week", 
                                               "Rarely <1/week ")))), 
         dairy_usage = ifelse(dairy == "dairy", food_dairy, food_dairy_v3), 
         dairy_food = ifelse(dairy_usage == 0, "Never", 
                        ifelse(dairy_usage == 1, "Weekly 1-7 times/week", 
                               ifelse(dairy_usage == 2, "Weekly 1-7 times/week", 
                                      ifelse(dairy_usage == 3, "Rarely <1/week", 
                                             "Rarely <1/week")))), 
         veg_usage = ifelse(dairy == "dairy", food_veg, food_veg_v3) , 
         veg = ifelse(veg_usage == 0, "Never", 
                        ifelse(veg_usage == 1, "Weekly 1-7 times/week", 
                               ifelse(veg_usage == 2, "Weekly 1-7 times/week", 
                                      ifelse(veg_usage == 3, "Rarely <1/week", 
                                             "Rarely <1/week")))),
         egg_usage = ifelse(dairy == "dairy", food_eggs, food_eggs_v3), 
         egg = ifelse(egg_usage == 0, "Never", 
                      ifelse(egg_usage == 1, "Weekly 1-7 times/week", 
                             ifelse(egg_usage == 2, "Weekly 1-7 times/week", 
                                    ifelse(egg_usage == 3, "Rarely <1/week", 
                                           "Rarely <1/week")))),
         beef_usage = ifelse(dairy == "dairy", food_beef, food_beef_v3), 
         beef = ifelse(beef_usage == 0, "Never", 
                      ifelse(beef_usage == 1, "Weekly 1-7 times/week", 
                             ifelse(beef_usage == 2, "Weekly 1-7 times/week", 
                                    ifelse(beef_usage == 3, "Rarely <1/week", 
                                           "Rarely <1/week")))),
         pork_usage = ifelse(dairy == "dairy", food_pork, food_pork_v3), 
         pork = ifelse(pork_usage == 0, "Never", 
                       ifelse(pork_usage == 1, "Weekly 1-7 times/week", 
                              ifelse(pork_usage == 2, "Weekly 1-7 times/week", 
                                     ifelse(pork_usage == 3, "Rarely <1/week", 
                                            "Rarely <1/week")))),
         chicken_usage = ifelse(dairy == "dairy", food_chicken, food_chicken_v3), 
         chicken = ifelse(chicken_usage == 0, "Never", 
                       ifelse(chicken_usage == 1, "Weekly 1-7 times/week", 
                              ifelse(chicken_usage == 2, "Weekly 1-7 times/week", 
                                     ifelse(chicken_usage == 3, "Rarely <1/week", 
                                            "Rarely <1/week")))),
         fish_usage = ifelse(dairy == "dairy", food_fish, food_fish_v3), 
         fish = ifelse(fish_usage == 0, "Never", 
                          ifelse(fish_usage == 1, "Weekly 1-7 times/week", 
                                 ifelse(fish_usage == 2, "Weekly 1-7 times/week", 
                                        ifelse(fish_usage == 3, "Rarely <1/week", 
                                               "Rarely <1/week")))),
         lamb_usage = ifelse(dairy == "dairy", food_lamb, food_lamb_v3), 
         lamb = ifelse(lamb_usage == 0, "Never", 
                       ifelse(lamb_usage == 1, "Weekly 1-7 times/week", 
                              ifelse(lamb_usage == 2, "Weekly 1-7 times/week", 
                                     ifelse(lamb_usage == 3, "Rarely <1/week", 
                                            "Rarely <1/week"))))) 

saveRDS(hdw_metadata_clean, "data/hdw_metadata_rademu.RDS")
my_data <- curatedMetagenomicData::sampleMetadata
hmp_done <- my_data %>% 
  filter(study_name == "HMP_2012")  %>% 
  filter(body_site == "stool")

hmp47_ids <- read.table("data/hmp47_ids.txt")
colnames(hmp47_ids) <- "sample_id"

hmp47_meta <- hmp_done %>% right_join(hmp47_ids)
saveRDS(hmp47_meta, "hmp47_meta.RDS")

manuscript_table1 <- hmp47_meta %>% 
  mutate(sex = gender) %>% 
  select(-gender) %>% 
  full_join(hdw_metadata_clean) %>% 
  mutate(group = ifelse(sample_id %in% c("HDW-C02-S104","HDW-C02-S106", "HDW-C02-S111",
                                         "HDW-C02-S112","HDW-C02-S113", "HDW-C02-S115"),"community", 
                        ifelse(sample_id %in% c("HDW-D03-S25",  "HDW-D03-S26",  "HDW-D03-S27",  
                                                "HDW-D03-S28",  "HDW-D03-S31",  "HDW-D03-S32",  
                                                "HDW-D03-S33",  "HDW-D03-S34", "HDW-D03-S35", "HDW-D03-S37"),"dairy","HMP")), 
         occupation = ifelse(group == "dairy", "dairy worker", 
                             ifelse(group == "community", "field worker", "unknown")), 
         study = ifelse(group %in% c("dairy","community"), "HDW","HMP"), 
         bmi_cat = ifelse(BMI <= 18.5 , "<18.5", 
                          ifelse(BMI > 18.5 & BMI <=25, "18.5-24.9", 
                                 ifelse(BMI > 25 & BMI <= 29.9, "25-29.9",
                                        ifelse(29.9 < BMI, "30.0+", NA))))) %>% 
  mutate(sample = ifelse(sample_id == "HDW-C02-S104", "C02_S104",
                         ifelse(sample_id == "HDW-C02-S106","C02_S106",
                                ifelse(sample_id == "HDW-C02-S111","C02_S111",
                                       ifelse(sample_id == "HDW-C02-S112", "C02_S112", 
                                              ifelse(sample_id == "HDW-C02-S113","C02_S113",
                                                     ifelse(sample_id == "HDW-C02-S115","C02_S115",
                                                            ifelse(sample_id == "HDW-D03-S25","D03_S25",
                                                                   ifelse(sample_id == "HDW-D03-S26","D03_S26",
                                                                          ifelse(sample_id == "HDW-D03-S27","D03_S27",
                                                                                 ifelse(sample_id == "HDW-D03-S28","D03_S28",
                                                                                        ifelse(sample_id == "HDW-D03-S31","D03_S31",
                                                                                               ifelse(sample_id == "HDW-D03-S32","D03_S32",
                                                                                                      ifelse(sample_id == "HDW-D03-S33","D03_S33",
                                                                                                             ifelse(sample_id == "HDW-D03-S34","D03_S34",
                                                                                                                    ifelse(sample_id == "HDW-D03-S35","D03_S35",
                                                                                                                           ifelse(sample_id == "HDW-D03-S37","D03_S37",sample_id))))))))))))))))) %>% 
  full_join(hdw_read_info) %>% 
  mutate(total_reads = ifelse(group %in% c("community","dairy"), number.of.pairs.analyzed, number_reads))
## let's look at age 
## let's look at grown up on a farm 



library(boot)
library(htmlTable)
library(table1)

manuscript_table1$veg <- factor(manuscript_table1$veg, levels = c("Rarely <1/week","Weekly 1-7 times/week"))
manuscript_table1$egg <- factor(manuscript_table1$egg, levels = c("Rarely <1/week","Weekly 1-7 times/week"))
manuscript_table1$beef <- factor(manuscript_table1$beef, levels = c("Rarely <1/week","Weekly 1-7 times/week"))
manuscript_table1$chicken <- factor(manuscript_table1$chicken, levels = c("Rarely <1/week","Weekly 1-7 times/week"))
manuscript_table1$lamb <- factor(manuscript_table1$lamb, levels = c("Never","Rarely <1/week"))
manuscript_table1$fish <- factor(manuscript_table1$fish, levels = c("Never","Rarely <1/week","Weekly 1-7 times/week"))
manuscript_table1$bmi_cat <- factor(manuscript_table1$bmi_cat, levels = c("18.5-24.9","25-29.9","30.0+"))


label(manuscript_table1$sex) <- "Sex"
label(manuscript_table1$age) <- "Age"
label(manuscript_table1$occupation) <- "Occupation"
label(manuscript_table1$BMI) <- "Body Mass Index"
label(manuscript_table1$smoke_now_table) <- "Current Smoker"
label(manuscript_table1$alcohol) <- "Alcohol Consumption"
label(manuscript_table1$dairy_food) <- "Dairy Consumption"
label(manuscript_table1$veg) <- "Vegetable Consumption"
label(manuscript_table1$egg) <- "Eggs Consumption"
label(manuscript_table1$beef) <- "Beef Consumption"
label(manuscript_table1$chicken) <- "Chicken Consumption"
label(manuscript_table1$lamb) <- "Lamb Consumption"
label(manuscript_table1$fish) <- "Fish Consumption"
label(manuscript_table1$work_time_thisfarm) <- "Years of Dairy Work"
label(manuscript_table1$bmi_cat) <- "BMI"

units(manuscript_table1$age) <- "years"
table1_obj <- table1(~ sex + age  + bmi_cat + total_reads | study, data = manuscript_table1,overall=F, extra.col=list(`P-value`=pvalue))
tab1_df <- as.data.frame(table1_obj)
write.csv(tab1_df,"table1_csv.csv", row.names = F)
saveRDS(manuscript_table1, "manuscript_table1.RDS")
manuscript_table1 <- readRDS("data/manuscript_table1.RDS")


pvalue <- function(x, ...) {
  # Construct vectors of data y, and groups (strata) g
  y <- unlist(x)
  g <- factor(rep(1:length(x), times=sapply(x, length)))
  if (is.numeric(y)) {
    # For numeric variables, perform a standard 2-sample t-test
    p <- t.test(y ~ g)$p.value
  } else {
    # For categorical variables, perform a chi-squared test of independence
    p <- chisq.test(table(y, g))$p.value
  }
  # Format the p-value, using an HTML entity for the less-than sign.
  # The initial empty string places the output on the line below the variable label.
  c("", sub("<", "&lt;", format.pval(p, digits=3, eps=0.001)))
}

hdw_table1 <- manuscript_table1 %>% 
  filter(group %in% c("community", "dairy"))

hdw_table1$veg <- factor(hdw_table1$veg, levels = c("Rarely <1/week","Weekly 1-7 times/week"))
hdw_table1$egg <- factor(hdw_table1$egg, levels = c("Rarely <1/week","Weekly 1-7 times/week"))
hdw_table1$beef <- factor(hdw_table1$beef, levels = c("Rarely <1/week","Weekly 1-7 times/week"))
hdw_table1$chicken <- factor(hdw_table1$chicken, levels = c("Rarely <1/week","Weekly 1-7 times/week"))
hdw_table1$lamb <- factor(hdw_table1$lamb, levels = c("Never","Rarely <1/week"))
hdw_table1$fish <- factor(hdw_table1$fish, levels = c("Never","Rarely <1/week","Weekly 1-7 times/week"))
hdw_table1$bmi_cat <- factor(hdw_table1$bmi_cat, levels = c("18.5-24.9","25-29.9","30.0+"))


label(hdw_table1$sex) <- "Sex"
label(hdw_table1$age) <- "Age"
label(hdw_table1$occupation) <- "Occupation"
label(hdw_table1$bmi_cat) <- "Body Mass Index"
label(hdw_table1$smoke_now_table) <- "Current Smoker"
label(hdw_table1$alcohol) <- "Alcohol Consumption"
label(hdw_table1$dairy_food) <- "Dairy Consumption"
label(hdw_table1$veg) <- "Vegetable Consumption"
label(hdw_table1$egg) <- "Eggs Consumption"
label(hdw_table1$beef) <- "Beef Consumption"
label(hdw_table1$chicken) <- "Chicken Consumption"
label(hdw_table1$lamb) <- "Lamb Consumption"
label(hdw_table1$fish) <- "Fish Consumption"

tbl1_hdw <- table1(~ age + bmi_cat  + smoke_now_table + alcohol + dairy_food + veg + egg + beef + chicken + 
         lamb + fish + occupation + total_reads| group, data = hdw_table1, overall=F, extra.col=list(`P-value`=pvalue))

tbl1_hdw_df <- as.data.frame(tbl1_hdw)
write.csv(tbl1_hdw_df, "data/tbl1_hdw_df.csv", row.names = F)
######################################################
# Differences in sequencing depths 
######################################################

hdw_read_info <- read.csv("data/HDW-QC-report.csv") %>% 
  select(sample, number.of.pairs.analyzed)

t.test(total_reads ~ group, data = hdw_table1)
t.test(age ~ group, data = hdw_table1)

chisq.test(table(hdw_table1$group,hdw_table1$bmi_cat))
chisq.test(table(hdw_table1$group,hdw_table1$alcohol))
chisq.test(table(hdw_table1$group,hdw_table1$smoke_now_table))
chisq.test(table(hdw_table1$group,hdw_table1$dairy_food))
chisq.test(table(hdw_table1$group,hdw_table1$veg))
chisq.test(table(hdw_table1$group,hdw_table1$egg))
chisq.test(table(hdw_table1$group,hdw_table1$beef))
chisq.test(table(hdw_table1$group,hdw_table1$chicken))
chisq.test(table(hdw_table1$group,hdw_table1$lamb))
chisq.test(table(hdw_table1$group,hdw_table1$fish))


min(hdw_read_info$number.of.pairs.analyzed)
max(hdw_read_info$number.of.pairs.analyzed)

mean(hdw_read_info$number.of.pairs.analyzed)
IQR(hdw_read_info$number.of.pairs.analyzed)


t.test(total_reads ~ study, data = manuscript_table1)
t.test(age ~ study, data = manuscript_table1)
chisq.test(table(manuscript_table1$group,manuscript_table1$bmi_cat))


#############################
## look at sequencing stats #
#############################
seq_stats <- read_csv("data/hdw_seq_stats.csv")
seq_stats$group <- factor(seq_stats$group, levels = c("community","dairy"))

table2 <- table1(~`number of pairs analyzed` + `total pairs filtered` + `total pairs passed` +  `host read pairs filtered` + `total pairs after host removal` + `number of contigs` + `num genes (prodigal)`| group, data = seq_stats)

table2_df <- as.data.frame(table2)
write.csv(table2_df, "table2_df.csv", row.names = F)


t.test(`number of pairs analyzed` ~ group, data = seq_stats) # 0.03547
t.test(`total pairs filtered`~ group, data = seq_stats) # p = 0.2011
t.test(`total pairs passed` ~ group, data = seq_stats) # p = 0.033
t.test(`host read pairs filtered` ~ group, data = seq_stats) # p = 0.49
t.test(`total pairs after host removal` ~ group, data = seq_stats) # p = 0.03264
t.test(`number of contigs` ~ group, data = seq_stats) # p = 0.0237
t.test(`num genes (prodigal)` ~ group, data = seq_stats) # p = 0.0164

#############################
# read in seq stats for HMP #
#############################

csv.files.hmp <- list.files(path="data/HMP_contigs_summary", full.names = TRUE) %>% 
  lapply(read_tsv)

reformatted_hmpdata <- list()
for (i in 1:length(csv.files.hmp)){
  reformatted_hmpdata[[i]] <- csv.files.hmp[[i]] %>% 
    as.data.frame() %>% 
    t() %>% 
    janitor::row_to_names(row_number = 1) %>% 
    as.data.frame() %>%
    rownames_to_column("sample") %>% 
    separate(sample, c("sample_id", "leftover"), sep = ".contigs") %>% 
    select(sample_id, `Num Contigs`, `Num Genes (prodigal)`) %>% 
    rename(`number of contigs` = `Num Contigs`, 
           `num genes (prodigal)` = `Num Genes (prodigal)`) %>% 
    mutate(group = "HMP", 
           `number of contigs` = as.numeric(`number of contigs`), 
           `num genes (prodigal)` = as.numeric(`num genes (prodigal)`))
}

sequencing_info <- hmp47_meta %>% select(sample_id, number_reads) %>% 
  rename(`total pairs after host removal` = number_reads)

hmp_contigs_data <- do.call("rbind",reformatted_hmpdata) %>% 
  as_tibble() %>% 
  right_join(sequencing_info) %>% 
  full_join(seq_stats) %>% 
  mutate(HMP = ifelse(group == "HMP", "HMP", "HDW"))
saveRDS(hmp_contigs_data, "data/hmp_contigs_data.RDS")

table2_HMP <- table1(~`number of pairs analyzed` + `total pairs filtered` + `total pairs passed` +  `host read pairs filtered` + `total pairs after host removal` + `number of contigs` + `num genes (prodigal)`| group, data = hmp_contigs_data)
table2_df <- as.data.frame(table2_HMP)

t.test(`number of contigs` ~ HMP, data = hmp_contigs_data) # p = 3.633e-09
t.test(`num genes (prodigal)` ~ HMP, data = hmp_contigs_data) # p = 0.368
t.test(`total pairs after host removal` ~ HMP, data = hmp_contigs_data) # p < 2.2e-16 


