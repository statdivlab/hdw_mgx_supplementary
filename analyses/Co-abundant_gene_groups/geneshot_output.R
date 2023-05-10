# GeneShot Results 
# I exported the pandas dataframes into csv's using Python 
# See geneshot_Singularity.txt for the code 
#install.packages("devtools")
#devtools::install_github("rstudio/reticulate", dependencies = TRUE)
library("reticulate")

# Use the path to Python on your system
use_python("/Users/paulinetrinh/Library/r-miniconda/envs/r-reticulate/bin/python")

# Create and use a new virtual enviornment named "r-reticulate"
virtualenv_create("r-reticulate")
use_virtualenv("r-reticulate") 

# Install packages required to read hdf5 data
py_install("pandas")
py_install("scipy")
py_install("pytables")
py_install("h5py")
py_install("rpy2")

# Import those packages into the environment
h5py <- import("h5py")
rpy2 <- import("rpy2")
#rpy2$robjects
#rpy2_ro <- import("rpy2.robjects")
#rpy2_pandas2ri <- import("rpy2.robjects.pandas2ri")
#rpy2_pandas2ri$py2rpy
pd <- import("pandas", convert = FALSE)
np <- import("numpy", convert = FALSE)

# Using the pandas read_hdf function (pd.read_hdf('file_path', key='your_group')

#Geneshot outputs files ("keys")
#'/manifest'                          '/ref/taxonomy'
#'/summary/all'                       '/ordination/pca'
#'/summary/breakaway'                 '/ordination/tsne'
#'/summary/experiment'                '/annot/gene/all'
#'/summary/genes_aligned'             '/annot/gene/cag'
#'/summary/genes_assembled'           '/annot/gene/eggnog'
#'/summary/readcount'                 '/annot/gene/tax'
#'/stats/cag/corncob'                 '/annot/cag/all'
#'/stats/cag/corncob_wide'
#'/abund/gene/wide'
#"/abund/cag/wide'

key <- pd$read_hdf("/Users/paulinetrinh/Documents/HDW/geneshot/output2/geneshot.results.hdf5", key = "/annot/gene/all", mode ='r+')
key_df <- py_to_r(key)
key_df

# Create an R data.frame for each of the geneshot outputs using the py_to_r function from the reticulate package
key <- pd$read_hdf("/Users/paulinetrinh/Documents/HDW/geneshot/output2/geneshot.results.hdf5", key = "key", mode ='r+')
key_df <- py_to_r(key)
key_df


corncob <- pd$read_hdf("/Users/paulinetrinh/Documents/HDW/geneshot/output2/geneshot.results.hdf5", key = "/stats/cag/corncob", mode ='r+')
corncob_df <- py_to_r(corncob)
corncob_df

library(tidyverse)
corncob_df %>% 
  filter(parameter == "dairy") %>% 
  arrange(q_value) -> corncob_DA
write_csv(corncob_DA,"/Users/paulinetrinh/Documents/HDW/geneshot/output2/tables/corncob_DA.csv")
corncob_DA <- read_csv("/Users/paulinetrinh/Documents/HDW/geneshot/output2/tables/corncob_DA.csv")
png("cag_volcano_plot.png", units="in", width=6, height=4, res=300)  
ggplot(data=corncob_DA, aes(x=estimate, y=-log10(p_value))) + 
  geom_point(size = 3, alpha = 0.6) + 
  theme_minimal() + 
  xlim(-2,20)
dev.off()

library(breakaway)
library(DivNet)
breakaway <- pd$read_hdf("/Users/paulinetrinh/Documents/HDW/geneshot/output2/geneshot.results.hdf5", key = "/summary/breakaway", mode ='r+')
breakaway_df <- py_to_r(breakaway)
breakaway_df %>% 
  mutate(dairy = ifelse(specimen %in% c("HDW-C02-S104","HDW-C02-S113","HDW-C02-S106",
                                        "HDW-C02-S115","HDW-C02-S112","HDW-C02-S111"),0,1)) -> breakaway_df
# I can test for differences in richness 
breakaway_gene_richness <- betta(chats = breakaway_df$estimate,
                      ses = breakaway_df$error,
                      X = model.matrix(~dairy, data = breakaway_df))
breakaway_gene_richness$table

# $table
#Estimates Standard Errors p-values
#(Intercept)  938242.9        53468.01    0.000
#dairy       -204418.8        67632.30    0.003
## There is a significant difference in protein coding gene richness where there is less gene richness in the dairy 
## worker metagenomes than the community members metagenomes
summary_readcount <- pd$read_hdf("/Users/paulinetrinh/Documents/HDW/geneshot/output2/geneshot.results.hdf5", key = "/summary/readcount", mode ='r+')

summary <- pd$read_hdf("/Users/paulinetrinh/Documents/HDW/geneshot/output2/geneshot.results.hdf5", key = "/summary/all", mode ='r+')
summary_df <- py_to_r(summary)
summary_ordered <- summary_df %>% 
  arrange(n_genes_assembled)

#cag_annot <- pd$read_hdf("/Users/paulinetrinh/Documents/HDW/geneshot/output2/geneshot.results.hdf5", key = "/annot/gene/cag", mode ='r+')
#cag_annot_df <- py_to_r(cag_annot)

#gene_eggnog <- pd$read_hdf("/Users/paulinetrinh/Documents/HDW/geneshot/output2/geneshot.results.hdf5", key = "annot/gene/eggnog", mode ='r+')
#gene_eggnog_df <- py_to_r(gene_eggnog)

#gene_tax <- pd$read_hdf("/Users/paulinetrinh/Documents/HDW/geneshot/output2/geneshot.results.hdf5", key = "/annot/gene/tax", mode ='r+')
#gene_tax_df <- py_to_r(gene_tax)

gene_all <- pd$read_hdf("/Users/paulinetrinh/Documents/HDW/geneshot/output2/geneshot.results.hdf5", key = "/annot/gene/all", mode ='r+')
gene_all_df <- py_to_r(gene_all)
saveRDS(gene_all_df,"gene_all_df.RDS")
gene_all_df %>% 
  filter(CAG == 28) -> CAG_28
## Hmm there are a lot of genes here, about 5941 genes 
dim(CAG_28)

summary <- CAG_28 %>% 
  group_by(tax_name,tax_rank) %>% 
  summarise(n = n()) %>% 
  arrange(-n)

summary_thresholds  <- summary %>% 
  filter(tax_rank == "species") %>% 
  filter(n > 30)

df_unlist<-function(df){
  df<-as.data.frame(df)
  nr<-nrow(df)
  c.names<-colnames(df)
  lscols<-as.vector(which(apply(df,2,is.list)==TRUE))
  if(length(lscols)!=0){
    for(i in lscols){
      temp<-as.vector(unlist(df[,i]))
      if(length(temp)!=nr){
        adj<-nr-length(temp)
        temp<-c(rep(0,adj),temp)
      }
      df[,i]<-temp
    } #end for
    df<-as.data.frame(df)
    colnames(df)<-c.names
  }
  return(df)
}

newDF <- df_unlist(gene_all_df)
for (i in 0:58){
  newDF %>% 
    filter(CAG == i) -> cag_data
  write_csv(cag_data,paste0("/Users/paulinetrinh/Documents/HDW/geneshot/output2/tables/CAG_",i,".csv"))
}
gene_all_df %>% 
  filter(CAG == 2)
gene_all_df <- as_tibble(gene_all_df)
write.table(gene_all_df,"/Users/paulinetrinh/Documents/HDW/geneshot/output2/tables/all_genes.txt")


#cag_all <- pd$read_hdf("/Users/paulinetrinh/Documents/HDW/geneshot/output2/geneshot.results.hdf5", key = "/annot/cag/all", mode ='r+')
#cag_all_df <- py_to_r(cag_all)
