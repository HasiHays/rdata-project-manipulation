## Script - Manupulate gene expression data from NCBI GEO datasets (GSE183947)
# ssetwd("~/Desktop/R/data-manipulation/Script")
# R version 4.1.2

# libraries
library(GEOquery)
library(tidyverse)
library(dplyr)

# Load the data
my_data <- read.csv(file = "../Data/GSE183947_fpkm.csv") 
dim(my_data)

# Get meta (clinical file that associate with gene expression data) data via GEOquery
gse_data <- getGEO(GEO = 'GSE183947', GSEMatrix = TRUE)
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 1000)

gse_data

metadata <- pData(phenoData(gse_data[[1]]))
head(metadata)

metadata_subset <- select(metadata, c(1,10,11,17))

#renaming a column
metadat_modified <- metadata %>% ###use pipe operator rename and change var
  select(1,10,11,17) %>%
  rename(tissue = characteristics_ch1) %>% #renamin a columnn
  rename(metastasis = characteristics_ch1.1) %>% #renamin a columnn
  mutate(tissue = gsub("tissue: ","", tissue)) %>%  # rename the column var (remove tissue:)
  mutate(metastasis = gsub("metastasis: ","", metastasis)) 
 

head(my_data)

# reshaping data
data_long <- my_data %>%
  rename(gene= X) %>%
  gather(key = 'samples', value = 'FPKM', -gene)

# join two data sets
 both_data <- data_long %>%
  left_join(.,metadat_modified, by = c("samples" = "description")) 

 #explore data with top 10 largest values (example)
 
 top10 <- both_data %>%
   arrange(desc(gene)) %>%
   slice(1:10)
 head()
 
 #explore data (compare two tumer related genes with tumer vs normal sample)

 arranged_data <- both_data %>%
   filter(gene == 'BRCA1' | gene== 'BRCA2') %>% ## select only two genes (| = or)
   group_by(gene,tissue) %>% #select tow column to cal mean
   summarize(mean_FPKM = mean(FPKM),
             median_FPKM = median(FPKM)) %>%  ## cal mean
   arrange(mean_FPKM) #arrange ascending order

#study dplyr pakage more for more handson experience
##########END##########  
 
