# Hrvatin Clustering 

# inDrop explanation here: https://www.illumina.com/science/sequencing-method-explorer/kits-and-arrays/indrop.html

# 10x seq https://medicine.yale.edu/keck/ycga/sequencing/10x/singcellsequencing/

~/test_datasets/Hrvatin2018


#libraries
library(Seurat)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(sampler)
library(randomForest)
library(rfUtilities)
library(tictoc)
library(caTools)
library(pROC)
library(stats)




#loading the data
# note it is too big to be stored on the github 
# claims: generating a final dataset containing 47,209 cells

Hrvatin2018_counts <- read.csv('~/test_datasets/Hrvatin2018/GSE102827_merged_all_raw.csv.gz', header = TRUE)
rownames(Hrvatin2018_counts) <- Hrvatin2018_counts$X
Hrvatin2018_counts <- Hrvatin2018_counts[,2:3531] #25187 genes, 3530 cells

#meta is for each sample, tells you sex and strain age etc
Hrvatin2018_meta <- read.csv('~/test_datasets/Hrvatin2018/SraRunTable.txt', header = TRUE)

#this contains the most important meta data
#65539 cells, 6 levels of details, mouse strain, major cell type, then subtypes
Hrvatin2018_celltypelabels <- read.csv('~/test_datasets/Hrvatin2018/GSE102827_cell_type_assignments.csv.gz', header = TRUE)


Hrvatin2018_MATRIX <- read.csv('~/test_datasets/Hrvatin2018/GSE102827_MATRIX.csv.gz', header = TRUE)

filtered_cell_labels <- Hrvatin2018_celltypelabels[ Hrvatin2018_celltypelabels$X %in% colnames(Hrvatin2018_counts),]

# data is spread across 40 datasets, but they are listed as bulk 
# https://www.ncbi.nlm.nih.gov/sra
