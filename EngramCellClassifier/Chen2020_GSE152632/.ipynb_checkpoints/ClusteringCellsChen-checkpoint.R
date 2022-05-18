# 
# 
# 
# This is the cell clustering an cleaning for Chen et al., (2020)
# 
# Chen, M. B., Jiang, X., Quake, S. R., & Südhof, T. C. (2020).
# Persistent transcriptional programmes are associated with remote memory.
# Nature, 587(7834), 437-442.

#Path on pavlab server
#setwd("/home/acampbell/PavLabEngrams/EngramCellClassifier/Chen2020_GSE152632")

#libraries
library(Seurat)
library(dplyr)
library(tidyverse)
library(ggplot2)


#Load the data and Metadata
chen2020_counts <- read.csv('GSE152632_GEO_mberkchen_TRAP2_counts.csv.gz', header = TRUE)
rownames(chen2020_counts) <- chen2020_counts$X
chen2020_counts <- chen2020_counts[,2:3531]
chen2020_meta <- read.csv( 'SraRunTable.txt', header = TRUE)

#add engram label
chen2020_meta$engram_label <-  as.factor(sapply(as.character(colnames(chen2020_counts)), function(y) if (grepl("_pos_", y, fixed=TRUE)) "tdT+" else "tdT-"))

#create the condition label
condition_label <- chen2020_meta$source_name %>%
  sapply( function(y) if (grepl("Homecage", y, fixed=TRUE)) "Homecage")

condition_label <- chen2020_meta$source_name
condition_label[str_detect(condition_label, "Homecage")] = "Homecage"
condition_label[str_detect(condition_label, "Context-Only")] = "Context-Only"
condition_label[str_detect(condition_label, "Fear-Only")] = "Fear-Only"
condition_label[str_detect(condition_label, "Fear-Recall")] = "Fear-Recall"
chen2020_meta$condition_label <- condition_label

#adding cell bacrcodes from coutn data to rows of metadata for seurat
chen2020_counts <- chen2020_counts[!c(1:23355) %in% underexpressedgenes,] #this produces an error
rownames(chen2020_meta) <- colnames(chen2020_counts)




#Normal Seurat Workflow
#from here: https://satijalab.org/seurat/articles/pbmc3k_tutorial.html

chen <- CreateSeuratObject(counts = hochgerner5k_2018_counts[, ])


