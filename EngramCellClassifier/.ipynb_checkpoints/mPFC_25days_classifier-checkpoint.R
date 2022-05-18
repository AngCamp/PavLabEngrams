#
#
#
# first attempt at a classifier
setwd("C:/Users/angus/Desktop/PavLabEngrams/EngramCellClassifier")

library(tidyverse)
library(GEOquery)
library(AnnotationDbi)
library(randomForest)
library(data.table)
library(reshape2)
library(FactoMineR)
library(factoextra)
library(Rtsne)
library(dplyr)
library(Seurat)
library(stringr)
library(patchwork)
library(metap)




# Chen 2020 ---------------------------------------------------------------


chen2020_counts <- read.csv('Chen2020_GSE152632/GSE152632_GEO_mberkchen_TRAP2_counts.csv.gz', header = TRUE)
rownames(chen2020_counts) <- chen2020_counts$X
chen2020_counts <- chen2020_counts[,2:3531]
chen2020_meta <- read.csv( 'Chen2020_GSE152632/SraRunTable.txt', header = TRUE)

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
chen2020_counts <- chen2020_counts[!c(1:23355) %in% underexpressedgenes,]
rownames(chen2020_meta) <- colnames(chen2020_counts)


#Make seurate object
chen2020 <- CreateSeuratObject(counts = chen2020_counts)

chen2020 <- FindVariableFeatures(chen2020, selection.method = "vst", nfeatures = 3530) #no idea how I chose 3530, 
# in the tutorial they chose 2000, I just chose to match the number of cells we had
chen2020nobatch <- ScaleData(chen2020)

chen2020 <- AddMetaData(chen2020, chen2020_meta)

#create index to split the chen data with
tdTpos.idx <- which(chen2020_meta$engram_label=="tdT+")
tdTneg.idx <- which(chen2020_meta$engram_label=="tdT-")

