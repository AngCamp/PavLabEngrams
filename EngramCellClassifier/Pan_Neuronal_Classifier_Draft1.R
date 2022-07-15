# PAN NEURONAL CLASSIFIER


# test set, subiculum data by cembrowski: https://www.sciencedirect.com/science/article/pii/S0092867418303118

# human dataset from Allen
# https://www.nature.com/articles/s41586-019-1506-7#data-availability


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
library(stringr)
library(Dict)
library(pheatmap)



# you need to load the data here, and filter genes accordingly

#Hrvatin2018_counts
Hrvatin2018_counts <- read.csv('~/test_datasets/Hrvatin2018/GSE102827_MATRIX.csv.gz', header = TRUE)
rownames(Hrvatin2018_counts) <- Hrvatin2018_counts$X
Hrvatin2018_counts <- Hrvatin2018_counts[,c(2:dim(Hrvatin2018_counts)[2])]

#meta is for each sample, tells you sex and strain age etc
Hrvatin2018_meta <- read.csv('~/test_datasets/Hrvatin2018/SraRunTable.txt', header = TRUE)

#this contains the most important meta data
#65539 cells, 6 levels of details, mouse strain, major cell type, then subtypes
Hrvatin2018_celltypelabels <- read.csv('~/test_datasets/Hrvatin2018/GSE102827_cell_type_assignments.csv.gz', header = TRUE)
rownames(Hrvatin2018_celltypelabels) <- Hrvatin2018_celltypelabels$X
Hrvatin2018_celltypelabels <- Hrvatin2018_celltypelabels[,c(2:dim(Hrvatin2018_celltypelabels)[2])]

filtered_cell_labels <- Hrvatin2018_celltypelabels[ Hrvatin2018_celltypelabels$X %in% colnames(Hrvatin2018_counts),]

#jeager
lacar2016_meta <- read.csv('Lacar2016_GSE77067/SraRunTable.txt', header = TRUE)
lacar2016_snHC_counts <- read.table('Lacar2016_GSE77067/GSE77067_sn_hc_counts.txt.gz')
lacar2016_snNE_counts <- read.table('Lacar2016_GSE77067/GSE77067_sn_ne_counts.txt.gz')
lacar2016_wc_counts <- read.table('Lacar2016_GSE77067/GSE77067_wc_counts.txt.gz')

#Loading Jeager data
#Jeagers meta rows are a little out of order wrt their counts, i.e. rows do no correspond to cells order we fix that in a bit
jeager2018_counts <- bind_cols(read.table('Jeager2018_GSE98679/GSE98679_count.txt.gz', header = TRUE, check.names = FALSE),
                               read.table('Jeager2018_GSE98679/GSE98679_v2_GSM3308862-GSM3309413_count.txt.gz', header = TRUE, check.names = FALSE))

jeager2018_meta <- read.csv('Jeager2018_GSE98679/SraRunTable.txt', header = TRUE)
jeager2018_meta = jeager2018_meta[c(1:46,599:912,47:598),] #here we fix the order
rownames(jeager2018_meta) <- c(1:912)

jeager2018_meta$predicted_cell_type <- as.character(lapply(jeager2018_meta$predicted_cell_type, function(x) if (x=="") {"DG"} else {x}))
jeager2018_meta$predicted_cell_type <- lapply(jeager2018_meta$predicted_cell_type, function(x) if (x=="") {"DG"} else {x})

jeager2018_meta$fos_status <- as.factor(sapply(as.character(jeager2018_meta$source_name), function(y) if (grepl("_F_", y, fixed=TRUE)) "Fos+" else "Fos-"  ))
jeager2018_meta$fos_status[361:dim(jeager2018_meta)[1]] <- "Fos+"

jeager2018_meta$Mouse_Number[c(361:912)] <- jeager2018_meta$mousingle_number[c(361:912)]
jeager2018_meta <- jeager2018_meta %>% 
  dplyr::select(-mousingle_number)

#filtering out cells as per instructions in Jeager et al., 2018)
under4k <- sapply(jeager2018_counts, function(y) sum(length(which(y>0))))

# as per the mthods section we remove cells with less than 4k gene expressed or 100k reads aligned
filtered.idx <- as.numeric(which((under4k>4000)&(jeager2018_meta$alignable>100000)))
filtered.idx <- order(c(filtered.idx,194))

jeager2018_counts <-jeager2018_counts[,filtered.idx]
jeager2018_counts[is.na(jeager2018_counts)] <- 0

jeager2018_meta <- jeager2018_meta[filtered.idx,]

# MERGING THE DATASETS 
#matching all genes and merging datasets
multi.intersect <- function(x) Reduce(intersect, x) #takes lists of lists, c() will not work

shared.genes <- multi.intersect(list(rownames(jeager2018_counts),
                                     rownames(lacar2016_wc_counts),
                                     rownames(lacar2016_snHC_counts),
                                     rownames(lacar2016_snNE_counts),
                                     rownames(Hrvatin2018_counts) )#closing list 
                                )#closing multi.intersect

jeager2018_counts$gene <- as.character(rownames(jeager2018_counts))
jeager2018_counts <- jeager2018_counts[rownames(jeager2018_counts) %in% shared.genes,]
#jeager2018_counts <- lastcol.to.firstcol(jeager2018_counts)

lacar2016_wc_counts$gene <- as.character(rownames(lacar2016_wc_counts))
lacar2016_wc_counts <- lacar2016_wc_counts[rownames(lacar2016_wc_counts) %in% shared.genes,]

#lacar2016_wc_counts <- lastcol.to.firstcol(lacar2016_wc_counts)

lacar2016_snHC_counts$gene <- as.character(rownames(lacar2016_snHC_counts))
lacar2016_snHC_counts <-lacar2016_snHC_counts[rownames(lacar2016_snHC_counts) %in% shared.genes,]
#lacar2016_snHC_counts <- lastcol.to.firstcol(lacar2016_snHC_counts)

lacar2016_snNE_counts$gene <- as.character(rownames(lacar2016_snNE_counts))
lacar2016_snNE_counts <-lacar2016_snNE_counts[rownames(lacar2016_snNE_counts) %in% shared.genes,]
#lacar2016_snNE_counts <- lastcol.to.firstcol(lacar2016_snNE_counts)

Hrvatin2018_counts$gene <- as.character(rownames(Hrvatin2018_counts))
Hrvatin2018_counts <- Hrvatin2018_counts[rownames(Hrvatin2018_counts) %in% shared.genes,]

# we will remove the PTZ treated cells as well before matching its genes
not.ptz <- which(lacar2016_meta$treatment != "PTZ")


#Match the gene sets by adding the gene names as rows, will strip later
DG.idx <- which(jeager2018_meta$predicted_cell_type=="DG")

# Hrvating celltypes
# > names(table(Hrvatin2018_celltypelabels$celltype)) 
# [1] "Astro"      "Endo_1"     "Endo_2"     "ExcL23"     "ExcL4"     
# [6] "ExcL5_1"    "ExcL5_2"    "ExcL5_3"    "ExcL6"      "Hip"       
# [11] "Int_Cck"    "Int_Npy"    "Int_Pv"     "Int_Sst_1"  "Int_Sst_2" 
# [16] "Int_Vip"    "Macrophage" "Micro_1"    "Micro_2"    "Olig_1"    
# [21] "Olig_2"     "Olig_3"     "Olig_4"     "Olig_5"     "Olig_6"    
# [26] "Olig_7"     "OPC_1"      "OPC_2"      "Pericyte"   "RSP"       
# [31] "SM_1"       "SM_2"       "Sub"


hrvatin.neuron.types <- c("ExcL23", "ExcL4", "ExcL5_1", "ExcL5_2", "ExcL5_3", "ExcL6",       
                          "Int_Cck", "Int_Npy", "Int_Pv", "Int_Sst_1", "Int_Sst_2", 
                          "Int_Vip")
hrvatin.neuron.idx <- which(Hrvatin2018_celltypelabels$celltype %in% hrvatin.neuron.types)

# > table(Hrvatin2018_celltypelabels$stim[hrvatin.neuron.idx])
# 
# 0h   1h   4h 
# 5470 3779 4917



# we must add 1 to the values of DG.idx and not.ptz to deal with the generow shifting the index 1
combined.counts <- jeager2018_counts[, c(DG.idx,862)] %>% 
  left_join(lacar2016_wc_counts[, c(not.ptz[not.ptz <= 82],83)], by = 'gene' , all.y = TRUE) %>%
  left_join(lacar2016_snHC_counts, by = 'gene', all.y = TRUE) %>%
  left_join(lacar2016_snNE_counts, by = 'gene', all.y = TRUE) %>%
  left_join(Hrvatin2018_counts[,c(hrvatin.neuron.idx, 65540)], by = 'gene', all.y = TRUE) 


#this join is possibly including the gene rows leading to mismatch number of cells later 

#give the combined.counts genes for rownames and get rid of that column
rownames(combined.counts) <- combined.counts$gene
combined.counts$gene <- NULL
combined.counts[is.na(combined.counts)] <- 0

#cleaning up the gene column from the other count data
jeager2018_counts$gene <- NULL
lacar2016_wc_counts$gene <- NULL
lacar2016_snHC_counts$gene <- NULL
lacar2016_snNE_counts$gene <- NULL
Hrvatin2018_counts$gene <- NULL

###   MAKING META-DATA FOR COMBINED COUTNS
#columns for which paper the cells are from
experiment.label <- c(replicate( length(DG.idx),"Jeager" ),
                      replicate( length(not.ptz),"Lacar" ))


#column for treatment
treatment <- c(jeager2018_meta$exposure[DG.idx],
               lacar2016_meta$treatment[not.ptz])
treatment <- as.character(lapply(treatment, function(x) if (x=="home-cage") {"HC"} else {x}))
treatment <- as.character(lapply(treatment, function(x) if (x=="novel environment") {"NE"} else {x}))

#fos status 
fos_status <-c(as.character(jeager2018_meta$fos_status[DG.idx]),
               lacar2016_meta$facs_sort[not.ptz])
fos_status <- as.character(lapply(fos_status, function(x) if (x=="Prox1+/Fos+") {"Fos+"} else {x}))
fos_status <- as.character(lapply(fos_status, function(x) if (x=="NeuN+/Prox1+") {"Fos-"} else {x}))
fos_status <- as.character(lapply(fos_status, function(x) if (x=="Prox1+/Fos-") {"Fos-"} else {x}))
fos_status <- as.character(lapply(fos_status, function(x) if (x=="GFP+") {"Fos-"} else {x}))
# need to check what GFP + means, I think it means it is fos positive


combined.meta <- data.frame(experiment.label,
                            treatment,
                            fos_status)
#this throws an error mismathc number of rows and genes most likely
rownames(combined.meta) <- colnames(combined.counts)




# functions (should consider making a module to call)
multi.intersect <- function(x) Reduce(intersect, x) #takes lists of lists, c() will not work
