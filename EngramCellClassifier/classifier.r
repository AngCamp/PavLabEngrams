#
#
#
#
#
# first attempt at a classifier
setwd("C:/Users/angus/Desktop/PavLabEngrams/EngramCellClassifier")

library(tidyverse)
library(GEOquery)
library(AnnotationDbi)


#jeager et al., (2018) meta data and counts
jeager2018_counts <- read.table('Jeager2018_GSE98679/GSE98679_count.txt.gz', header = TRUE, check.names = FALSE)
jeager2018_v2_counts <- read.table('Jeager2018_GSE98679/GSE98679_v2_GSM3308862-GSM3309413_count.txt.gz', header = TRUE, check.names = FALSE)

jeager2018_meta <- read.csv('Jeager2018_GSE98679/SraRunTable.txt', header = TRUE)


#Lacar et al., (2016)
lacar2016_meta <- read.csv('Lacar2016_GSE77067/SraRunTable.txt', header = TRUE)
lacar2016_snHC_counts <- read.table('Lacar2016_GSE77067/GSE77067_sn_hc_counts.txt.gz')
lacar2016_snNE_counts <- read.table('Lacar2016_GSE77067/GSE77067_sn_ne_counts.txt.gz')



