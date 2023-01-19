# cleaning Kwon 2021


# Clustering Kwon in seurat

library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(stringr)
library(purrr)

# ## Libraries
# library(randomForest)
# library(rfUtilities)
# library(stringr)
# library(sampler)
# library(caTools)
# library(pROC)
# library(stats)
# library(Dict)
# library(pheatmap)
# library(caret)
# library(data.table)

setwd("~/test_datasets/Kwon2021_GSE145970")


# OG Paper
# Kwon, D. Y., Xu, B., Hu, P., Zhao, Y. T., Beagan, J. A., Nofziger, J. H., ... &
# Zhou, Z. (2022). Neuronal Yin Yang1 in the prefrontal cortex regulates 
# transcriptional and behavioral responses to chronic stress in mice. Nature 
# communications, 13(1), 1-19.
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8748737/

Ctrl_Rep30 <- read.table("GSM4340462_ZhouLab_cortex_Ctrl_Rep30_DGE.txt.gz")
Ctrl_Rep31 <- read.table("GSM4340463_ZhouLab_cortex_Ctrl_Rep31_DGE.txt.gz")
Ctrl_Rep33 <- read.table("GSM4340464_ZhouLab_cortex_Ctrl_Rep33_DGE.txt.gz")
Ctrl_Rep34 <- read.table("GSM4340465_ZhouLab_cortex_Ctrl_Rep34_DGE.txt.gz")
Stress_Rep45 <- read.table("GSM4340466_ZhouLab_cortex_Stress_Rep45_DGE.txt.gz")
Stress_Rep46 <- read.table("GSM4340467_ZhouLab_cortex_Stress_Rep46_DGE.txt.gz")
Stress_Rep48 <- read.table("GSM4340468_ZhouLab_cortex_Stress_Rep48_DGE.txt.gz")
Stress_Rep49 <- read.table("GSM4340469_ZhouLab_cortex_Stress_Rep49_DGE.txt.gz")

change_colnames <- function(df.in, s, idx = c(1:dim(df.in)[2]) ){
  colnames(df.in) <- df.in[1,]
  df.in <- df.in[c(2:dim(df.in)[1]),]
  strlist <- colnames(df.in)
  for(i in idx ){
    strlist[i] <- paste(strlist[i],s)
  }
  colnames(df.in) <- strlist
  return(df.in)
}
cells <- c(2:5001)
Ctrl_Rep30 <- change_colnames(Ctrl_Rep30, "_Ctrl_Rep30", idx = cells)
Ctrl_Rep31 <-change_colnames(Ctrl_Rep31, "_Ctrl_Rep31",idx = cells) 
Ctrl_Rep33 <- change_colnames(Ctrl_Rep33, "_Ctrl_Rep33", idx = cells)
Ctrl_Rep34 <- change_colnames(Ctrl_Rep34, "_Ctrl_Rep34", idx = cells) 
Stress_Rep45 <- change_colnames(Stress_Rep45, "_Stress_Rep45", idx = cells)
Stress_Rep46 <- change_colnames(Stress_Rep46, "_Stress_Rep46", idx = cells)
Stress_Rep48 <- change_colnames(Stress_Rep48, "_Stress_Rep48", idx = cells)
Stress_Rep49 <- change_colnames(Stress_Rep49, "_Stress_Rep49", idx = cells)

dflist <- list(Ctrl_Rep30, Ctrl_Rep31, Ctrl_Rep33, Ctrl_Rep34,
               Stress_Rep45, Stress_Rep46, Stress_Rep48, Stress_Rep49)

kwon2021_counts <- purrr::reduce(dflist, dplyr::left_join, by = 'GENE')

write_csv(kwon2021_counts, "kwon2021_dropseq_counts.csv")