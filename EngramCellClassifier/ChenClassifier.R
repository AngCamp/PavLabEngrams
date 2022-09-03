# Chen Classification



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
library(caret)


#Loading data

setwd("/home/acampbell/PavLabEngrams/EngramCellClassifier")

#Load the data and Metadata
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

# comes from ClusteringCellsChen.R in chen folder
cellmarkers.col <- read.table("Chen2020_GSE152632/Chen2020_cellclusterlabels.csv")


#adding cell bacrcodes from coutn data to rows of metadata for seurat
rownames(chen2020_meta) <- colnames(chen2020_counts)


# our first attempt will focus on within cell types I will try my normalization 
# done within cell types

#normalization functions, log.norm calls logplusone
# mean center scale then log(x+1) for normalizing
logplusone <- function(x){
  return( log(x+1) )
}

log.norm <- function(df.in){
  #performs the lognorm transform and scales the data, removes NA's first
  if( sum(is.na(df.in)) ){
    df.in[is.na(df.in)] <- 0
  }
  df.out <- apply(df.in,
                  MARGIN = 1,
                  FUN = logplusone
  )
  df.out <- scale( t(df.out) ) #we need it transposed so that the scaling is done per gene not cell
  df.out <- data.frame( t(df.out) )
  colnames(df.out) <- rownames(df.in)
  return(df.out)
}


celltype.lognorm <-function(countsdata, celltype.labels){
  #log normalizes within cell types in counts data
  #celltype labels and colnames of countsdata must have same order
  
  #retunrs a transposed and normalize dataframe 
  
  print("Normalizing cell type...")
  
  celltypes <- unique(celltype.labels)
  df.out <- data.frame(gene = rownames(countsdata))
  #df.out <- t(df.out)
  #colnames(df.out) <- rownames(countsdata)
  
  cell_names <- colnames(countsdata) # keep this for reorganizing later
  df.out.rownames <- c()
  for(type in celltypes){
    print(type[1])
    normalized.within.type <- log.norm(countsdata[,celltype.labels==type])
    normalized.within.type <- t(normalized.within.type ) # lognomr flips its data
    normalized.within.type <- data.frame(normalized.within.type)
    normalized.within.type <- rownames_to_column(normalized.within.type, var ="gene")
    df.out <- left_join( df.out, normalized.within.type, by = 'gene' )
  }
  
  df.out <- df.out[,2:dim(df.out)[2]] # drops gene column
  df.out <- df.out[cell_names] # to keep original order
  df.out <- t(df.out)
  return( data.frame(df.out) ) 
}


# for making controls
gene.shuffle <-function(dat){
  #shuffles the genes valus within cells,
  #intended to be applied before normalization
  for ( cell in c(1:ncol(dat)) ){
    rand <- sample(nrow(dat)) # generates random vector of numbers from colls
    dat[,cell] <- dat[rand,cell] # shuffles the genes within cells
    
  }# end of loop over cells
  return(dat)
}


# balanced cross validation random forest

balanced.crossvalidated.rf <- function(){
  
}





