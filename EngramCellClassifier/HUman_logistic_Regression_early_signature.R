## logistic regression attempt




## Libraries
library(randomForest)
library(rfUtilities)
library(Seurat)
library(stringr)
library(sampler)
library(caTools)
library(pROC)
library(ggplot2)
library(stats)
library(Dict)
library(pheatmap)
library(caret)
library(data.table)
library(dplyr)
library(readr) # for cross validated logistic regression

#Ayhan 2018
# Ayhan, F., Kulkarni, A., Berto, S., Sivaprakasam, K., Douglas, C., Lega, B. C., &
# Konopka, G. (2021). Resolving cellular and molecular diversity along the 
# hippocampal anterior-to-posterior axis in humans. Neuron, 109(13), 2091-2105.
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8273123/
# 

# Loading data


setwd("/home/acampbell/PavLabEngrams/EngramCellClassifier")
lacar2016_meta <- read.csv('Lacar2016_GSE77067/SraRunTable.txt', header = TRUE)
lacar2016_snHC_counts <- read.table('Lacar2016_GSE77067/GSE77067_sn_hc_counts.txt.gz')
lacar2016_snHC_counts[is.na(lacar2016_snHC_counts)] <- 0
lacar2016_snNE_counts <- read.table('Lacar2016_GSE77067/GSE77067_sn_ne_counts.txt.gz')
lacar2016_snNE_counts[is.na(lacar2016_snNE_counts)] <- 0
lacar2016_wc_counts <- read.table('Lacar2016_GSE77067/GSE77067_wc_counts.txt.gz')
lacar2016_wc_counts[is.na(lacar2016_wc_counts)] <- 0

#Loading Jeager data
#Jeagers meta rows are a little out of order w.r.t. their count data, i.e. rows do no correspond to cells order we fix that in a bit
jeager2018_counts <- bind_cols(read.table('Jeager2018_GSE98679/GSE98679_count.txt.gz', header = TRUE, check.names = FALSE),
                               read.table('Jeager2018_GSE98679/GSE98679_v2_GSM3308862-GSM3309413_count.txt.gz', header = TRUE, check.names = FALSE))
jeager2018_counts[is.na(jeager2018_counts)] <- 0

jeager2018_meta <- read.csv('Jeager2018_GSE98679/SraRunTable.txt', header = TRUE)
jeager2018_meta$umi <- colSums(jeager2018_counts)

jeager2018_meta <- jeager2018_meta[c(1:46,599:912,47:598),] #here we fix the order
rownames(jeager2018_meta) <- c(1:912)
jeager2018_meta$CellID <- colnames(jeager2018_counts)

arclabels <- read.csv("Jeager2018_GSE98679/Jaeger2018_meta_arclabels.csv", header = TRUE)
jeager2018_meta <- left_join(jeager2018_meta,
                             arclabels %>% dplyr::select(Title, ArcStatus),
                             by = c("CellID" = "Title"))




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
                                     rownames(lacar2016_snNE_counts)
)#closing list 
)#closing multi.intersect

jeager2018_counts$gene <- as.character(rownames(jeager2018_counts))
jeager2018_counts <- jeager2018_counts[rownames(jeager2018_counts) %in% shared.genes,]
#jeager2018_counts <- lastcol.to.firstcol(jeager2018_counts)

wc_umi <- colSums(lacar2016_wc_counts)
lacar2016_wc_counts$gene <- as.character(rownames(lacar2016_wc_counts))
lacar2016_wc_counts <- lacar2016_wc_counts[rownames(lacar2016_wc_counts) %in% shared.genes,]

#lacar2016_wc_counts <- lastcol.to.firstcol(lacar2016_wc_counts)

snHC_umi <- colSums(lacar2016_snHC_counts)
lacar2016_snHC_counts$gene <- as.character(rownames(lacar2016_snHC_counts))
lacar2016_snHC_counts <-lacar2016_snHC_counts[rownames(lacar2016_snHC_counts) %in% shared.genes,]
#lacar2016_snHC_counts <- lastcol.to.firstcol(lacar2016_snHC_counts)

snNE_umi <- colSums(lacar2016_snNE_counts)
lacar2016_snNE_counts$gene <- as.character(rownames(lacar2016_snNE_counts))
lacar2016_snNE_counts <-lacar2016_snNE_counts[rownames(lacar2016_snNE_counts) %in% shared.genes,]
#lacar2016_snNE_counts <- lastcol.to.firstcol(lacar2016_snNE_counts)

# we will remove the PTZ treated cells as well before matching its genes
not.ptz <- which(lacar2016_meta$treatment != "PTZ")


#Match the gene sets by adding the gene names as rows, will strip later
DG.idx <- which(jeager2018_meta$predicted_cell_type=="DG")

# we must add 1 to the values of DG.idx and not.ptz to deal with the generow shifting the index 1
combined.counts <- jeager2018_counts[, c(DG.idx,which(colnames(jeager2018_counts)=='gene'))] %>% 
  left_join(lacar2016_wc_counts[, c(not.ptz[not.ptz <= 82],which(colnames(lacar2016_wc_counts)=='gene'))], by = 'gene' , all.y = TRUE) %>%
  left_join(lacar2016_snHC_counts, by = 'gene', all.y = TRUE) %>%
  left_join(lacar2016_snNE_counts, by = 'gene', all.y = TRUE) #%>% 

# making the umi counts
combined_umi <- c(jeager2018_meta$umi[DG.idx],wc_umi[not.ptz[not.ptz <= 82]],
                  snHC_umi, snNE_umi  )
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
facs_sort <-c(as.character(jeager2018_meta$fos_status[DG.idx]),
              lacar2016_meta$facs_sort[not.ptz])
fos_status <- facs_sort # make fos_status from 
fos_status <-c(as.character(jeager2018_meta$fos_status[DG.idx]),
               lacar2016_meta$facs_sort[not.ptz])
fos_status <- as.character(lapply(fos_status, function(x) if (x=="Prox1+/Fos+") {"Fos+"} else {x}))
fos_status <- as.character(lapply(fos_status, function(x) if (x=="NeuN+/Prox1+") {"Fos-"} else {x}))
fos_status <- as.character(lapply(fos_status, function(x) if (x=="Prox1+/Fos-") {"Fos-"} else {x}))
fos_status <- as.character(lapply(fos_status, function(x) if (x=="GFP+") {"Fos-"} else {x}))

ActivityStatus <- fos_status
ActivityStatus <- as.character(lapply(ActivityStatus, function(x) if (x=="Fos-") {"Inactive"} else {x}))
ActivityStatus <- as.character(lapply(ActivityStatus, function(x) if (x=="Fos+") {"Active"} else {x}))
# need to check what GFP + means, I think it means it is fos positive


combined.meta <- data.frame(experiment.label,
                            treatment,
                            facs_sort,
                            fos_status,
                            ActivityStatus,
                            combined_umi)

#this throws an error mismatch number of rows and genes most likely
rownames(combined.meta) <- colnames(combined.counts)
combined.meta$CellID <- rownames(combined.meta)

arclabels <- read.csv("Jeager2018_GSE98679/Jaeger2018_meta_arclabels.csv", header = TRUE)

combined.meta <- left_join(combined.meta,
                           arclabels %>% dplyr::select(Title, ArcStatus),
                           by = c("CellID" = "Title"))

combined.meta <- left_join(combined.meta,
                           arclabels %>% dplyr::select(Title, FosStatus),
                           by = c("CellID" = "Title"))

# I am unsure if the unlablled 5hr cells are arc+ or not
#combined.meta$ArcStatus[is.na(combined.meta$ArcStatus)&combined.meta$treatment=="5hr"] <- "pos"
combined.meta$ArcStatus[is.na(combined.meta$ArcStatus)&combined.meta$treatment=="4hr"] <- "pos"
#combined.meta$ArcStatus[is.na(combined.meta$ArcStatus)] <- "Unknown_or_neg"
# the reason its unknown or negative is that in some rows they do not
# even stain for it for instance at 1hr points or in lacar
# at the 5hr time points they only included 

#write.csv(combined.counts, "jeager_lacar_combinedcounts.csv")
#write.csv(combined.meta, "jeager_lacar_combinedmeta.csv")



#please note that really fos- and arc+ cells should be considered
# sencondary response, this part still needs some work as the labels are kinda messed up
combined.meta$Activity_class <- rep("fillme", dim(combined.meta)[1])
combined.meta$Activity_class[combined.meta$ActivityStatus=="Inactive"] <- "Inactive"
combined.meta$Activity_class[combined.meta$FosStatus=="pos" 
                             & combined.meta$ArcStatus=="pos"] <- "Reactivated"
combined.meta$Activity_class[combined.meta$ActivityStatus=="Active" 
                             & combined.meta$treatment=="1hr"] <- "EarlySignature"
combined.meta$Activity_class[combined.meta$ActivityStatus=="Active" 
                             & combined.meta$treatment=="NE"] <- "EarlySignature"
combined.meta$Activity_class[combined.meta$treatment=="4hr"|
                               combined.meta$treatment=="5hr"|
                               combined.meta$treatment=="A>C"] <- "LateSignature"
combined.meta$Activity_class[combined.meta$treatment=="A>A"&
                               combined.meta$FosStatus=="neg"] <- "LateSignature"
#HOCHGERNER DATA
testsetpath <- "/home/acampbell/test_datasets" # needs to be changed for pavlab server

hochgerner5k_2018_counts <- read.table(paste(testsetpath,"/Hochgerner2018/GSE95315_10X_expression_data_v2.tab.gz", sep=""))

colnames(hochgerner5k_2018_counts) <- hochgerner5k_2018_counts[1,]
rownames(hochgerner5k_2018_counts) <- hochgerner5k_2018_counts[,1]

hochgerner5k_2018_meta <- hochgerner5k_2018_counts %>% 
  dplyr::slice(c(1:3)) %>%
  t() %>%
  data.frame %>%
  dplyr::slice(-1) %>% 
  dplyr::select(-cellid)

hochgerner5k_2018_counts <- hochgerner5k_2018_counts %>% 
  dplyr::select(-cellid) %>% 
  dplyr::slice(-c(1:3))
hochgerner5k_2018_counts[is.na(hochgerner5k_2018_counts)] <- 0

# changing to numeric removing nas and calculating umi column for metadata for nromalizating later
hochcells <- colnames(hochgerner5k_2018_counts)
hochgenes <- rownames(hochgerner5k_2018_counts)
hochgerner5k_2018_counts <- data.frame(lapply(hochgerner5k_2018_counts, as.numeric))
colnames(hochgerner5k_2018_counts) <- hochcells
rownames(hochgerner5k_2018_counts) <- hochgenes

hochgerner5k_2018_counts[is.na(hochgerner5k_2018_counts)] <- 0
# getting umi
hochgerner5k_2018_meta$umi <- colSums(hochgerner5k_2018_counts)

#get index locations of adult p35 DGCs
hoch5k.GC_Adult.p35.idx <- (hochgerner5k_2018_meta$age.days.=="35") | (hochgerner5k_2018_meta$age.days.=="35*")
hoch5k.GC_Adult.p35.idx <- (hoch5k.GC_Adult.p35.idx) & (hochgerner5k_2018_meta$cluster_name == "Granule-mature")
hoch5k.GC_Adult.p35.idx <- which(hoch5k.GC_Adult.p35.idx)


hochgerner5k.DGC.idx <- hochgerner5k_2018_meta$cluster_name == "Granule-mature"
# we use this later
hochgernerDGC_counts <-hochgerner5k_2018_counts[, hochgerner5k.DGC.idx]


#### Functions

# paul prefers a log base that's easy to do headmath with so no eulers numebr
# normalization functions
pseudocount_log2p1_transform <- function(x, scale_factor = 10^6, UMI.provided = NULL){
  if(is.null(UMI.provided)){
    counts <- sum(x)}else{
      counts <- UMI.provided
    }
  x <- (x+1)/counts
  x <- x/scale_factor
  return(log2(x))
}

pavlab.normalize <- function(df, UMI = NULL){
  df.cols <- colnames(df)
  df.rows <- rownames(df)
  if( is.null(UMI)){
    df <- data.frame(apply(df,  MARGIN = 2, pseudocount_log2p1_transform))
  }else{
    #
    df[] <- Map(pseudocount_log2p1_transform, df, UMI.provided = UMI)
    
  }
  colnames(df) <- df.cols
  rownames(df)<- df.rows
  return(df)
}


#normalization functions, log.norm calls logplusone
seurat_log1p_transform <- function(x, scale_factor = 10000, UMI.provided = NULL){
  if(is.null(UMI.provided)){
    counts <- sum(x)}else{
      counts <- UMI.provided
    }
  x <- (x+1)/counts
  x <- x/scale_factor
  return(log(x))
}

seurat.normalize <- function(df, UMI = NULL){
  df.cols <- colnames(df)
  df.rows <- rownames(df)
  if( is.null(UMI)){
    df <- data.frame(apply(df,  MARGIN = 2, seurat_log1p_transform))
  }else{
    #
    df[] <- Map(seurat_log1p_transform, df, UMI.provided = UMI)
  }
  colnames(df) <- df.cols
  rownames(df)<- df.rows
  return(df)
}


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


resample.randomForest <-function( df.in,
                                  under_represented_class,
                                  over_represented_class,
                                  proportion,
                                  batches, 
                                  trees){
  #NOTE: df.in should have a column called engram cell with the class labels i.e. postive or negative
  
  #this function resamples from our samples and retrains new models then combines them
  # this is too prevent over fitting on cells
  trees.per.batch <- as.integer(trees/batches)
  n.cells <- trunc( sum(df.in$Engramcell==under_represented_class)*proportion)
  batches <- c(1:batches)
  for( batch in batches){
    resample.set <- rbind(sample(which(df.in$Engramcell==under_represented_class), size = n.cells),
                          sample(which(df.in$Engramcell==over_represented_class), size = n.cells)
    )
    resample.set <- df.in[resample.set,]
    
    # creates rf.model
    if(batch==1){
      rf.model <- randomForest(x = resample.set[,1:(length(resample.set)-1)],
                               y = resample.set$Engramcell,
                               ntree = trees.per.batch)
    }
    #trains new models in rf.fit and combines tham with rf.model
    if(batch>1){
      rf.fit = randomForest(x = resample.set[,1:(length(resample.set)-1)],
                            y = resample.set$Engramcell,
                            ntree = trees.per.batch)
      rf.model <- randomForest::combine(rf.fit, rf.model)
    }
  }#end of for loop over batches
  
  return(rf.model)
}


make.predictions.df <- function(classifier.object, 
                                test_df,
                                meta.data.label.column,
                                label = c("Active","Inactive")
){
  #generate predictions for making classifier summary
  predictions <- as.data.frame(predict(classifier.object, test_df[,1:(length(test_df))], type = "prob"))
  predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)] #1:2 for the number of classes
  predictions$observed <- meta.data.label.column #this should be changed if you want to make this functions more modular
  colnames(predictions)[1:2] <- c("label_pos","label_neg")
  predictions$engramobserved <- ifelse(predictions$observed==label[1], 1, 0)
  predictions$inactiveobserved <- ifelse(predictions$observed==label[2], 1, 0)
  return(predictions)
}


# assess a single run of resampled.randomforest
assessment <- function(predictions.df, 
                       label = c("Active","Inactive") 
){
  # returns a vector of assessments to be used to make dataframe summarizing classifiers performance
  # can be used to make df of all classifiers trained in a single run
  TP <- sum((predictions.df$predict == label[1])&(predictions.df$observed == label[1]))
  TN <- sum((predictions.df$predict == label[2])&(predictions.df$observed == label[2]))
  FN <- sum((predictions.df$predict == label[2])&(predictions.df$observed == label[1]))
  FP <- sum((predictions.df$predict == label[1])&(predictions.df$observed == label[2]))
  
  #precision and recall as well as sumamry stats F1
  precision <- TP/(TP+FP)
  recall <- TP/(TP+FN)
  F1.score = 2 * (precision * recall) / (precision + recall)
  FPR <- FP/(TN+FP)
  FNR <- FN/(TP+FN)
  
  #getting auc
  roc.engramcell <- roc(predictions.df$engramobserved, as.numeric(predictions.df$label_pos) )
  AUC <- auc(roc.engramcell)
  
  return( c(F1.score, AUC, precision, recall, FPR, FNR,
            TP, FN, TN, FP) )
}



resampled.randomForest.crossvalidated <-function(data,
                                                 under.represented.class,
                                                 over.represented.class,
                                                 folds,
                                                 trees.total,
                                                 proportion.each.batch=0.8,
                                                 batches.per.fold=20){
  # takes a data frame with a label column assumed to be named Engramcell, data$Engramcell
  # returns a model that has been k-fold cross validated, with an attribute called Assessment
  # assessment has the performance metrics of all the folds and a column of means and SD's for each
  # metric
  #NOTE: ROC curve needs to be implemented
  
  folds.obj <- createFolds(data$Engramcell, k = folds)
  loops <- c(1:folds)
  for( i in loops ){
    #create indices
    test.idx <- folds.obj[[i]]
    # needs to be a list so it can act as an index
    train.idx <- which(!(rownames(data) %in% test.idx) )
    
    #split data for this fold
    training_set <- data[train.idx,]
    testing_set <- data[test.idx,]
    
    # divvies up number of trees
    trees.in.the.fold = as.integer(trees.total/folds)
    if ( ( trees.total%%(batches.per.fold*folds) )>0  ){ 
      stop("Number of trees does not devide evenly by batches and folds.")
    }
    # we still need to settle on stuff to 
    rf.this_fold <- resample.randomForest(df.in = training_set,
                                          under_represented_class = under.represented.class,
                                          over_represented_class = over.represented.class,
                                          proportion= proportion.each.batch,
                                          batches = batches.per.fold, 
                                          trees = trees.in.the.fold)
    
    if(i == 1){
      rf.out <- rf.this_fold
      pred <- make.predictions.df(rf.this_fold, testing_set[1:(length(testing_set)-1)], testing_set$Engramcell)
      assess <- assessment( pred ) 
      fold.performance <- data.frame(assess )
      rownames(fold.performance) <- c("F1 Score", "AUC", "Precision", "Recall",
                                      "FPR", "FNR", "True Positives", "False Negatives", 
                                      "True Negatives", "False Positives")
    }else{
      rf.out <- randomForest::combine(rf.out, rf.this_fold)
      # we need votes for all cells to calculate
      pred <- make.predictions.df(rf.this_fold, testing_set[1:(length(testing_set)-1)], testing_set$Engramcell)
      assess <- assessment( pred ) 
      fold.performance[,ncol(fold.performance) + 1] <- assess
    }
    
  }# end of for loop
  colnames(fold.performance) <- names(folds.obj)
  fold.performance$Mean <- apply(fold.performance,MARGIN=1,  FUN = mean)
  fold.performance$SigDiff <- apply(fold.performance,MARGIN=1,  FUN = sd)
  rf.out$Assessment <- fold.performance
  
  #votes needs to be updated to make roc curve
  rf.out$votes <- predict(object = rf.out, newdata = data, type = 'vote', norm.votes = FALSE)
  return(rf.out)
}

resample.regularizedRF <- function( df.in,
                                    under_represented_class,
                                    over_represented_class,
                                    proportion,
                                    batches, 
                                    trees){
  #NOTE: df.in should have a column called engram cell with the class labels i.e. postive or negative
  
  #this function resamples from our samples and retrains new models then combines them
  # this is too prevent over fitting on cells
  trees.per.batch <- as.integer(trees/batches)
  n.cells <- trunc( sum(df.in$Engramcell==under_represented_class)*proportion)
  batches <- c(1:batches)
  for( batch in batches){
    resample.set <- rbind(sample(which(df.in$Engramcell==under_represented_class), size = n.cells),
                          sample(which(df.in$Engramcell==over_represented_class), size = n.cells)
    )
    resample.set <- df.in[resample.set,]
    
    # creates rf.model
    if(batch==1){
      rf.model <- RRF(x = resample.set[,1:(length(resample.set)-1)],
                      y = resample.set$Engramcell,
                      ntree = trees.per.batch)
    }
    #trains new models in rf.fit and combines tham with rf.model
    if(batch>1){
      rf.fit = RRF(x = resample.set[,1:(length(resample.set)-1)],
                   y = resample.set$Engramcell,
                   ntree = trees.per.batch)
      rf.model <- RRF::combine(rf.fit, rf.model)
    }
  }#end of for loop over batches
  
  return(rf.model)
}

#
resampled.regularizedRF.crossvalidated <-function(data,
                                                  under.represented.class,
                                                  over.represented.class,
                                                  folds,
                                                  trees.total,
                                                  proportion.each.batch=0.8,
                                                  batches.per.fold=20){
  # takes a data frame with a label column assumed to be named Engramcell, data$Engramcell
  # returns a model that has been k-fold cross validated, with an attribute called Assessment
  # assessment has the performance metrics of all the folds and a column of means and SD's for each
  # metric
  #NOTE: ROC curve needs to be implemented
  
  folds.obj <- createFolds(data$Engramcell, k = folds)
  loops <- c(1:folds)
  for( i in loops ){
    #create indices
    test.idx <- folds.obj[[i]]
    # needs to be a list so it can act as an index
    train.idx <- which(!(rownames(data) %in% test.idx) )
    
    #split data for this fold
    training_set <- data[train.idx,]
    testing_set <- data[test.idx,]
    
    # divvies up number of trees
    trees.in.the.fold = as.integer(trees.total/folds)
    if ( ( trees.total%%(batches.per.fold*folds) )>0  ){ 
      stop("Number of trees does not devide evenly by batches and folds.")
    }
    # we still need to settle on stuff to 
    rf.this_fold <- resample.regularizedRF(df.in = training_set,
                                           under_represented_class = under.represented.class,
                                           over_represented_class = over.represented.class,
                                           proportion= proportion.each.batch,
                                           batches = batches.per.fold, 
                                           trees = trees.in.the.fold)
    
    if(i == 1){
      rf.out <- rf.this_fold
      pred <- make.predictions.df(rf.this_fold, testing_set[1:(length(testing_set)-1)], testing_set$Engramcell)
      assess <- assessment( pred ) 
      fold.performance <- data.frame(assess )
      rownames(fold.performance) <- c("F1 Score", "AUC", "Precision", "Recall",
                                      "FPR", "FNR", "True Positives", "False Negatives", 
                                      "True Negatives", "False Positives")
    }else{
      rf.out <- RRF::combine(rf.out, rf.this_fold)
      # we need votes for all cells to calculate
      pred <- make.predictions.df(rf.this_fold, testing_set[1:(length(testing_set)-1)], testing_set$Engramcell)
      assess <- assessment( pred ) 
      fold.performance[,ncol(fold.performance) + 1] <- assess
    }
    
  }# end of for loop
  colnames(fold.performance) <- names(folds.obj)
  fold.performance$Mean <- apply(fold.performance,MARGIN=1,  FUN = mean)
  fold.performance$SigDiff <- apply(fold.performance,MARGIN=1,  FUN = sd)
  rf.out$Assessment <- fold.performance
  
  #votes needs to be updated to make roc curve
  rf.out$votes <- predict(object = rf.out, newdata = data, type = 'vote', norm.votes = FALSE)
  return(rf.out)
}





gene.shuffle <-function(dat){
  #shuffles the genes valus within cells,
  #intended to be applied before normalization
  for ( cell in c(1:ncol(dat)) ){
    rand <- sample(nrow(dat)) # generates random vector of numbers from colls
    dat[,cell] <- dat[rand,cell] # shuffles the genes within cells
    
  }# end of loop over cells
  return(dat)
}





#####
###  FINDING CELLS IN OTHER MICE  (HOCHGERNER) -gene merging takes place here
####


# Actually running the code with cross validation

labeled.data.hochgerner2018_genes <-  pavlab.normalize(combined.counts,
                                                       UMI= combined.meta$combined_umi)
genes <- rownames(labeled.data.hochgerner2018_genes)
labeled.data.hochgerner2018_genes <- transpose(labeled.data.hochgerner2018_genes)
colnames(labeled.data.hochgerner2018_genes) <- genes 

hoch5k.adultDGCs.lognorm <- pavlab.normalize(hochgerner5k_2018_counts[, hoch5k.GC_Adult.p35.idx],
                                             UMI = hochgerner5k_2018_meta$umi[hoch5k.GC_Adult.p35.idx])
genes <- rownames(hoch5k.adultDGCs.lognorm)
hoch5k.adultDGCs.lognorm <- transpose(hoch5k.adultDGCs.lognorm)
colnames(hoch5k.adultDGCs.lognorm) <- genes

#gene matching, there were extra commas here, possibly the cause of the errors
shared.genes <- intersect( colnames(hoch5k.adultDGCs.lognorm), colnames(labeled.data.hochgerner2018_genes) )
hoch5k.adultDGCs.lognorm  <- hoch5k.adultDGCs.lognorm[,colnames(hoch5k.adultDGCs.lognorm) %in% shared.genes] # 
labeled.data.hochgerner2018_genes <- labeled.data.hochgerner2018_genes[,colnames(labeled.data.hochgerner2018_genes) %in% shared.genes]

labeled.data.hochgerner2018_genes$Engramcell <- as.factor(combined.meta$ActivityStatus)
labeled.data.hochgerner2018_genes$Engramcell <- as.factor(labeled.data.hochgerner2018_genes$Engramcell=='Active')

### logistic regression calssifier




ctrlspecs <- trainControl(method="cv", 
                          number=10, 
                          savePredictions="all",
                          classProbs=TRUE)
set.seed(1985)

# Specify logistic regression model to be estimated using training data
# and k-fold cross-validation process
model1 <- train(Engramcell ~ ., data=labeled.data.hochgerner2018_genes, 
                method="glm", 
                family=binomial, 
                trControl=ctrlspecs)



# AUC plot 
# https://stackoverflow.com/questions/69246553/plot-the-average-cross-validated-auc-from-caret-package


# predicting new data
predict(gbmFit1, newdata = head(testing), type = "prob")







# note that the third kind of DGCs they found were ignored
ayhanDGC_counts <- read.csv("~/test_datasets/Ayhan2021_GSE160189/ayhanDGC_counts.csv", header = TRUE)
rownames(ayhanDGC_counts) <- ayhanDGC_counts$gene
# the OG data has na but we switched them to zeros
ayhanDGC_counts <- ayhanDGC_counts[, c(1:dim(ayhanDGC_counts)[2]-1)]
ayhanDGC_meta <- read.csv("~/test_datasets/Ayhan2021_GSE160189/ayhanDGC_meta.csv", header = TRUE)

# Orthologue matching
# from alex's ortholog mapping

hg_to_mm <- read.table("hg_mm_1to1_ortho_genes_DIOPT-v8.tsv", sep = '\t', header = TRUE)
# just run left join onto the appropriate column for each dataset


present.orthologs.idx <- (hg_to_mm$Symbol_hg %in% rownames(ayhanDGC_counts))&(hg_to_mm$Symbol_mm %in% rownames(combined.counts))&(hg_to_mm$Symbol_mm %in% rownames(hochgernerDGC_counts) )
hg_to_mm <- hg_to_mm[present.orthologs.idx,] # filter for matches, ~ 14403 present between Ayhan and Jeager/Lacar

# filter using left join and swithc out mm for hg in jeager/lacar data
combined.counts.hg <- combined.counts
combined.counts.hg$Symbol_mm <- rownames(combined.counts.hg)

# filtering for orthologs present in both data sets
combined.counts.hg <- left_join(x = hg_to_mm, y = combined.counts.hg, by = "Symbol_mm" )
rownames(combined.counts.hg) <- combined.counts.hg$Symbol_hg
combined.counts.hg <- combined.counts.hg[,c(4:(dim(combined.counts.hg)[2]) )]

ayhanDGC_counts$Symbol_hg <- rownames(ayhanDGC_counts)
ayhanDGC_counts <- left_join(x = hg_to_mm, y = ayhanDGC_counts, by = "Symbol_hg" )
rownames(ayhanDGC_counts) <- ayhanDGC_counts$Symbol_hg
ayhanDGC_counts <- ayhanDGC_counts[,c(4:dim(ayhanDGC_counts)[2])]

hochgernerDGC_counts.hg <- hochgernerDGC_counts
hochgernerDGC_counts.hg$Symbol_mm <- rownames(hochgernerDGC_counts.hg)
hochgernerDGC_counts.hg <- left_join(x = hg_to_mm, y = hochgernerDGC_counts.hg, by = "Symbol_mm" )
rownames(hochgernerDGC_counts.hg) <- hochgernerDGC_counts.hg$Symbol_hg
hochgernerDGC_counts.hg <- hochgernerDGC_counts.hg[,c(4:(dim(hochgernerDGC_counts.hg)[2]) )]

# gene filtering for matching genes



# normalize data, tranpose and add engram label
combined.counts.hg.normed <-  pavlab.normalize(combined.counts.hg,
                                               UMI= combined.meta$combined_umi)
genes <- rownames(combined.counts.hg.normed )
combined.counts.hg.normed  <- transpose(combined.counts.hg.normed)
colnames(combined.counts.hg.normed ) <- genes 
combined.counts.hg.normed$Engramcell <- as.factor(combined.meta$ActivityStatus)
#combined.counts.hg.normed$Engramcell <- as.factor(combined.counts.hg.normed$Engramcell=='Active')
  
# trian model
ctrlspecs <- trainControl(method="cv", 
                          number=5, 
                          savePredictions="all",
                          classProbs=TRUE)
set.seed(1984)

logreg.earlysignature.hg <- train(Engramcell ~ ., data=combined.counts.hg.normed, 
                method="glm", 
                family=binomial, 
                trControl=ctrlspecs)

test <- data.frame(coefs = logreg.earlysignature.hg$finalModel$coefficients, 
                   genes=names(logreg.earlysignature.hg$finalModel$coefficients))

test <- arrange(test, coefs, decreasing = TRUE)






classifier.hochgerner2018_genes <- resampled.randomForest.crossvalidated(data= combined.counts.hg.normed,
                                                                         under.represented.class = "Inactive",
                                                                         over.represented.class = "Active",
                                                                         trees.total = 1000,
                                                                         folds = 5,
                                                                         proportion.each.batch=0.8,
                                                                         batches.per.fold=5)

# Importance
#classifier.hochgerner2018_genes <- classifier.hoch.normal
importance.df.hoch<- data.frame(gene = as.character( rownames(classifier.hochgerner2018_genes$importance) ),
                                importance_score = as.numeric(classifier.hochgerner2018_genes$importance ) ) %>%
  arrange(desc(importance_score))

head(importance.df.hoch, 1000) # don't forget to run this before gettng  the regularized






