## this will be the seconds draft of the jeager data which will A run properly and B implement cross validation
# I will also include the controls here, this script can be used to generate the figures I will put into
# the final version in R markdown for the results which will just call the figures and present the code as needed

# consider running this: https://datacornering.com/how-to-run-r-script-from-another-r-script-and-use-as-a-source/#:~:text=You%20can%20execute%20R%20script,the%20file%20path%20contains%20space.
# this way we could turn loading each of these datasets into a single line of code
# or create a module that does this

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
library(caret)
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

# we use this later
hochgernerDGC_counts <-hochgerner5k_2018_counts[, hochgerner5k_2018_meta$cluster_name == "Granule-mature"]
hochgernerDGC_counts <- data.frame( lapply(hochgernerDGC_counts,as.numeric) ) # this data is for some reason all strings
colnames(hochgernerDGC_counts) <- colnames(hochgerner5k_2018_counts)[hochgerner5k_2018_meta$cluster_name == "Granule-mature"]

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

#modified version of randomForest::combine() which can handle different sized forests
# without this it throws and error when combining the votes
my_combine <- function (...) 
{
  pad0 <- function(x, len) c(x, rep(0, len - length(x)))
  padm0 <- function(x, len) rbind(x, matrix(0, nrow = len - 
                                              nrow(x), ncol = ncol(x)))
  rflist <- list(...)
  areForest <- sapply(rflist, function(x) inherits(x, "randomForest"))
  if (any(!areForest)) 
    stop("Argument must be a list of randomForest objects")
  rf <- rflist[[1]]
  classRF <- rf$type == "classification"
  trees <- sapply(rflist, function(x) x$ntree)
  ntree <- sum(trees)
  rf$ntree <- ntree
  nforest <- length(rflist)
  haveTest <- !any(sapply(rflist, function(x) is.null(x$test)))
  vlist <- lapply(rflist, function(x) rownames(importance(x)))
  numvars <- sapply(vlist, length)
  if (!all(numvars[1] == numvars[-1])) 
    stop("Unequal number of predictor variables in the randomForest objects.")
  for (i in seq_along(vlist)) {
    if (!all(vlist[[i]] == vlist[[1]])) 
      stop("Predictor variables are different in the randomForest objects.")
  }
  haveForest <- sapply(rflist, function(x) !is.null(x$forest))
  if (all(haveForest)) {
    nrnodes <- max(sapply(rflist, function(x) x$forest$nrnodes))
    rf$forest$nrnodes <- nrnodes
    rf$forest$ndbigtree <- unlist(sapply(rflist, function(x) x$forest$ndbigtree))
    rf$forest$nodestatus <- do.call("cbind", lapply(rflist, 
                                                    function(x) padm0(x$forest$nodestatus, nrnodes)))
    rf$forest$bestvar <- do.call("cbind", lapply(rflist, 
                                                 function(x) padm0(x$forest$bestvar, nrnodes)))
    rf$forest$xbestsplit <- do.call("cbind", lapply(rflist, 
                                                    function(x) padm0(x$forest$xbestsplit, nrnodes)))
    rf$forest$nodepred <- do.call("cbind", lapply(rflist, 
                                                  function(x) padm0(x$forest$nodepred, nrnodes)))
    tree.dim <- dim(rf$forest$treemap)
    if (classRF) {
      rf$forest$treemap <- array(unlist(lapply(rflist, 
                                               function(x) apply(x$forest$treemap, 2:3, pad0, 
                                                                 nrnodes))), c(nrnodes, 2, ntree))
    }
    else {
      rf$forest$leftDaughter <- do.call("cbind", lapply(rflist, 
                                                        function(x) padm0(x$forest$leftDaughter, nrnodes)))
      rf$forest$rightDaughter <- do.call("cbind", lapply(rflist, 
                                                         function(x) padm0(x$forest$rightDaughter, nrnodes)))
    }
    rf$forest$ntree <- ntree
    if (classRF) 
      rf$forest$cutoff <- rflist[[1]]$forest$cutoff
  }
  else {
    rf$forest <- NULL
  }
  #
  #Tons of stuff removed here...
  #
  if (classRF) {
    rf$confusion <- NULL
    rf$err.rate <- NULL
    if (haveTest) {
      rf$test$confusion <- NULL
      rf$err.rate <- NULL
    }
  }
  else {
    rf$mse <- rf$rsq <- NULL
    if (haveTest) 
      rf$test$mse <- rf$test$rsq <- NULL
  }
  rf
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
      print(dim(testing_set))
      print('rf.this_fold votes length...')
      print(length(rf.this_fold$votes))
      print('rf.out votes length...')
      print(length(rf.out$votes))
      #rf.out <- randomForest::combine(rf.out, rf.this_fold) #old code can't handle the 
      rf.out <- my_combine(rf.out, rf.this_fold)
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
      # we need votes for all cells to calculate
      pred <- make.predictions.df(rf.this_fold, testing_set[1:(length(testing_set)-1)], testing_set$Engramcell)
      assess <- assessment( pred ) 
      fold.performance[,ncol(fold.performance) + 1] <- assess
    }
    
  }# end of for loop
  colnames(fold.performance) <- names(folds.obj)
  fold.performance$Mean <- apply(fold.performance, MARGIN=1,  FUN = mean)
  fold.performance$SigDiff <- apply(fold.performance, MARGIN=1,  FUN = sd)
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

# old data
labeled.data.hochgerner2018_genes <-  pavlab.normalize(combined.counts,
                                                       UMI= combined.meta$combined_umi)
genes <- rownames(labeled.data.hochgerner2018_genes)
labeled.data.hochgerner2018_genes <- transpose(labeled.data.hochgerner2018_genes)
colnames(labeled.data.hochgerner2018_genes) <- genes 

# discovery data
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

# RRFclassifier.hoch <- resampled.regularizedRF.crossvalidated
# classifier.hoch.normal <- resampled.randomForest.crossvalidated 
classifier.hochgerner2018_genes <- resampled.randomForest.crossvalidated(data= labeled.data.hochgerner2018_genes,
                                                     under.represented.class = "Inactive",
                                                     over.represented.class = "Active",
                                                     trees.total = 100,
                                                     folds = 10,
                                                     proportion.each.batch=0.8,
                                                     batches.per.fold=20)

# testing resampling method
test <- resampled.randomForest.crossvalidated(data=labeled.data.hochgerner2018_genes,
                                              under.represented.class = "Inactive",
                                              over.represented.class = "Active",
                                              trees.total = 100,
                                              folds = 10,
                                              proportion.each.batch=0.8,
                                              batches.per.fold=5)


# Importance
#classifier.hochgerner2018_genes <- classifier.hoch.normal
importance.df.hoch<- data.frame(gene = as.character( rownames(classifier.hochgerner2018_genes$importance) ),
                                importance_score = as.numeric(classifier.hochgerner2018_genes$importance ) ) %>%
  arrange(desc(importance_score))

head(importance.df.hoch, 1000) # don't forget to run this before gettng  the regularized


write_csv(importance.df.hoch, "JeagerNoReg_Importance_matchedtoHochgerner.csv")
# Importance without regularization
# gene importance_score
# 1   Lingo1        1.9641651
# 2    Inhba        1.9381943
# 3      Arc        1.4624443
# 4    Synpo        1.4025699
# 5   Adgrl3        1.2974050
# 6     Plk2        1.2607054
# 7    Spry2        1.1584135
# 8    Mapk4        1.1555400
# 9    Pcdh8        1.1433587
# 10    Bdnf        1.1265872
# 11     Abr        1.1164898
# 12   Ptgs2        1.1157099
# 13  H2-T23        1.1053527
# 14   Nptx2        1.0128181
# 15   Fmnl1        0.9381681
# 16    Chgb        0.8938230
# 17    Nrn1        0.6753464
# 18 Rtn4rl1        0.6647070
# 19  Epha10        0.6534157
# 20  Shank1        0.6147683

# Sans regularization performance
# > classifier.hochgerner2018_genes$Assessment[,c(11,12)]
# Mean      SigDiff
# F1 Score         0.95706375 0.0199746947
# AUC              0.99964236 0.0005634191
# Precision        1.00000000 0.0000000000
# Recall           0.91836056 0.0364454166
# FPR              0.00000000 0.0000000000
# FNR              0.08163944 0.0364454166
# True Positives  52.90000000 2.2113344387
# False Negatives  4.70000000 2.1000000000
# True Negatives  17.40000000 0.4898979486
# False Positives  0.00000000 0.0000000000

# importance with regularization
# > head(importance.df.hoch, 20)
# gene importance_score
# 1  Lingo1        7.8737083
# 2   Inhba        5.0797660
# 3   Synpo        4.2159676
# 4   Nptx2        4.1663379
# 5   Mapk4        4.0700417
# 6     Arc        3.6015545
# 7  Adgrl3        3.0335451
# 8   Fmnl1        2.5432460
# 9   Ptgs2        2.3351065
# 10   Plk2        2.2557475
# 11   Chgb        1.8660928
# 12 Brinp1        1.4450246
# 13 Epha10        1.2512286
# 14   Bdnf        1.0985766
# 15  Pcdh8        1.0265999
# 16   Acan        0.9969457
# 17    Npy        0.9793275
# 18    Abr        0.9715940
# 19 Hpcal4        0.8472109
# 20  Tmem2        0.8265166

# Regularized Performance
#                        Mean      SigDiff
# F1 Score         0.95489177 0.0294354607
# AUC              0.99964992 0.0007140612
# Precision        1.00000000 0.0000000000
# Recall           0.91515426 0.0524665641
# FPR              0.00000000 0.0000000000
# FNR              0.08484574 0.0524665641
# True Positives  52.70000000 2.8301943396
# False Negatives  4.90000000 3.0479501308
# True Negatives  17.40000000 0.4898979486
# False Positives  0.00000000 0.0000000000



roc.engramcell = roc(labeled.data.hochgerner2018_genes$Engramcell,
                     classifier.hochgerner2018_genes$votes[,2], 
                     plot=TRUE, legacy.axes=TRUE, percent=TRUE,
                     xlab="False Positive Percentage", ylab="True Postive Percentage", 
                     col="firebrick4", lwd=4, print.auc=TRUE)


dev.off()
jpeg("ROCBinarized.jpg", width = 700, height = 700)
plot(roc.engramcell, main = "ROC of RF Classifier")
dev.off()

#heatmap of hochgerner cells gene expression

# this plot needs the cell_df stuff done with the human data

# make the image file
dev.off()
jpeg("Penk_vs_EngramProbability.jpg", width = 500, height = "500")
pheatmap(t(vis.matrix), main = "Hochgerner Cells Activity State",
         cluster_rows = F, cluster_cols=F, 
         annotation_col = cell_df, annotation_col = cell_df,
         show_colnames = F, annotation_names_col = F, annotation_names_row = F)
dev.off()

# Histogram of Engram Probability
# http://www.sthda.com/english/wiki/ggplot2-histogram-plot-quick-start-guide-r-software-and-data-visualization

# make predictions on hochgerner2018 DGCs
bogus.factor <- labeled.data.hochgerner2018_genes$Engramcell
bogus.factor[751:1014] <- labeled.data.hochgerner2018_genes$Engramcell[1:264]
bogus.factor[1:1014] <- levels(bogus.factor)[1]

on.hoch5k <- make.predictions.df(classifier.hochgerner2018_genes,
                                 hoch5k.adultDGCs.lognorm,
                                 meta.data.label.column = bogus.factor)

test <- predict(classifier.hochgerner2018_genes, hoch5k.adultDGCs.lognorm, type = 'prob')
# getting quantile thresholds
ninetyfive= as.numeric( quantile(on.hoch5k$label_pos,0.95) )
ninetysevenpointfive = as.numeric( quantile(on.hoch5k$label_pos,0.975))

df <- on.hoch5k[,c(1,3)] #counts of probability
colnames(df) <- c("label_pos","Predicted")
p <- ggplot(data = df, aes(x=label_pos) )

#giving some weird error on the server
dev.off()
jpeg("Penk_vs_EngramProbabilityCV.jpg", width = 700, height = 700)
p + geom_histogram(color = "darkgreen", fill = "lightgreen") + theme_classic() +
  geom_vline(data=, aes( xintercept=ninetysevenpointfive, color="orange"),
             linetype="dashed") +
  geom_vline(data=, aes( xintercept=ninetyfive, color="red"),
             linetype="dashed") +
  xlab("Probability of being an Engram Cell")+
  ylab("Counts") +
  scale_color_discrete(name = "Thresholds", labels= c("0.975", "0.95") )
dev.off()



## Controls  -- SHUFFLING GENES
# don't forget to use the unshuffled data labels when doing prediction
# otherwise it looks like it's still doing perfect

#run classifer on shuffled hochgerner dataset also train on shuffled combined.counts
# then try shuffling cell ID in combined.counts and training the classifier and running it on hochgerner

## Shuffled Genes in training data


combined.shuffledgenes <- gene.shuffle(  combined.counts ) 
combined.shuffledgenes <- combined.shuffledgenes[rownames(combined.shuffledgenes) %in% shared.genes,]
combined.shuffledgenes <- gene.shuffle(combined.shuffledgenes )
combined.shuffledgenes <- log.norm(combined.shuffledgenes)

# creating our training and validation data, we will take 30% of the
# baseline cells and and equivalent number of the engram cells to maintain balance in the
# testing data

# we do not regenerate our training and test set but
combined.shuffledgenes$Engramcell <- combined.meta$ActivityStatus
combined.shuffledgenes$Engramcell <- as.factor(combined.shuffledgenes$Engramcell)
training_set.shuffled <- combined.shuffledgenes[which( !(combined.meta$idx %in% df.temp$idx) ), ]
validation_set.shuffled <- combined.shuffledgenes[which(combined.meta$idx %in% df.temp$idx), ]


training_set.shuffled$treatment <- combined.meta$treatment[which( !(combined.meta$idx %in% df.temp$idx) ) ]

tic()
shuffledgenes.classifier <- resample.randomForest.crossvalidated( df.in = training_set.shuffled,
                                                                  under.represented.class = "Inactive",
                                                                  over.represented.class = "Active",
                                                                  trees.total = 1000,
                                                                  folds = 10,
                                                                  proportion.each.batch=0.8,
                                                                  batches.per.fold=20)
toc()


importance.df.shuffledgenes <- data.frame(gene = as.character( rownames(shuffledgenes.classifier$importance) ),
                                          importance_score = as.numeric( shuffledgenes.classifier$importance ) ) %>%
  arrange(desc(importance_score))

head(importance.df.shuffledgenes, 10) # few of these genes make any sense, and the performance is low

# > head(importance.df.shuffledgenes, 10)
# gene importance_score
# 1          Dusp9       0.08831070
# 2  6330408A02Rik       0.08071013
# 3          Tm2d1       0.07955446
# 4            Msc       0.07487376
# 5        Ccdc157       0.06941980
# 6          Mylip       0.06885648
# 7           Dtd2       0.06633179
# 8           Tjp1       0.06617151
# 9         Gtpbp8       0.06613198
# 10         Synj2       0.06610811


#roc curve
levels(test.predictions.shuffled$engramobserved) <- c(0,1)
#Plotting ROC...
roc.engramcell <- roc(test.predictions.shuffled$engramobserved, 
                      as.numeric(test.predictions.shuffled$Fos_pos) )

#roc.inactive <- roc(predictions$inactiveobserved, as.numeric(predictions$Fos_neg) )

dev.off()
jpeg("ROC_lognorm.jpg", width = 350, height = "350")
plot(roc.engramcell, col = "red", main = "ROC of RF Classifier")
dev.off()


## Shuffle genes in discovery set
# Normally we see a long tailed distribution of the activity score

hoch5k.shuffledgenes <- apply(hochgerner5k_2018_counts[, hoch5k.GC_Adult.p35.idx], 
                              MARGIN = 1, 
                              FUN = as.integer
)
hoch5k.shuffledgenes <- t(hoch5k.shuffledgenes)
hoch5k.shuffledgenes  <- hoch5k.shuffledgenes[,colnames(hoch5k.shuffledgenes) %in% shared.genes,]
hoch5k.shuffledgenes <- gene.shuffle(hoch5k.shuffledgenes)
hoch5k.shuffledgenes <- log.norm( hoch5k.shuffledgenes )

hoch5k.shuffledgenes$Engramcell <- rep("Inactive", dim(hoch5k.shuffledgenes)[1])

# run with a classifier that works on regular data
on.hoch5kshuffled <- make.predictions.df(classifier.hochgerner2018_genes,
                                         hoch5k.shuffledgenes)

#HISTOGRAM OF SHUFFLED HOCHGERN ENGRAM CELL PROBABILITY
df.shuffled <- on.hoch5kshuffled[,c(1,3)]
p <- ggplot(data = df.shuffled, aes(x=Fos_pos) )

dev.off()
jpeg("Penk_vs_EngramProbability.jpg", width = 350, height = "350")
p + geom_histogram(color = "darkgreen", fill = "lightgreen") + theme_classic() +
  xlab("Probability of being an Engram Cell") +
  ylab("Counts") 
dev.off()


##  Shuffle cell IDs


combined.shuffledcells <- combined.counts 
combined.shuffledcells <- combined.shuffledcells[rownames(combined.shuffledcells) %in% shared.genes, ]
combined.shuffledcells <- log.norm( combined.shuffledcells )
combined.shuffledcells <- combined.shuffledcells[ sample( c(1:nrow(combined.shuffledcells)) ), ]




# creating our training and validation data, we will take 30% of the
# baseline cells and and equivalent number of the engram cells to maintain balance in the
# testing data

# we do not regenerate our training and test set but
combined.shuffledcells$Engramcell <- combined.meta$ActivityStatus
combined.shuffledcells$Engramcell <- as.factor(combined.shuffledcells$Engramcell)
training_set.shuffledcells <- combined.shuffledcells[which( !(combined.meta$idx %in% df.temp$idx) ), ]
validation_set.shuffledcells <- combined.shuffledcells[which(combined.meta$idx %in% df.temp$idx), ]


training_set.shuffledcells$treatment <- combined.meta$treatment[which( !(combined.meta$idx %in% df.temp$idx) ) ]

tic()
shuffledcells.classifier <- resample.randomForest( df.in = training_set.shuffledcells, 
                                                   proportion = 0.8, 
                                                   batches = 20, 
                                                   trees = 1000)
toc()



test.predictions.shuffledcells <- make.predictions.df(shuffledcells.classifier , validation_set.shuffledcells)

rf.performances$shuffled_cells <- assessment(test.predictions.shuffledcells) 

importance.df.shuffledcells <- data.frame(gene = as.character( rownames(shuffledcells.classifier$importance) ),
                                          importance_score = as.numeric( shuffledcells.classifier$importance ) ) %>%
  arrange(desc(importance_score))

head(importance.df.shuffledcells, 10) # few of these genes make any sense, and the performance is low

# > head(importance.df.shuffledcells, 10) # few of these genes make any sense, and the performance is low
# gene importance_score
# 1     Tpx2       0.06932673
# 2    Rsrc2       0.06425252
# 3     Ostc       0.06281410
# 4     Urm1       0.06144273
# 5    Nudt5       0.05987852
# 6  Heatr5a       0.05898927
# 7  Tmem214       0.05861042
# 8  Eif2s3y       0.05778846
# 9     Dtd1       0.05724997
# 10    Pax6       0.05715121


#roc curve
levels(test.predictions.shuffledcells$engramobserved) <- c(0,1)
#Plotting ROC...
roc.engramcell <- roc(test.predictions.shuffledcells$engramobserved, 
                      as.numeric(test.predictions.shuffledcells$Fos_pos) )

#roc.inactive <- roc(predictions$inactiveobserved, as.numeric(predictions$Fos_neg) )

dev.off()
jpeg("ROC_lognorm.jpg", width = 350, height = "350")
plot(roc.engramcell, col = "red", main = "ROC of RF Classifier")
dev.off()




########################################
### APPLYING THIS TO HUMAN DATA
##################################################

library(Seurat)

# 
# Ayhan, F., Kulkarni, A., Berto, S., Sivaprakasam, K., Douglas, C., Lega, B. C.,
# & Konopka, G. (2021). Resolving cellular and molecular diversity along the 
# hippocampal anterior-to-posterior axis in humans. Neuron, 109(13), 2091-2105.
# URL:https://pubmed.ncbi.nlm.nih.gov/34051145/
# GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE160189
#
# git here https://github.com/konopkalab/10x_scRNAseq_HippoAxisSeq
# meta data is available here: https://cells.ucsc.edu/?ds=human-hippo-axis#
# got to Info & Download button in top left of UMAP cell expression plot

##############
# GETTING DGCS OUT OF AYHAN, COMMENTED OUT SO SCRIPT CAN RUN FAST
##########

# ayhan2021_counts <- read.csv("~/test_datasets/Ayhan2021_GSE160189/GSE160189_Hippo_Counts.csv.gz")
# rownames(ayhan2021_counts) <- ayhan2021_counts$gene
# ayhan2021_counts[is.na(ayhan2021_counts)] <- 0
# ayhan2021_counts <- ayhan2021_counts[, c(2:dim(ayhan2021_counts)[2])]
# # 

# ayhan2021_meta <- read.csv("~/test_datasets/Ayhan2021_GSE160189/meta.tsv",
#                            sep = '\t', header = TRUE)


#note that genes are in the first column, might as well just leave it for when we
# need to merge it with jeager anyway.
# subject id's are in the cell_id i.e. "P57_AAAGTAGGTCCAGTAT" "P57_AACCATGGTAAACACA"


# note as per Ayhan et al., 2021 we do not want Den.Gyr3 as it is mostly from a single subject
# ayhanDGC.idx <- ayhan2021_meta$Cell[(ayhan2021_meta$Cluster == "Den.Gyr2")|(ayhan2021_meta$Cluster == "Den.Gyr1")]
# ayhanDGC.idx <- colnames(ayhan2021_counts) %in% ayhanDGC.idx

# ayhanDGC_counts <- ayhan2021_counts[,ayhanDGC.idx]
# 
# write_csv(ayhanDGC_counts, "~/test_datasets/Ayhan2021_GSE160189/ayhanDGC_counts.csv")
# 
# write_csv(ayhan2021_meta[ayhanDGC.idx,],"~/test_datasets/Ayhan2021_GSE160189/ayhanDGC_meta.csv")


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


present.orthologs.idx <- (hg_to_mm$Symbol_hg %in% rownames(ayhanDGC_counts))&(hg_to_mm$Symbol_mm %in% rownames(combined.counts) )
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


# Making IEGs plots for Ayhan

# > which(rownames(ayhanDGC_counts)=="ARC")
# [1] 219
# > 
#   > which(rownames(ayhanDGC_counts)=="FOS")
# [1] 1267
# > 
#   > which(rownames(ayhanDGC_counts)=="INHBA")
# [1] 1800
# > 
#   > which(rownames(ayhanDGC_counts)=="NPTX2")
# [1] 4789

df <- t(ayhanDGC_counts[c(219, 1267, 1800, 4789), ])
df <- apply(df, 2,as.numeric)
df <- data.frame(df)

p.fos <- ggplot(data = df, aes(x=FOS) )
p.arc <- ggplot(data = df, aes(x=ARC) )
p.inhba <- ggplot(data = df, aes(x=INHBA) )
p.nptx2 <- ggplot(data = df, aes(x=NPTX2) )

dev.off()
jpeg("ayhan2021DGCs_foscounts.jpg", width = 700, height = 700)
p.fos + geom_histogram(color = "darkgreen", fill = "lightgreen") + theme_classic() +
  xlab("Reads") + ylab ("Number of Cells") + ggtitle("FOS") +
  theme(axis.text.x=element_text(size=15, face = "bold")) + 
  theme(axis.text.y=element_text(size=15, face = "bold")) +
  theme(axis.title.x=element_text(size=20, face = "bold")) + 
  theme(axis.title.y=element_text(size=20, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.5, size=22, face = "bold"))
dev.off()


dev.off()
jpeg("ayhan2021DGCs_ARCcounts.jpeg", width = 700, height = 700)
p.arc + geom_histogram(color = "darkgreen", fill = "lightgreen") + theme_classic() +
  xlab("Reads") + ylab ("Number of Cells") + ggtitle("ARC") +
  theme(axis.text.x=element_text(size=15, face = "bold")) + 
  theme(axis.text.y=element_text(size=15, face = "bold")) +
  theme(axis.title.x=element_text(size=20, face = "bold")) + 
  theme(axis.title.y=element_text(size=20, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.5, size=22, face = "bold"))
dev.off()


dev.off()
jpeg("ayhan2021DGCs_Inhbacounts.jpeg", width = 700, height = 700)
p.inhba + geom_histogram(color = "darkgreen", fill = "lightgreen") + theme_classic() +
  xlab("Reads") + ylab ("Number of Cells") + ggtitle("INHBA") +
  theme(axis.text.x=element_text(size=15, face = "bold")) + 
  theme(axis.text.y=element_text(size=15, face = "bold")) +
  theme(axis.title.x=element_text(size=20, face = "bold")) + 
  theme(axis.title.y=element_text(size=20, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.5, size=22, face = "bold"))
dev.off()


dev.off()
jpeg("ayhan2021DGCs_Nptx2counts.jpeg", width = 700, height = 700)
p.nptx2 + geom_histogram(color = "darkgreen", fill = "lightgreen") + theme_classic() +
  xlab("Reads") + ylab ("Number of Cells") + ggtitle("NPTX2") +
  theme(axis.text.x=element_text(size=15, face = "bold")) + 
  theme(axis.text.y=element_text(size=15, face = "bold")) +
  theme(axis.title.x=element_text(size=20, face = "bold")) + 
  theme(axis.title.y=element_text(size=20, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.5, size=22, face = "bold"))
dev.off()

# trying with -human-mouse TF one to one orthologs


mouse_tfs <- read.table("/home/acampbell/PavLabEngrams/EngramCellClassifier/mouse_tfs.tsv",
                        sep = "\t", header = TRUE)

tf_hg_to_mm <- hg_to_mm[hg_to_mm$Symbol_mm %in% mouse_tfs$Symbol,]

tf.present.orthologs.idx <- (tf_hg_to_mm$Symbol_hg %in% rownames(ayhanDGC_counts))&(tf_hg_to_mm$Symbol_mm %in% rownames(combined.counts) )
tf_hg_to_mm <- tf_hg_to_mm[tf.present.orthologs.idx,] # filter for matches, ~ 14403 present between Ayhan and Jeager/Lacar

# filter using left join and swithc out mm for hg in jeager/lacar data
combined.normed <- seurat.normalize(combined.counts)
combined.normed.tf_hg <- combined.normed
combined.normed.tf_hg$Symbol_mm <- rownames(combined.normed.tf_hg)

# filtering for orthologs present in both data sets
combined.normed.tf_hg <- left_join(x = tf_hg_to_mm, y = combined.normed.tf_hg, by = "Symbol_mm" )
rownames(combined.normed.tf_hg) <- combined.normed.tf_hg$Symbol_hg
combined.normed.tf_hg <- combined.normed.tf_hg[,c(4:(dim(combined.normed.tf_hg)[2]) )]


ayhanDGC_normed.tf_hg <- seurat.normalize(ayhanDGC_counts)
ayhanDGC_normed.tf_hg$Symbol_hg <- rownames(ayhanDGC_normed.tf_hg)
ayhanDGC_normed.tf_hg <- left_join(x = tf_hg_to_mm, y = ayhanDGC_normed.tf_hg, by = "Symbol_hg" )
rownames(ayhanDGC_normed.tf_hg) <- ayhanDGC_normed.tf_hg$Symbol_hg
ayhanDGC_normed.tf_hg <- ayhanDGC_normed.tf_hg[,c(4:dim(ayhanDGC_normed.tf_hg)[2])]

# shared genes
shared.tf_hg_mm_orthologs <- multi.intersect(c(rownames(ayhanDGC_normed.tf_hg),
                                             rownames(combined.normed.tf_hg) ) )






########## this works but it is not reproducible
########## this works but it is not reproducible
########## this works but it is not reproducible

# first we do gene matching between all datasets
combined.counts.hg <- combined.counts
combined.counts.hg$Symbol_mm <- rownames(combined.counts.hg)
combined.counts.hg <- left_join(x = hg_to_mm, y = combined.counts.hg, by = "Symbol_mm" )
rownames(combined.counts.hg) <- combined.counts.hg$Symbol_hg
combined.counts.hg <- combined.counts.hg[,c(4:(dim(combined.counts.hg)[2]) )]


hochgernerDGC_counts.hg <- hochgernerDGC_counts
colnames(hochgernerDGC_counts.hg) <- colnames(hochgerner5k_2018_counts)[hochgerner5k_2018_meta$cluster_name == "Granule-mature"]
hochgernerDGC_counts.hg$Symbol_mm <- rownames(hochgernerDGC_counts.hg)
hochgernerDGC_counts.hg <- left_join(x = hg_to_mm, y = hochgernerDGC_counts.hg, by = "Symbol_mm" )
rownames(hochgernerDGC_counts.hg) <- hochgernerDGC_counts.hg$Symbol_hg
hochgernerDGC_counts.hg <- hochgernerDGC_counts.hg[,c(4:dim(hochgernerDGC_counts.hg)[2])]

ayhanDGC_counts$Symbol_hg <- rownames(ayhanDGC_counts) # for bringing to 
ayhanDGC_counts <- left_join(x = hg_to_mm, y = ayhanDGC_counts, by = "Symbol_hg" )
rownames(ayhanDGC_counts) <- ayhanDGC_counts$Symbol_hg
ayhanDGC_counts <- ayhanDGC_counts[,c(4:dim(ayhanDGC_counts)[2])]

gene.list <- list(rownames(hochgernerDGC_counts.hg),
                  rownames(hochgernerDGC_counts.hg),
                  rownames(ayhanDGC_counts) )
shared.genes <- multi.intersect(gene.list )
rm(gene.list)

# all the datasets share genes now no need for a script that matches genes
ayhanDGC_counts.integrated <- ayhanDGC_counts %>% 
  rownames_to_column() %>%
  filter(rowname %in% shared.genes) %>%
  select( -rowname)


#Feature selection based on shared variability
# in a thousand classifiers tried and failed there is a way to do it with
# hochgerner instead

# species label, we are going to integrate using hochgerner and 
all.cells <- cbind(ayhanDGC_counts, hochgernerDGC_counts.hg)
all.cells <- cbind(all.cells, combined.counts.hg)

# making meta data and label to split data by
species.idx <- rep("human", dim(ayhanDGC_counts)[2]) 
species.idx <- c(species.idx, rep("mouse", dim(combined.counts.hg)[2]) )
species.idx <- c(species.idx, rep("mouse", dim(hochgernerDGC_counts.hg)[2]) )

experiment <- rep("ayhan2021", dim(ayhanDGC_counts)[2])
experiment <- c( experiment, rep("jeager2018", dim(combined.counts.hg)[2]) )
experiment <- c( experiment, rep("hochgerner2018", dim(hochgernerDGC_counts.hg)[2]) )


#activity <- rep("unlabelled", dim(ayhanDGC_counts)[2]) 
#activity <- c(activity, combined.meta$ActivityStatus)



all.cells.meta <- data.frame(experiment)
all.cells.meta$species <- species.idx
#all.cells.meta$activity <- activity

#rownames(all.cells.meta) <- c(colnames(ayhanDGC_counts), colnames(hochgernerDGC_counts.hg) )
rownames(all.cells.meta) <- c(colnames(ayhanDGC_counts), colnames(combined.counts.hg), colnames(hochgernerDGC_counts.hg) )
all.cells <- cbind(ayhanDGC_counts, hochgernerDGC_counts.hg)

# making a seurat object to normalize using their anchors,
# could potentially do clustering later in order to observe if cells are clustering together
integration_obj <- CreateSeuratObject(counts = all.cells, 
                                      min.cells = 3, 
                                      min.features = 100, 
                                      meta.data = all.cells.meta)

# split the dataset into a list of two seurat objects (by species)
DGC.list <- SplitObject(integration_obj, split.by = "species")

# normalize and identify variable features for each dataset independently
# may have to play with number of genes to include, so we will include more
# genes than seurat recommends because we may find them of use
# we will use all the genes 
DGC.list <- lapply(X = DGC.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = dim(hg_to_mm)[1])
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures( object.list = DGC.list, nfeatures = 2000 )

########## this works but it is not reproducible
########## this works but it is not reproducible
########## this works but it is not reproducible

# training classifier on jeager2018, we only have 2000 genes though this contains the most important
# genes that previous classifiers have gripped on to still the low number of genes may decrease
# the classifiers performance
  
# filtering for orthologs present in all data sets
int.feats.df <- data.frame(Symbol_hg = features)
int.feats.df <- left_join(x = int.feats.df, y = hg_to_mm, by = "Symbol_hg" )
  

combinedcounts.integrated <- combined.counts.hg
combinedcounts.integrated$Symbol_hg <- rownames(combinedcounts.integrated)
combinedcounts.integrated <- left_join(x = int.feats.df, y = combinedcounts.integrated, by = "Symbol_hg" )
rownames(combinedcounts.integrated) <- combinedcounts.integrated$Symbol_hg
combinedcounts.integrated <- combinedcounts.integrated[,c(4:(dim(combinedcounts.integrated)[2]) )]


hochgernerDGC_counts.integrated <-  hochgernerDGC_counts.hg
hochgernerDGC_counts.integrated$Symbol_hg <- rownames(hochgernerDGC_counts.integrated)
hochgernerDGC_counts.integrated <- left_join(x = int.feats.df, y = hochgernerDGC_counts.integrated, by = "Symbol_hg" )
rownames(hochgernerDGC_counts.integrated) <- hochgernerDGC_counts.integrated$Symbol_hg
hochgernerDGC_counts.integrated <- hochgernerDGC_counts.integrated[,c(4:(dim(hochgernerDGC_counts.integrated)[2]) )]


ayhan2021.integrated <- ayhanDGC_counts
ayhan2021.integrated$Symbol_hg <- rownames(ayhan2021.integrated) # for bringing to 
ayhan2021.integrated <- left_join(x = int.feats.df, y = ayhan2021.integrated, by = "Symbol_hg" )
rownames(ayhan2021.integrated) <- ayhan2021.integrated$Symbol_hg
ayhan2021.integrated <- ayhan2021.integrated[,c(4:dim(ayhan2021.integrated)[2])]

capture.output(rownames(combined.counts), 'jeager_dgc_gene_names.csv')

capture.output(rownames(ayhanDGC_counts), 'ayhan_dgc_gene_names.csv')

## training classifier for human data
# add

#normalizations

# vst transform from seurat
# https://cran.r-project.org/web/packages/sctransform/sctransform.pdf
# https://github.com/satijalab/sctransform
# Publications on this method: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02584-9
# https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1874-1#Sec8




# top quantile
# has serious issues as the counts are a serious problem 
topquantile.binarize <-function(scRNAseq_data, percentile, reference_scRNAseq_data = 'None'){
  #takes in data that is gene by cell
  # returns dataframe of cell by gene
  # calculates quantile for each gene
  # if a gene in a cell is in the top xth percentile (example above 90th percentiles)
  # that gene is accepted
  # refernce scRNAseq data may be used 
  if(reference_scRNAseq_data!= 'None'){
    baseline_data <- reference_scRNAseq_data
  }else{
    baseline_data <- scRNAseq_data
  }
  
  binarized <- t(+(scRNAseq_data > matrixStats::rowQuantiles(as.matrix(baseline_data), probs=percentile)))
  binarized <- data.frame(binarized)
  
  return(binarized)
}

# ayhan2021.integrated <- topquantile.binarize(ayhan2021.integrated, threshold)
# combinedcounts.integrated <- opquantile.binarize(combinedcounts.integrated)
# combinedcounts.integrated$Engramcell <- as.factor(combined.meta$ActivityStatus)



# Log - Norma Normalisations
combinedcounts.integrated <- log.norm(combinedcounts.integrated)
combinedcounts.integrated$Engramcell <- as.factor(combined.meta$ActivityStatus)

ayhan2021.integrated <- log.norm(ayhan2021.integrated)

# train classifier
# RRFclassifier.integrated <- resampled.regularizedRF.crossvalidated
# classifier.integrated.normal <- resampled.randomForest.crossvalidated 
classifier.integrated  <- resampled.regularizedRF.crossvalidated( 
                                                        data = combinedcounts.integrated,
                                                        under.represented.class = "Inactive",
                                                        over.represented.class = "Active",
                                                        trees.total = 1000,
                                                        folds = 10,
                                                        proportion.each.batch=0.8,
                                                        batches.per.fold=20)

RRFclassifier.integrated <- classifier.integrated
# Importance
importance.integrated <- data.frame(gene = as.character( rownames(classifier.integrated$importance) ),
                               importance_score = as.numeric(classifier.integrated$importance ) ) %>%
  arrange(desc(importance_score))

head(importance.integrated, 20)

# we get an amazingly high AUC as usual but the important genes are very bizzare
roc.engramcell = roc(combinedcounts.integrated$Engramcell,
                     classifier.integrated$votes[,2], 
                     plot=TRUE, legacy.axes=TRUE, percent=TRUE,
                     xlab="False Positive Percentage", ylab="True Postive Percentage", 
                     col="firebrick4", lwd=4, print.auc=TRUE)

#
dev.off()
jpeg("ROC_IntegrationFeaturesTransform.jpg", width = 700, height = 700)
roc(combinedcounts.integrated$Engramcell,
    classifier.integrated$votes[,2], 
    plot=TRUE, legacy.axes=TRUE, percent=TRUE,
    xlab="False Positive Percentage", ylab="True Postive Percentage", 
    col="firebrick4", lwd=4, print.auc=TRUE)
dev.off()

#Applying it to hochgerner


#applying it to Ayhan
# this bogus factor is just here so make.predictions.df function works
bogus.factor <- combinedcounts.integrated$Engramcell
bogus.factor[751:1014] <- combinedcounts.integrated$Engramcell[1:264]
bogus.factor[1:dim(ayhan2021.integrated)[2]] <- levels(bogus.factor)[1]
ayhan2021.integrated$Engramcell <- bogus.factor

ayhan.integrated.predictions <- make.predictions.df(classifier.integrated,
                                            ayhan2021.integrated,
                                            meta.data.label.column = bogus.factor)

test <- predict(classifier.integrated, ayhan2021.integrated, type = 'prob')

### MAKING FIGURES ON HUMAN DATA

# Histogram of Engram Probability
# http://www.sthda.com/english/wiki/ggplot2-histogram-plot-quick-start-guide-r-software-and-data-visualization

ayhan.prob.count <- data.frame(predict(classifier.integrated, ayhan2021.integrated, type = 'prob') )
ayhan.prob.count$predict <- predict(classifier.integrated, ayhan2021.integrated)

#counts of probability

ninetyfive= as.numeric( quantile(ayhan.prob.count$Active,0.95) )
ninetysevenpointfive = as.numeric( quantile(ayhan.prob.count$Active,0.975))
p <- ggplot(data = ayhan.prob.count, aes(x=Active) )

dev.off()
jpeg("HumanDGC_ENGRAM_Prob_distribution_seuratintegration.jpg", width = 350, height = "350")
p + geom_histogram(color = "darkgreen", fill = "lightgreen") + theme_classic() +
  geom_vline(data=, aes( xintercept=ninetysevenpointfive, color="orange"),
             linetype="dashed") +
  geom_vline(data=, aes( xintercept=ninetyfive, color="red"),
             linetype="dashed") +
  xlab("Probability of Human DGC being an Engram Cell (Seurat Integration Features)")+
  ylab("Counts") +
  scale_color_discrete(name = "Thresholds", labels= c("0.975", "0.95") )
dev.off()

mean(apply(ayhan2021_counts, MARGIN = 2, sum))

dev.off()
jpeg("~/PavLabEngrams/HumanDGC_EGRAMProb_distribution.jpg", width = 350, height = "350")
p + geom_histogram(color = "darkgreen", fill = "lightgreen") + theme_classic() +
  xlab("Probability of Human DGC being an Engram Cell (Seurat integration)")+
  ylab("Counts")
dev.off()


# heatmaps of engram cells vs regular cells


#picking cell indeces for the heatmap
#putative.engram.cells <- which(ayhan.integrated.predictions$predict=="Active")
#putative.engram.cells <- sample( putative.engram.cells, 10, replace = FALSE)
putative.engram.cells <- with(ayhan.integrated.predictions,order(-label_pos))[1:20]
rand.cells <- which(ayhan.integrated.predictions$label_pos<0.33)
rand.cells <- sample( rand.cells, 20, replace = FALSE)
cells.idx <- c(putative.engram.cells, rand.cells)


# picking gene's for the heatmap, important genes vs random ones, 10 important 10 random
top.genes <- importance.integrated$gene[1:20]
iegs <- c("FOS",   "ARC",   "INHBA", "NPAS4", "JUN",   
          "FOSB", "BDNF",  "MAPK4", "JUN","SORCS3")

#convert our gene/cells to a matrix

#convert to a matrix for pheatmap
vis.matrix.importance <- as.matrix(ayhan2021.integrated[cells.idx, top.genes])

vis.matrix <- as.matrix(ayhan2021.integrated[cells.idx, iegs])
vis.matrix <- scale(vis.matrix)

# label the cells
cell_df <- data.frame ("Cells" = c(rep("Putative Engram Cell", 20), rep("Random Cell",20))
)
rownames(cell_df) <- rownames(vis.matrix)
cell_df$Cells <- as.factor(cell_df$Cells)

# make the image file
dev.off()
jpeg("Penk_vs_EngramProbability.jpg", width = 500, height = "500")
pheatmap(t(vis.matrix), main = "Human DGC Activity State",
         cluster_rows = F, cluster_cols=F, 
         annotation_col = cell_df,
         show_colnames = F, annotation_names_col = F, annotation_names_row = F)
dev.off()


# make the image file
dev.off()
jpeg("Penk_vs_EngramProbability.jpg", width = 500, height = "500")
pheatmap(t(vis.matrix), main = "Human DGC Activity State",,
         cluster_rows = F, cluster_cols=F, annotation_names_col = F,
         annotation_col = cell_df, show_colnames = F)
dev.off()

# still got issues, the heatmat is not that clean
# could try more cells
# https://htmlcolorcodes.com/color-picker/


# plotting heatmap, scaling colors
pheatmap(vis.matrix, cluster_rows = F, cluster_cols=F, annotation_row = cell_df,
         annotation_names_col = F, scale = "column", color = color, 
         annotation_colors = list(Cells = c("Putative Engram Cell" = "#2AFE00", 
                                            "Random Cell" = "#ACACAC")),
         show_rownames = F)

pheatmap(vis.matrix, cluster_rows = F, cluster_cols=F, annotation_row = cell_df,
         annotation_names_col = F, color = color, 
         annotation_colors = list(Cells = c("Putative Engram Cell" = "#2AFE00", 
                                            "Random Cell" = "#ACACAC")),
         show_rownames = F)






################
#  NEGATIVE CONTROLS
#########

# cell id shuffle is failing


#, shuffled labels and shuffled genes
# this is giving issues, shuffling ids

# shuffle cell ids

combinedcounts.integrated.shuffled <- transform(combinedcounts.integrated, Engramcell = sample(Engramcell))
#combinedcounts.integrated.shuffled$Engramcell <- sample( c(1:nrow(combinedcounts.integrated))
# kinda a dumb way to do it, literally could just shuffle the fos status instead 

# cell id and gene shuffle
combined.shuffledgenes <- gene.shuffle(  combinedcounts.integrated[,c(1:2000)] )
combinedcounts.integrated.shuffled$Engramcell <- as.factor(sample(combined.meta$ActivityStatus) )

# train classifier
classifier.integrated.shuffled <- resampled.randomForest.crossvalidated( data = combinedcounts.integrated.shuffled,
                                                                under.represented.class = "Inactive",
                                                                over.represented.class = "Active",
                                                                trees.total = 1000,
                                                                folds = 5,
                                                                proportion.each.batch=0.8,
                                                                batches.per.fold=20)

# Importance
importance.integrated.shuffled <- data.frame(gene = as.character( rownames(classifier.integrated.shuffled$importance) ),
                                    importance_score = as.numeric(classifier.integrated.shuffled$importance ) ) %>%
  arrange(desc(importance_score))

head( importance.integrated.shuffled, 10)
# > head(importance.integrated.shuffled, 10)
# gene importance_score
# 1   SLC30A7        0.3161887
# 2     SEL1L        0.2977922
# 3     MYO5B        0.2950458
# 4     NDRG3        0.2941058
# 5      CHGB        0.2808794
# 6  SLC39A10        0.2765453
# 7      NEFM        0.2714302
# 8     GDAP1        0.2671715
# 9     KCNK1        0.2469788
# 10    RCOR1        0.2457436

# testing roc curve
roc.engramcell.shuffled = roc(combinedcounts.integrated$Engramcell,
                              classifier.integrated.shuffled$votes[,2],
                     plot=TRUE, legacy.axes=TRUE, percent=TRUE,
                     xlab="False Positive Percentage", ylab="True Postive Percentage", 
                     col="firebrick4", lwd=4, print.auc=TRUE)


integrated.shuffled.predictions <- make.predictions.df(classifier.integrated.shuffled,
                                                       ayhan2021.integrated,
                                                       meta.data.label.column = bogus.factor)

prob.test <-extractProb(
  models = list(classifier.integrated.shuffled),
  unkX = ayhan2021.integrated)

human.shuffled.proability <- integrated.shuffled.predictions[,2:3] #counts of probability

ninetyfive.shuff = as.numeric( quantile(integrated.shuffled.predictions$label_pos,0.95) )
ninetysevenpointfive.shuff = as.numeric( quantile(integrated.shuffled.predictions$label_pos,0.975))
p <- ggplot(data = human.shuffled.proability , aes(x=label_pos) )

dev.off()
jpeg("HumanDGC_ENGRAM_shuffled_negativecontrol_probdistribution.jpg", width = 350, height = "350")
p + geom_histogram(color = "darkgreen", fill = "lightgreen") + theme_classic() +
  geom_vline(data=, aes( xintercept=ninetysevenpointfive.shuff, color="orange"),
             linetype="dashed") +
  geom_vline(data=, aes( xintercept=ninetyfive.shuff, color="red"),
             linetype="dashed") +
  xlab("Probability of Human DGC being an Engram Cell (Seurat Integration Features)")+
  ylab("Counts") +
  scale_color_discrete(name = "Thresholds", labels= c("0.975", "0.95") )
dev.off()


##############################
#############################

#### REACTIVATION SIGNATURE, 2 class model

############

# we will do arc+/fos+ vs inactive and arc+/fos- cells
arclabels <- read.csv("Jeager2018_GSE98679/Jaeger2018_meta_arclabels.csv", header = TRUE)










# New method for gene selection:  Minimal Wasserstein between two negative binomial distributions

# explanation fo the wasserstein distance
# https://rpubs.com/FJRubio/NWD

# Relevant describes bounds on negative binomials for acurate Wassersteing dis prediciton:
#   Barbour, A. D., Gan, H. L., & Xia, A. (2015). Stein factors for negative binomial
# approximation in Wasserstein distance. Bernoulli, 21(2), 1002-1013.

# Possible libraries and fucnitons:
#   Wassersteing distance calculations:
#   https://www.rdocumentation.org/packages/transport/versions/0.12-4/topics/wasserstein
#   - this can work with
#   https://search.r-project.org/CRAN/refmans/transport/html/wasserstein1d.html
#   -this second one can work with just two vectors of samples
# 
#   Negatvie binomial calcuation:
#     1) fit negative binomial
#     2) sample from a distribution with the parameters fit above:
#   https://stat.ethz.ch/R-manual/R-devel/library/MASS/html/rnegbin.html
#     3) take wasserstein distance using:  https://search.r-project.org/CRAN/refmans/transport/html/wasserstein1d.html
    
# Alternatively this function may be able to represent samples, it describes them as unweighted masses.

ayhantest <- log.norm(ayhanDGC_counts) # needs to be lognormed

sum(is.na(hochgernerDGC_counts))
hochtest <- log.norm(hochgernerDGC_counts) # neeeds to be lonormed???? check

# Run this on allthe genes in hochgerner and ayhan
#wasserstein1d(a, b, p = 1, wa = NULL, wb = NULL)
# we should do wasserstien on transcripts per million normalizations

# a single cell cross species comaprisson tool with references that are useful
# https://www.sciencedirect.com/science/article/pii/S2405471219301991#bib4


# Ding, H., Blair, A., Yang, Y., & Stuart, J. M. (2019). Biological process 
# activity transformation of single cell gene expression for cross-species alignment.
# Nature communications, 10(1), 1-6.
# https://www.nature.com/articles/s41467-019-12924-w

# Seurats addModuleScore?

# Review of interneurons discusses methods of calculating conservation
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7064158/
# this paper has a method for calculating conservation in scRNAseq
# Tosches, M. A., Yamawaki, T. M., Naumann, R. K., Jacobi, A. A., Tushev, G., &
#   Laurent, G. (2018). Evolution of pallium, hippocampus, and cortical cell types
# revealed by single-cell transcriptomics in reptiles. Science, 360(6391), 881-888.
# https://www.science.org/doi/pdf/10.1126/science.aar4237 (turtles to mice)

# another paper 
# Hodge, R. D., Bakken, T. E., Miller, J. A., Smith, K. A., Barkan, E. R., 
# Graybuck, L. T., ... & Lein, E. S. (2018). Conserved cell types with divergent 
# features between human and mouse cortex. BioRxiv, 384826.
# https://www.biorxiv.org/content/10.1101/384826v1.abstract

# Liu, Y., Wang, T., Zhou, B., & Zheng, D. (2021). Robust integration of multiple 
# single-cell RNA sequencing datasets using a single reference space. Nature 
# biotechnology, 39(7), 877-884.

# Hu, Z., Ahmed, A. A., & Yau, C. (2021). CIDER: an interpretable meta-clustering
# framework for single-cell RNA-seq data integration and evaluation. Genome Biology,
# 22(1), 1-21.
# -this paper lists several methods for doing this
# "To test the accuracy of identifying populations, we benchmarked CIDER against 
# other 12 workflows: nine workflows that combined integration approaches and clustering
# (Seurat-CCA [9], fastMNN [6], Scanorama [8], Harmony [11], LIGER [12], Combat [13], 
# Monocle3 [7], Conos [14], and RPCA [10]) and three single-cell clustering approaches
# (Seurat v3-Louvain [5], SC3 [3], and RaceID [4])."

# https://www.biorxiv.org/content/10.1101/164889v1.full


# Lastly there is simply using CCA to find the gene x gene components that explain most of the varience.
# https://cran.r-project.org/web/packages/CCA/CCA.pdf
# https://stats.oarc.ucla.edu/r/dae/canonical-correlation-analysis/
# I could alternatively pass the human and mouse data through the cannonical correlation
# then train the classifier on this data and apply it to the transformed human data.
# use hochgerner as a reference set

##### Testing stuff for downa nd upsampling with cross validation entirely with caret
# to run a randomfroest through caret
feature_1 <- runif(100)
feature_2 <- runif(100)
feature_3 <- c( rnorm(75, mean = 1), rnorm(25, mean = 10) )

labels <- c( rep("A", 75), rep("B", 25) )

test <- data.frame("feature_1" = feature_1, "feature_2" = feature_2, "feature_3" = feature_3)
test$class <- as.factor(labels)

nfolds = 3
folds.test<- createFolds(test$class , k = nfolds)

for (i in c(1:nfolds)){
  
  val.idx <- rownames(test)[ folds.test[[i]] ] 
  train.idx <- which(!(rownames(test) %in% val.idx ) )
  train.data <- test[train.idx,]
  validation.data <- test[val.idx,]
  print(dim(train.data))
  print(dim(validation.data))
  
  
}


# this appears to be the solution you were looking for lol
# https://topepo.github.io/caret/subsampling-for-class-imbalances.html#subsampling-during-resampling
# https://stackoverflow.com/questions/45250252/how-to-downsample-using-r-caret

# Regularized random forest R pacakge
# https://cran.r-project.org/web/packages/RRF/RRF.pdf
# cites this method: https://arxiv.org/pdf/1201.1587.pdf
# installed correctly

# If i am reading the xgboost documentation correctly they are 
# regularizing as well
# https://xgboost.readthedocs.io/en/stable/tutorials/model.html

# https://rdrr.io/cran/xgboost/man/xgb.cv.html
# how to implement cv and subsampling into xgboost
# there is a paramenter called folds
# is is a list where each element in the lsit is a 
# vector of folds so we can downsample  in each fold to match 
# the number of fos- and fos+ cells