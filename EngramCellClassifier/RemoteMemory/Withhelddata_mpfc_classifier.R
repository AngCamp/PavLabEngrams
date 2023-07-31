setwd("~/PavLabEngrams/EngramCellClassifier/RemoteMemory")

# libraries used
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
library(aod) #for logistic regression
#library(lme4) apparently lme4 cannot be loaded with seurat


# we will have to make plots in a seperate R script
# just add relevant figures here and reference the script

#kwon 
kwon2021_counts <- read.table("~/test_datasets/Kwon2021_GSE145970/kwon2021_mpfc_neurons_counts.csv.gz",
                              header = TRUE, sep =",")
rownames(kwon2021_counts) <- kwon2021_counts$GENE

combined.meta <- read.csv("~/test_datasets/Kwon2021_GSE145970/kwon_mpfc_neurons_meta.csv",
                          header = TRUE)

# we now need to filter out cells that were removed during seurat processing
rownames(combined.meta) <- combined.meta$CellID

kwon2021_meta <- combined.meta[combined.meta$CellID %in% colnames(kwon2021_counts),]
kwon2021_counts <- kwon2021_counts[,colnames(kwon2021_counts) %in% kwon2021_meta$CellID]
kwon2021_meta <- combined.meta[colnames(kwon2021_counts),]


# there was a discrepancy between the names in kwon meta dat cell id and the counts
# this commented out code is here to fix it
 fun <- function(x){
    out <- str_c(
      str_sub(x, start = 1L, end = 12L), 
      "_", 
      str_sub(x, start = 13L)
      )
    return(out)}
 
kwon2021_meta$CellID <- lapply(kwon2021_meta$CellID, FUN = fun)

# filter kwon2021_coutns to match cells in the metadata
kwon2021_counts <- kwon2021_counts[,colnames(kwon2021_counts) %in% kwon2021_meta$CellID ]
kwon2021_meta <- kwon2021_meta[kwon2021_meta$CellID %in% colnames(kwon2021_counts),]


#Chen
chen2020_counts <- read.csv('~/PavLabEngrams/EngramCellClassifier/Chen2020_GSE152632/GSE152632_GEO_mberkchen_TRAP2_counts.csv.gz', header = TRUE)

rownames(chen2020_counts) <- chen2020_counts$X
chen2020_counts <- chen2020_counts[,2:3531]

# pull the meta data
chen2020_meta <- read.csv( '~/PavLabEngrams/EngramCellClassifier/Chen2020_GSE152632/SraRunTable.txt', header = TRUE)

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
rownames(chen2020_meta) <- colnames(chen2020_counts)

chen2020_meta <- cbind(chen2020_meta, read.csv("~/PavLabEngrams/EngramCellClassifier/Chen2020_GSE152632/Chen2020_ClusterMarkers.csv") )

# add the umi counts to the meta data
chen2020_meta$umi <- as.numeric(colSums(chen2020_counts[,c(1:(dim(chen2020_counts)[2]) )]))


# FUNCTIONS
# normalization functions


#matching all genes and merging datasets
multi.intersect <- function(x) Reduce(intersect, x) #takes lists of lists, c() will not work


pseudocount_log2p1_transform <- function(x, scale_factor = 10^6, UMI.provided = NULL){
  # Almost as Seurat::NormalizeData but we use log2 rather than natural log
  # from here https://satijalab.org/seurat/reference/normalizedata
  if(is.null(UMI.provided)){
    counts <- sum(x)}else{
      counts <- UMI.provided
    }
  x <- (x)/counts # Feature counts for each cell are divided by the total counts for that cell...
  x <- x*scale_factor # and multiplied by the scale.factor. 
  # the we log2 plus 1 rather than natural log plus 1 seurat uses
  return(log2(x+1))
}

pavlab.normalize <- function(df, UMI = NULL, median.scale=FALSE, scaleby = 10000){
  df.cols <- colnames(df)
  df.rows <- rownames(df)
  if(median.scale){ scaleby = median(UMI)}
  if( is.null(UMI)){
    df <- data.frame(apply(df,  MARGIN = 2, pseudocount_log2p1_transform))
  }else{
    #
    df[] <- Map(pseudocount_log2p1_transform, df, scale_factor = scaleby, UMI.provided = UMI)
    
  }
  colnames(df) <- df.cols
  rownames(df)<- df.rows
  return(df)
}


#normalization functions, log.norm calls logplusone
seurat_log1p_transform <- function(x, scale_factor = 10000, UMI.provided = NULL){
  # as per the LogNormalize option in Seurat::NormalizeData
  # from this URL: https://satijalab.org/seurat/reference/normalizedata
  if(is.null(UMI.provided)){
    counts <- sum(x)}else{
      counts <- UMI.provided
    }
  x <- (x)/counts # Feature counts for each cell are divided by the total counts for that cell... 
  x <- x*scale_factor # and multiplied by the scale.factor. 
  # This is then natural-log transformed using log1p.
  return(log(x+1))
}

seurat.normalize <- function(df, UMI = NULL, median.scale=FALSE, scaleby = 10000){
  df.cols <- colnames(df)
  df.rows <- rownames(df)
  if(median.scale){ scaleby = median(UMI)}
  if( is.null(UMI)){
    df <- data.frame(apply(df,  MARGIN = 2, seurat_log1p_transform))
  }else{
    #
    df[] <- Map(seurat_log1p_transform, df, UMI.provided = UMI, scale_factor = scaleby)
  }
  colnames(df) <- df.cols
  rownames(df)<- df.rows
  return(df)
}

# to write pheatmap images to file as a .png
save_pheatmap <- function(x, filename, width=480, height=960) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  png(filename,width = width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

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

scale.and.rotate<-function(df.in){
  incols <- colnames(df.in)
  inrows <- rownames(df.in)
  df.in <- t(df.in)
  df.in <- scale(df.in)
  df.in <- data.frame(df.in)
  rownames(df.in) <- incols
  colnames(df.in) <- inrows
  return(df.in)
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
                                label = c("Inactive","Active")
){
  #generate predictions for making classifier summary
  predictions <- as.data.frame(predict(classifier.object, test_df[,1:(length(test_df))], type = "prob"))
  predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)] #1:2 for the number of classes
  predictions$observed <- meta.data.label.column #this should be changed if you want to make this functions more modular
  colnames(predictions)[1:2] <- c("label_neg","label_pos")
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
    #roc.list <- c() # roc lsit for ggroc for per fold roc curves
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
      # making roc list
      pred.fold.votes <- predict(rf.out,
                                 testing_set[1:(length(testing_set)-1)], 
                                 type = 'vote')
      head(pred.fold.votes)
      roc.list <- list( roc(testing_set$Engramcell, pred.fold.votes[,2]) )
    }else{
      print(dim(testing_set))
      # making roc list
      pred.fold.votes <- predict(rf.this_fold,
                                 testing_set[1:(length(testing_set)-1)], 
                                 type = 'vote')
      head(pred.fold.votes)
      roc.list <- append(roc.list, list( roc(testing_set$Engramcell, pred.fold.votes[,2]) ) )
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
  rf.out$roc_perfold <- roc.list
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

# we're removing chen cells from this step

# doing it with the pavlab log2p1
chen_normed <- pavlab.normalize(chen2020_counts, UMI = chen2020_meta$umi, median.scale=TRUE)

chen_normed <- data.frame(t(chen_normed))
colnames(chen_normed) <- rownames(chen2020_counts)
rownames(chen_normed) <- colnames(chen2020_counts)

# drop the labeling transcript to prevent teh classifier having an advantage it shouldn't have
chen_normed <- chen_normed[,!(names(chen_normed) %in% "TdTom-transgene")]


###
# First with Recall activate neurons, Fear-Recall condition
##

chen_normed$Engramcell <- chen2020_meta$condition_label == "Fear-Recall" & chen2020_meta$engram_label == "tdT+"

glut_chen_normed <- chen_normed[chen2020_meta$BroadCellTypes =='Excitatory',]
glut_chen_normed$Engramcell <- as.factor(glut_chen_normed$Engramcell)
levels(glut_chen_normed$Engramcell) <- c("Inactive", "Active")

# split the data into training and test
split <- sample.split(glut_chen_normed$Engramcell, SplitRatio = 0.7)

training_set = subset(glut_chen_normed, split == TRUE)
validation_set = subset(glut_chen_normed, split == FALSE)


split.classifier.Excitatory_mpfc <- resample.randomForest( training_set,
                                  under_represented_class = "Active",
                                  over_represented_class = "Inactive",
                                  proportion=0.8,
                                  batches=20, 
                                  trees=500)

saveRDS(split.classifier.Excitatory_mpfc, "FR_rf.RDS")


# Importance
#classifier.hochgerner2018_genes <- classifier.hoch.normal
importance.df <- data.frame(gene = as.character( rownames(split.classifier.Excitatory_mpfc$importance) ),
                                importance_score = as.numeric(split.classifier.Excitatory_mpfc$importance ) ) %>%
  arrange(desc(importance_score))

head(importance.df, 20) # don't forget to run this before getting the regularized one as well


write.csv(importance.df, "FR_split_fr_importantgenes.csv")


# getting roc curve on witheld data
validation_set$predicted <- as.factor(predict(split.classifier.Excitatory_mpfc, 
                                           validation_set,
                                           type="response"))

ypred <- predict(split.classifier.Excitatory_mpfc, 
                 validation_set,
                 type = "prob")

yscore <- data.frame(ypred)
rdb <- cbind(validation_set$Engramcell,yscore)
colnames(rdb) = c('Engramcell','ActivityScore')


jpeg('Chen_ROC_curve_plots/fearrecall_engram_roc.jpeg',width = 600, height=600)
roc(rdb$Engramcell,
    rdb$ActivityScore, 
    plot=TRUE, legacy.axes=TRUE, percent=TRUE, 
    plot.window(xlim=c(100, 0),
                ylim=c(0, 100)),
    main = 'Fear Recall Td+ Excitatory Neuron',
    xlab="False Positive Percentage", ylab="True Postive Percentage", 
    col="firebrick4", lwd=1, print.auc=TRUE)
dev.off()



###
# Now with No Recall, Fear-Only group
###



chen_normed$Engramcell <- chen2020_meta$condition_label == "Fear-Only" & chen2020_meta$engram_label == "tdT+"

glut_chen_normed <- chen_normed[chen2020_meta$BroadCellTypes =='Excitatory',]
glut_chen_normed$Engramcell <- as.factor(glut_chen_normed$Engramcell)
levels(glut_chen_normed$Engramcell) <- c("Inactive", "Active")

# split the data into training and test
split <- sample.split(glut_chen_normed$Engramcell, SplitRatio = 0.7)

training_set = subset(glut_chen_normed, split == TRUE)
validation_set = subset(glut_chen_normed, split == FALSE)


split.classifier.Excitatory_mpfc <- resample.randomForest( training_set,
                                                           under_represented_class = "Active",
                                                           over_represented_class = "Inactive",
                                                           proportion=0.8,
                                                           batches=20, 
                                                           trees=500)

saveRDS(split.classifier.Excitatory_mpfc, "NR_rf.RDS")

# Importance
#classifier.hochgerner2018_genes <- classifier.hoch.normal
importance.df <- data.frame(gene = as.character( rownames(split.classifier.Excitatory_mpfc$importance) ),
                            importance_score = as.numeric(split.classifier.Excitatory_mpfc$importance ) ) %>%
  arrange(desc(importance_score))

head(importance.df, 20) # don't forget to run this before getting the regularized one as well


write.csv(importance.df, "NR_split_fr_importantgenes.csv")

# getting roc curve on witheld data
validation_set$predicted <- as.factor(predict(split.classifier.Excitatory_mpfc, 
                                              validation_set,
                                              type="response"))

ypred <- predict(split.classifier.Excitatory_mpfc, 
                 validation_set,
                 type = "prob")

yscore <- data.frame(ypred)
rdb <- cbind(validation_set$Engramcell,yscore)
colnames(rdb) = c('Engramcell','ActivityScore')


jpeg('Chen_ROC_curve_plots/norecall_control_roc.jpeg',width = 600, height=600)
roc(rdb$Engramcell,
    rdb$ActivityScore, 
    plot=TRUE, legacy.axes=TRUE, percent=TRUE, 
    plot.window(xlim=c(100, 0),
                ylim=c(0, 100)),
    main = 'No Recall Td+ Excitatory Neuron',
    xlab="False Positive Percentage", ylab="True Postive Percentage", 
    col="firebrick4", lwd=1, print.auc=TRUE)
dev.off()



###
# Now with No Fear, Context-Only group
###


chen_normed$Engramcell <- chen2020_meta$condition_label == "Context-Only" & chen2020_meta$engram_label == "tdT+"

glut_chen_normed <- chen_normed[chen2020_meta$BroadCellTypes =='Excitatory',]
glut_chen_normed$Engramcell <- as.factor(glut_chen_normed$Engramcell)
levels(glut_chen_normed$Engramcell) <- c("Inactive", "Active")

# split the data into training and test
split <- sample.split(glut_chen_normed$Engramcell, SplitRatio = 0.7)

training_set = subset(glut_chen_normed, split == TRUE)
validation_set = subset(glut_chen_normed, split == FALSE)


split.classifier.Excitatory_mpfc <- resample.randomForest( training_set,
                                                           under_represented_class = "Active",
                                                           over_represented_class = "Inactive",
                                                           proportion=0.8,
                                                           batches=20, 
                                                           trees=500)

saveRDS(split.classifier.Excitatory_mpfc, "NF_rf.RDS")

# Importance
#classifier.hochgerner2018_genes <- classifier.hoch.normal
importance.df <- data.frame(gene = as.character( rownames(split.classifier.Excitatory_mpfc$importance) ),
                            importance_score = as.numeric(split.classifier.Excitatory_mpfc$importance ) ) %>%
  arrange(desc(importance_score))

head(importance.df, 20) # don't forget to run this before getting the regularized one as well


write.csv(importance.df, "NF_split_fr_importantgenes.csv")


# getting roc curve on witheld data
validation_set$predicted <- as.factor(predict(split.classifier.Excitatory_mpfc, 
                                              validation_set,
                                              type="response"))

ypred <- predict(split.classifier.Excitatory_mpfc, 
                 validation_set,
                 type = "prob")

yscore <- data.frame(ypred)
rdb <- cbind(validation_set$Engramcell,yscore)
colnames(rdb) = c('Engramcell','ActivityScore')


jpeg('Chen_ROC_curve_plots/contextonly_control_roc.jpeg',width = 600, height=600)
roc(rdb$Engramcell,
    rdb$ActivityScore, 
    plot=TRUE, legacy.axes=TRUE, percent=TRUE, 
    plot.window(xlim=c(100, 0),
                ylim=c(0, 100)),
    main = 'No Fear Td+ Excitatory Neuron',
    xlab="False Positive Percentage", ylab="True Postive Percentage", 
    col="firebrick4", lwd=1, print.auc=TRUE)
dev.off()



###
# Now with HC control, Homecage group
###


chen_normed$Engramcell <- chen2020_meta$condition_label == "Homecage" & chen2020_meta$engram_label == "tdT+"

glut_chen_normed <- chen_normed[chen2020_meta$BroadCellTypes =='Excitatory',]
glut_chen_normed$Engramcell <- as.factor(glut_chen_normed$Engramcell)
levels(glut_chen_normed$Engramcell) <- c("Inactive", "Active")

# split the data into training and test
split <- sample.split(glut_chen_normed$Engramcell, SplitRatio = 0.7)

training_set = subset(glut_chen_normed, split == TRUE)
validation_set = subset(glut_chen_normed, split == FALSE)


split.classifier.Excitatory_mpfc <- resample.randomForest( training_set,
                                                           under_represented_class = "Active",
                                                           over_represented_class = "Inactive",
                                                           proportion=0.8,
                                                           batches=20, 
                                                           trees=500)

saveRDS(split.classifier.Excitatory_mpfc, "HC_rf.RDS")

# Importance
importance.df <- data.frame(gene = as.character( rownames(split.classifier.Excitatory_mpfc$importance) ),
                            importance_score = as.numeric(split.classifier.Excitatory_mpfc$importance ) ) %>%
  arrange(desc(importance_score))

head(importance.df, 20) # don't forget to run this before getting the regularized one as well


write.csv(importance.df, "HC_split_fr_importantgenes.csv")

# getting roc curve on witheld data
validation_set$predicted <- as.factor(predict(split.classifier.Excitatory_mpfc, 
                                              validation_set,
                                              type="response"))

ypred <- predict(split.classifier.Excitatory_mpfc, 
                 validation_set,
                 type = "prob")

yscore <- data.frame(ypred)
rdb <- cbind(validation_set$Engramcell,yscore)
colnames(rdb) = c('Engramcell','ActivityScore')


jpeg('Chen_ROC_curve_plots/hc_control_roc.jpeg',width = 600, height=600)
roc(rdb$Engramcell,
    rdb$ActivityScore, 
    plot=TRUE, legacy.axes=TRUE, percent=TRUE, 
    plot.window(xlim=c(100, 0),
                ylim=c(0, 100)),
    main = 'Homecage Td+ Excitatory Neuron',
    xlab="False Positive Percentage", ylab="True Postive Percentage", 
    col="firebrick4", lwd=1, print.auc=TRUE)
dev.off()



###
# Now any positively labelled cell
###


chen_normed$Engramcell <- chen2020_meta$engram_label == "tdT+"

glut_chen_normed <- chen_normed[chen2020_meta$BroadCellTypes =='Excitatory',]
glut_chen_normed$Engramcell <- as.factor(glut_chen_normed$Engramcell)
levels(glut_chen_normed$Engramcell) <- c("Inactive", "Active")

# split the data into training and test
split <- sample.split(glut_chen_normed$Engramcell, SplitRatio = 0.7)

training_set = subset(glut_chen_normed, split == TRUE)
validation_set = subset(glut_chen_normed, split == FALSE)


split.classifier.Excitatory_mpfc <- resample.randomForest( training_set,
                                                           under_represented_class = "Active",
                                                           over_represented_class = "Inactive",
                                                           proportion=0.8,
                                                           batches=20, 
                                                           trees=500)

saveRDS(split.classifier.Excitatory_mpfc, "anypositive_rf.RDS")

# Importance
#classifier.hochgerner2018_genes <- classifier.hoch.normal
importance.df <- data.frame(gene = as.character( rownames(split.classifier.Excitatory_mpfc$importance) ),
                            importance_score = as.numeric(split.classifier.Excitatory_mpfc$importance ) ) %>%
  arrange(desc(importance_score))

head(importance.df, 20) # don't forget to run this before getting the regularized one as well


write.csv(importance.df, "anypositive_split_fr_importantgenes.csv")


# getting roc curve on witheld data
validation_set$predicted <- as.factor(predict(split.classifier.Excitatory_mpfc, 
                                              validation_set,
                                              type="response"))

ypred <- predict(split.classifier.Excitatory_mpfc, 
                 validation_set,
                 type = "prob")

yscore <- data.frame(ypred)
rdb <- cbind(validation_set$Engramcell,yscore)
colnames(rdb) = c('Engramcell','ActivityScore')


jpeg('Chen_ROC_curve_plots/anypositive_control_roc.jpeg',width = 600, height=600)
roc(rdb$Engramcell,
    rdb$ActivityScore, 
    plot=TRUE, legacy.axes=TRUE, percent=TRUE, 
    plot.window(xlim=c(100, 0),
                ylim=c(0, 100)),
    main = 'Any Td+ Excitatory Neuron',
    xlab="False Positive Percentage", ylab="True Postive Percentage", 
    col="firebrick4", lwd=1, print.auc=TRUE)
dev.off()


##
# Comparing important genes
##


fr.impgenes.df <- read.csv("FR_split_fr_importantgenes.csv")
nr.impgenes.df <- read.csv("NR_split_fr_importantgenes.csv")
nf.impgenes.df <- read.csv("NF_split_fr_importantgenes.csv")
hc.impgenes.df <- read.csv("HC_split_fr_importantgenes.csv")
anypos.impgenes.df <- read.csv("anypositive_split_fr_importantgenes.csv")

n <- 1000

multi.intersect(list(fr.impgenes.df$gene[1:n],
                     nr.impgenes.df$gene[1:n],
                     nf.impgenes.df$gene[1:n],
                     hc.impgenes.df$gene[1:n],
                     anypos.impgenes.df$gene[1:n]))

multi.intersect(list(fr.impgenes.df$gene[1:n],
                     nr.impgenes.df$gene[1:n],
                     nf.impgenes.df$gene[1:n],
                     hc.impgenes.df$gene[1:n],
                     anypos.impgenes.df$gene[1:n]))

##
# Percent of cells in each condition
##
total.poscells <- sum(chen2020_meta$engram_label == "tdT+")
sum(chen2020_meta$condition_label == "Fear-Recall" & chen2020_meta$engram_label == "tdT+")/total.poscells
sum(chen2020_meta$condition_label == "Fear-Only" & chen2020_meta$engram_label == "tdT+")/total.poscells
sum(chen2020_meta$condition_label == "Context-Only" & chen2020_meta$engram_label == "tdT+")/total.poscells
sum(chen2020_meta$condition_label == "Homecage" & chen2020_meta$engram_label == "tdT+")/total.poscells
sum(chen2020_meta$engram_label == "tdT+")/total.poscells




glut_chen_normed <- chen_normed[chen2020_meta$BroadCellTypes =='Excitatory',]







#####

# Null Models

#####

# condition shuffle

# here we pull the glutematergic neurons first 
# the shuffling of various label may laed to loss of data if we included the inhibitory
# neurons in the shuffling process
glut.idx <- chen2020_meta$BroadCellTypes =='Excitatory'

glut_chen_normed <- chen_normed[glut.idx ,]

cols.i.want <- c('condition_label','engram_label')
glut_chen_meta <- chen2020_meta[glut.idx,cols.i.want]

# shuffle the meta data, we use sample() to do this
shuffledcondition.glut_chen_meta <- data.frame(engram_label = glut_chen_meta$engram_label,
                                               condition_label = sample(glut_chen_meta$condition_label)
                                              )

 
 
# create labels based on shuffled meta data
# this is done differently here so we don't randomly lose excitatory cells
# when we did the shuffling
glut_chen_normed$Engramcell <- shuffledcondition.glut_chen_meta$condition_label == "Fear-Recall" & shuffledcondition.glut_chen_meta$engram_label == "tdT+"
glut_chen_normed$Engramcell <- as.factor(glut_chen_normed$Engramcell)
levels(glut_chen_normed$Engramcell) <- c("Inactive", "Active")

# split the data into training and test
split <- sample.split(glut_chen_normed$Engramcell, SplitRatio = 0.7)

training_set = subset(glut_chen_normed, split == TRUE)
validation_set = subset(glut_chen_normed, split == FALSE)


conditionshuffle.classifier.Excitatory_mpfc <- resample.randomForest( training_set,
                                                           under_represented_class = "Active",
                                                           over_represented_class = "Inactive",
                                                           proportion=0.8,
                                                           batches=20, 
                                                           trees=500)

saveRDS(conditionshuffle.classifier.Excitatory_mpfc, "conditionshuffle.Excitatory_mpfc_rf.RDS")

# Importance
#classifier.hochgerner2018_genes <- classifier.hoch.normal
importance.df <- data.frame(gene = as.character( rownames(conditionshuffle.classifier.Excitatory_mpfc$importance) ),
                            importance_score = as.numeric(conditionshuffle.classifier.Excitatory_mpfc$importance ) ) %>%
  arrange(desc(importance_score))

head(importance.df, 20) # don't forget to run this before getting the regularized one as well


write.csv(importance.df, "conditionshuffle_split_fr_importantgenes.csv")


# getting roc curve on witheld data
# getting roc curve on witheld data
# we use the object names to find the validation set from that run
val.idx <- !(rownames(chen_normed)  %in% rownames(conditionshuffle.classifier.Excitatory_mpfc$votes))
conditionshuffle.validation_set <-chen_normed[val.idx,]
conditionshuffle.validation_set$predicted <- as.factor(predict(conditionshuffle.classifier.Excitatory_mpfc, 
                                                               conditionshuffle.validation_set,
                                                           type="response"))

ypred <- predict(conditionshuffle.classifier.Excitatory_mpfc, 
                 conditionshuffle.validation_set,
                 type = "prob")

yscore <- data.frame(ypred)
rdb.conditionshuffle <- cbind(conditionshuffle.validation_set$Engramcell,yscore)
colnames(rdb.conditionshuffle) = c('Engramcell','ActivityScore')


conditionshuffle.roc.obj <- roc(rdb.conditionshuffle$Engramcell,
                            rdb.conditionshuffle$ActivityScore)

jpeg('Chen_ROC_curve_plots/conditionshuffle_roc.jpeg',width = 800, height=800)
ggroc(conditionshuffle.roc.obj,  alpha = 1,
      colour = "red", linetype = 'solid',
      size = 3, legacy.axes = TRUE) + 
  theme_classic() + 
  ggtitle("Condition Shuffled") + 
  xlab("FPR") + ylab("TPR") + 
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), linewidth = 2,
               color="darkgrey", linetype="dashed") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text=element_text(size=20, face = "bold", colour="black"),
        axis.title =element_text(size=20, face = "bold", colour="black"),
        panel.border = element_rect(color = "black",
                                    fill = NA,
                                    linewidth = 2)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0)), breaks = c(0,0.25,0.5,0.75,1)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0)),  breaks = c(0.25,0.5,0.75,1))
dev.off()




# label shuffle

# neurons in the shuffling process
glut.idx <- chen2020_meta$BroadCellTypes =='Excitatory'

glut_chen_normed <- chen_normed[glut.idx ,]

cols.i.want <- c('condition_label','engram_label')
glut_chen_meta <- chen2020_meta[glut.idx,cols.i.want]

# shuffle the meta data, we use sample() to do this
shuffledlabel.glut_chen_meta <- data.frame(engram_label = sample(glut_chen_meta$engram_label),
                                               condition_label = glut_chen_meta$condition_label)



# create labels based on shuffled meta data
# this is done differently here so we don't randomly lose excitatory cells
# when we did the shuffling
glut_chen_normed$Engramcell <- shuffledlabel.glut_chen_meta$condition_label == "Fear-Recall" & shuffledlabel.glut_chen_meta$engram_label == "tdT+"
glut_chen_normed$Engramcell <- as.factor(glut_chen_normed$Engramcell)
levels(glut_chen_normed$Engramcell) <- c("Inactive", "Active")

# split the data into training and test
split <- sample.split(glut_chen_normed$Engramcell, SplitRatio = 0.7)

training_set = subset(glut_chen_normed, split == TRUE)
validation_set = subset(glut_chen_normed, split == FALSE)


labelshuffle.classifier.Excitatory_mpfc <- resample.randomForest( training_set,
                                                                      under_represented_class = "Active",
                                                                      over_represented_class = "Inactive",
                                                                      proportion=0.8,
                                                                      batches=20, 
                                                                      trees=500)

saveRDS(labelshuffle.classifier.Excitatory_mpfc, "labelshuffle.classifier.Excitatory_mpfc_rf.RDS")

# Importance
#classifier.hochgerner2018_genes <- classifier.hoch.normal
importance.df <- data.frame(gene = as.character( rownames(labelshuffle.classifier.Excitatory_mpfc$importance) ),
                            importance_score = as.numeric(labelshuffle.classifier.Excitatory_mpfc$importance ) ) %>%
  arrange(desc(importance_score))

head(importance.df, 20) # don't forget to run this before getting the regularized one as well


write.csv(importance.df, "labelshuffle_split_fr_importantgenes.csv")


# getting roc curve on witheld data
# we use the object names to find the validation set from that run
val.idx <- !(rownames(chen_normed)  %in% rownames(labelshuffle.classifier.Excitatory_mpfc$votes))
labelshuffle.validation_set <-chen_normed[val.idx,]
labelshuffle.validation_set$predicted <- as.factor(predict(labelshuffle.classifier.Excitatory_mpfc, 
                                                           labelshuffle.validation_set,
                                                 type="response"))

ypred <- predict(labelshuffle.classifier.Excitatory_mpfc, 
                 labelshuffle.validation_set,
                 type = "prob")

yscore <- data.frame(ypred)
rdb.labelshuffle <- cbind(labelshuffle.validation_set$Engramcell,yscore)
colnames(rdb.labelshuffle) = c('Engramcell','ActivityScore')


labelshuffle.roc.obj <- roc(rdb.labelshuffle$Engramcell,
                  rdb.labelshuffle$ActivityScore)

jpeg('Chen_ROC_curve_plots/labelshuffle_roc.jpeg',width = 800, height=800)
ggroc(labelshuffle.roc.obj,  alpha = 1,
      colour = "red", linetype = 'solid',
      size = 3, legacy.axes = TRUE) + 
  theme_classic() + 
  ggtitle("Activity Label Shuffled") + 
  xlab("FPR") + ylab("TPR") + 
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), linewidth = 2,
               color="darkgrey", linetype="dashed") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text=element_text(size=20, face = "bold", colour="black"),
        axis.title =element_text(size=20, face = "bold", colour="black"),
        panel.border = element_rect(color = "black",
                                    fill = NA,
                                    linewidth = 2)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0)), breaks = c(0,0.25,0.5,0.75,1)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0)),  breaks = c(0.25,0.5,0.75,1))
dev.off()





# Train it on cell type
# when we did the shuffling



# restricted engram label
# try training it on cell type








# Full gene shuffle



###
# Making Prettier plots
###

# use pROC package with the ggroc() function
# ggroc docs here: https://rdrr.io/cran/pROC/man/ggroc.html

# loading the rdds files

# fear recall
fr.classifier <- readRDS("FR_rf.RDS")

#recreating the engram label
chen_normed$Engramcell <- chen2020_meta$condition_label == "Fear-Recall" & chen2020_meta$engram_label == "tdT+"
chen_normed$Engramcell <- as.factor(chen_normed$Engramcell)
levels(chen_normed$Engramcell) <- c("Inactive", "Active")

val.idx <- !(rownames(chen_normed)  %in% rownames(fr.classifier$votes))
fr.validation_set <-chen_normed[val.idx,]
fr.validation_set$predicted <- as.factor(predict(fr.classifier, 
                                              fr.validation_set,
                                              type="response"))

ypred.fr <- predict(fr.classifier, 
                 fr.validation_set,
                 type = "prob")

yscore.fr <- data.frame(ypred.fr)
rdb.fr <- cbind(fr.validation_set$Engramcell,yscore.fr)
colnames(rdb.fr) = c('Engramcell','ActivityScore')


# trying with ggplot

fr.roc.obj <- roc(rdb.fr$Engramcell,
                  rdb.fr$ActivityScore)

jpeg('Chen_ROC_curve_plots/fearrecall_roc.jpeg',width = 800, height=800)
ggroc(fr.roc.obj,  alpha = 1,
      colour = "red", linetype = 'solid',
      size = 3, legacy.axes = TRUE) + 
  theme_classic() + 
  ggtitle("Fear Recall Td+ Excitatory Neurons") + 
  xlab("FPR") + ylab("TPR") + 
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), linewidth = 2,
               color="darkgrey", linetype="dashed") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text=element_text(size=20, face = "bold", colour="black"),
        axis.title =element_text(size=20, face = "bold", colour="black"),
  panel.border = element_rect(color = "black",
                              fill = NA,
                              linewidth = 2)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0)), breaks = c(0,0.25,0.5,0.75,1)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0)),  breaks = c(0.25,0.5,0.75,1))
dev.off()



#plot.margin = margin(t = 0, r = 0, b = 0, l = 0),


# no recall
nr.classifier <- readRDS("NR_rf.RDS")

#recreating the engram label
chen_normed$Engramcell <- chen2020_meta$condition_label == "Fear-Only" & chen2020_meta$engram_label == "tdT+"
chen_normed$Engramcell <- as.factor(chen_normed$Engramcell)
levels(chen_normed$Engramcell) <- c("Inactive", "Active")

val.idx <- !(rownames(chen_normed)  %in% rownames(nr.classifier$votes))
nr.validation_set <-chen_normed[val.idx,]
nr.validation_set$predicted <- as.factor(predict(nr.classifier, 
                                                 nr.validation_set,
                                                 type="response"))

ypred.nr <- predict(nr.classifier, 
                    nr.validation_set,
                    type = "prob")

yscore.nr <- data.frame(ypred.nr)
rdb.nr <- cbind(nr.validation_set$Engramcell,yscore.nr)
colnames(rdb.nr) = c('Engramcell','ActivityScore')

nr.roc.obj <- roc(rdb.nr$Engramcell,
                  rdb.nr$ActivityScore)

jpeg('Chen_ROC_curve_plots/norecall_control_roc.jpeg',width = 800, height=800)
ggroc(nr.roc.obj,  alpha = 1,
      colour = "red", linetype = 'solid',
      size = 3, legacy.axes = TRUE) + 
  theme_classic() + 
  ggtitle("No Recall Td+ Excitatory Neurons") + 
  xlab("FPR") + ylab("TPR") + 
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), linewidth = 2,
               color="darkgrey", linetype="dashed") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text=element_text(size=20, face = "bold", colour="black"),
        axis.title =element_text(size=20, face = "bold", colour="black"),
        panel.border = element_rect(color = "black",
                                    fill = NA,
                                    linewidth = 2),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10) ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0)), breaks = c(0,0.25,0.5,0.75,1)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0)),  breaks = c(0.25,0.5,0.75,1))
dev.off()


# no fear
nf.classifier <- readRDS("NF_rf.RDS")

#recreating the engram label
chen_normed$Engramcell <- chen2020_meta$condition_label == "Context-Only" & chen2020_meta$engram_label == "tdT+"
chen_normed$Engramcell <- as.factor(chen_normed$Engramcell)
levels(chen_normed$Engramcell) <- c("Inactive", "Active")

val.idx <- !(rownames(chen_normed)  %in% rownames(nf.classifier$votes))
nf.validation_set <-chen_normed[val.idx,]
nf.validation_set$predicted <- as.factor(predict(nf.classifier, 
                                                 nf.validation_set,
                                                 type="response"))

ypred.nf <- predict(nf.classifier, 
                    nf.validation_set,
                    type = "prob")

yscore.nf <- data.frame(ypred.nf)
rdb.nf <- cbind(nf.validation_set$Engramcell,yscore.nf)
colnames(rdb.nf) = c('Engramcell','ActivityScore')

nf.roc.obj <- roc(rdb.nf$Engramcell,
                  rdb.nf$ActivityScore)

jpeg('Chen_ROC_curve_plots/contextonly_control_roc.jpeg',width = 600, height=600)
ggroc(nf.roc.obj,  alpha = 1,
      colour = "red", linetype = 'solid',
      size = 3, legacy.axes = TRUE) + 
  theme_classic() + 
  ggtitle("No Fear Td+ Excitatory Neurons") + 
  xlab("FPR") + ylab("TPR") + 
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), linewidth = 2,
               color="darkgrey", linetype="dashed") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text=element_text(size=20, face = "bold", colour="black"),
        axis.title =element_text(size=20, face = "bold", colour="black"),
        panel.border = element_rect(color = "black",
                                    fill = NA,
                                    linewidth = 2),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10) ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0)), breaks = c(0,0.25,0.5,0.75,1)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0)),  breaks = c(0.25,0.5,0.75,1))
dev.off()


# homecage
hc.classifier <- readRDS("HC_rf.RDS")

#recreating the engram label
chen_normed$Engramcell <- chen2020_meta$condition_label == "Homecage" & chen2020_meta$engram_label == "tdT+"
chen_normed$Engramcell <- as.factor(chen_normed$Engramcell)
levels(chen_normed$Engramcell) <- c("Inactive", "Active")

val.idx <- !(rownames(chen_normed)  %in% rownames(hc.classifier$votes))
hc.validation_set <-chen_normed[val.idx,]
hc.validation_set$predicted <- as.factor(predict(hc.classifier, 
                                                 hc.validation_set,
                                                 type="response"))

ypred.hc <- predict(hc.classifier, 
                    hc.validation_set,
                    type = "prob")

yscore.hc <- data.frame(ypred.hc)
rdb.hc <- cbind(hc.validation_set$Engramcell,yscore.hc)
colnames(rdb.hc) = c('Engramcell','ActivityScore')

hc.roc.obj <- roc(rdb.hc$Engramcell,
                  rdb.hc$ActivityScore)

jpeg('Chen_ROC_curve_plots/hc_control_roc.jpeg',width = 600, height=600)
ggroc(hc.roc.obj,  alpha = 1,
      colour = "red", linetype = 'solid',
      size = 3, legacy.axes = TRUE) + 
  theme_classic() + 
  ggtitle("Homecage Td+ Excitatory Neurons") + 
  xlab("FPR") + ylab("TPR") + 
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), linewidth = 2,
               color="darkgrey", linetype="dashed") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text=element_text(size=20, face = "bold", colour="black"),
        axis.title =element_text(size=20, face = "bold", colour="black"),
        panel.border = element_rect(color = "black",
                                    fill = NA,
                                    linewidth = 2),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10) ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0)), breaks = c(0,0.25,0.5,0.75,1)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0)),  breaks = c(0.25,0.5,0.75,1))
dev.off()

# any positive
anypos.classifier <- readRDS("anypositive_rf.RDS")

#recreating the engram label
chen_normed$Engramcell <- chen2020_meta$engram_label == "tdT+"
chen_normed$Engramcell <- as.factor(chen_normed$Engramcell)
levels(chen_normed$Engramcell) <- c("Inactive", "Active")

val.idx <- !(rownames(chen_normed)  %in% rownames(anypos.classifier$votes))
anypos.validation_set <-chen_normed[val.idx,]
anypos.validation_set$predicted <- as.factor(predict(anypos.classifier, 
                                                     anypos.validation_set,
                                                     type="response"))

ypred.anypos <- predict(anypos.classifier, 
                        anypos.validation_set,
                        type = "prob")

yscore.anypos <- data.frame(ypred.anypos)
rdb.anypos <- cbind(anypos.validation_set$Engramcell,yscore.anypos)
colnames(rdb.anypos) = c('Engramcell','ActivityScore')

anypos.roc.obj <- roc(rdb.anypos$Engramcell,
                  rdb.anypos$ActivityScore)

jpeg('Chen_ROC_curve_plots/anypositive_control_roc.jpeg',width = 600, height=600)
ggroc(anypos.roc.obj,  alpha = 1,
      colour = "red", linetype = 'solid',
      size = 3, legacy.axes = TRUE) + 
  theme_classic() + 
  ggtitle("Any Td+ Excitatory Neurons") + 
  xlab("FPR") + ylab("TPR") + 
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), linewidth = 2,
               color="darkgrey", linetype="dashed") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text=element_text(size=20, face = "bold", colour="black"),
        axis.title =element_text(size=20, face = "bold", colour="black"),
        panel.border = element_rect(color = "black",
                                    fill = NA,
                                    linewidth = 2),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10) ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0)), breaks = c(0,0.25,0.5,0.75,1)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0)),  breaks = c(0.25,0.5,0.75,1))
dev.off()

####
# Applying it to Kwon
###



#filter counts by genes shared wiht kwon
genes_shared_with_kwon <- intersect(rownames(chen2020_counts), rownames(kwon2021_counts))
chen2020_counts.kwongenes <- chen2020_counts[rownames(chen2020_counts) %in% genes_shared_with_kwon,]

# normalize and rotate
chen_normed.kwongenes <- pavlab.normalize(chen2020_counts.kwongenes,
                                UMI = chen2020_meta$umi,
                                median.scale=TRUE)

chen_normed.kwongenes <- data.frame(t(chen_normed.kwongenes))
colnames(chen_normed.kwongenes) <- rownames(chen2020_counts.kwongenes)
rownames(chen_normed.kwongenes) <- colnames(chen2020_counts.kwongenes)

chen_normed.kwongenes <- chen_normed.kwongenes[,!(names(chen_normed.kwongenes) %in% "TdTom-transgene")]


chen_normed.kwongenes$Engramcell <- chen2020_meta$condition_label == "Fear-Recall" & chen2020_meta$engram_label == "tdT+"

glut_chen_normed.kwongenes <- chen_normed.kwongenes[chen2020_meta$BroadCellTypes =='Excitatory',]
glut_chen_normed.kwongenes$Engramcell <- as.factor(glut_chen_normed.kwongenes$Engramcell)
levels(glut_chen_normed.kwongenes$Engramcell) <- c("Inactive", "Active")

# split the data into training and test
split <- sample.split(glut_chen_normed.kwongenes$Engramcell, SplitRatio = 0.7)

training_set.kwongenes = subset(glut_chen_normed.kwongenes, split == TRUE)
validation_set.kwongenes = subset(glut_chen_normed.kwongenes, split == FALSE)


split.classifier.Excitatory_mpfc.kwon <- resample.randomForest( training_set.kwongenes,
                                                           under_represented_class = "Active",
                                                           over_represented_class = "Inactive",
                                                           proportion=0.8,
                                                           batches=20, 
                                                           trees=500)



classifier.Excitatory_mpfc <- resampled.randomForest.crossvalidated(
  data= glut_chen_normed.kwongenes,
  under.represented.class = "Active",
  over.represented.class = "Inactive",
  trees.total = 1000,
  folds = 10,
  proportion.each.batch=0.8,
  batches.per.fold=20)

# save the file
saveRDS(classifier.Excitatory_mpfc, "Excitatory_mpfc_remotememory_rfv2.RDS")


classifier.Excitatory_mpfc <- readRDS("Excitatory_mpfc_remotememory_rfv2.RDS")

# filter, normalize, and rotate kwon with the matched genes i.e. kwongenes
#filter....
kwon2021_counts.kwongenes <- kwon2021_counts[rownames(kwon2021_counts) %in% genes_shared_with_kwon,]

# normalize, use UMI from full counts matrix
kwon2021_normed.kwongenes <- pavlab.normalize(kwon2021_counts.kwongenes, 
                                              UMI = colSums(kwon2021_counts),
                                              median.scale = TRUE )                                      )

# rotate t() turns it into a matrix so convert back to dataframe and use counts to
# recover the row and col names which now mcuh be switched row ->col and vice versa
kwon2021_normed.kwongenes <- data.frame(t(kwon2021_normed.kwongenes))
colnames(kwon2021_normed.kwongenes) <- rownames(kwon2021_counts.kwongenes)
rownames(kwon2021_normed.kwongenes) <- colnames(kwon2021_counts.kwongenes)

kwon.controls.normed <- kwon2021_normed.kwongenes[kwon2021_meta$condition=='Control',]

# test <- predict(regular.resampRF.model, hoch5k.adultDGCs.lognorm, type = 'prob')
test.kwon.controls <- as.data.frame(predict(classifier.Excitatory_mpfc,
                                            kwon2021_normed.kwongenes, 
                                                   type = 'prob'))
ninetyfive= as.numeric( quantile(test.kwon.controls$Active,0.95) )
ninetysevenpointfive = as.numeric( quantile(test.kwon.controls$Active,0.975))


# df <- on.hoch5k[,c(2,3)] #counts of probability
# colnames(df) <- c("label_pos","Predicted")
# p <- ggplot(data = df, aes(x=label_pos) )
p <- ggplot(data = test.kwon.controls, aes(x=Active) )

jpeg('Chen_ROC_curve_plots/kwon_engram_score.jpeg',width = 800, height=800)
p + geom_histogram(color = "darkgreen",
                   fill = "lightgreen",
                   binwidth=0.001) + 
  theme_classic() +
  geom_vline(data=, aes( xintercept=ninetysevenpointfive, color="orange"),
             linetype="dashed") +
  geom_vline(data=, aes( xintercept=ninetyfive, color="red"),
             linetype="dashed") +
  geom_vline(data=, aes( xintercept=0.5, color="grey"),
             linetype="dashed") +
  xlab("Probability of being an Engram Cell")+
  ylab("Counts") +
  ggtitle('Kwon et al., (2021) Engram Score Distribution') +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1)) +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text=element_text(size=20, face = "bold", colour="black"),
        axis.title =element_text(size=20, face = "bold", colour="black"),
        panel.border = element_rect(color = "black",
                                    fill = NA,
                                    linewidth = 2),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10) )+
  scale_color_discrete(name = "Thresholds", labels= c("0.975", "0.95") )
dev.off()





