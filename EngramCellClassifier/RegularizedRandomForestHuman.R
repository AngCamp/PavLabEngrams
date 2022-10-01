## Trying a regularized random forest  giving a single script here to avoid clutter

library(RRF)

## include the processing that gives us combined counts
# orthologs etce tcet
# include human data

# rewritten RF classifier


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

split <- sample.split(labeled.data.hg$Engramcell, SplitRatio = 0.7)

training_set = subset(labeled.data.hg, split == TRUE)
validation_set = subset(labeled.data.hg, split == FALSE)


test.rrf <- resample.regularizedRF(df.in = training_set,
                       under_represented_class = "Fos-",
                       over_represented_class = "Fos+",
                       proportion= 0.8,
                       batches = 5, 
                       trees = 100)

importance.test.rrf <- data.frame(gene = as.character( rownames(test.rrf$importance) ),
                                             importance_score = as.numeric(test.rrf$importance ) ) %>%
  arrange(desc(importance_score))

predtest <- make.predictions.df(test.rrf , validation_set[1:(length(validation_set)-1)], validation_set$Engramcell)
assesstest <- assessment( predtest )


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


###########  APPLYING THE CLASSIFIERS

test.cvreg <- resampled.regularizedRF.crossvalidated( data = labeled.data.hg,
                                       under.represented.class = "Inactive",
                                       over.represented.class = "Active",
                                       trees.total = 1000,
                                       folds = 5,
                                       proportion.each.batch=0.8,
                                       batches.per.fold=20)


importance.test.cvreg <- data.frame(gene = as.character( rownames(test.cvreg$importance) ),
                                  importance_score = as.numeric(test.cvreg$importance ) ) %>%
  arrange(desc(importance_score))

head(importance.test.cvreg,10)

# > head(importance.test.cvreg,10)
# gene importance_score
# 1   INHBA         9.509352
# 2     ARC         7.031772
# 3  LINGO1         4.919619
# 4    BDNF         4.584055
# 5  ADGRL3         4.208608
# 6   FMNL1         3.888185
# 7   MAPK4         3.785846
# 8   HLA-E         2.691648
# 9  SHANK1         2.570328
# 10  SYNPO         2.556971

 roc(labeled.data.hg$Engramcell,
     test.cvreg$votes[,2], 
     plot=TRUE, legacy.axes=TRUE, percent=TRUE,
     xlab="False Positive Percentage", ylab="True Postive Percentage", 
     col="firebrick4", lwd=4, print.auc=TRUE)

#predicting on human data

regRF.cv.pred <- predict(test.cvreg, 
                         ayhanDGC.lognorm, 
                          type = 'vote', norm.votes = FALSE)

regRF.cv.pred <- data.frame(regRF.cv.pred)
regRF.cv.pred$predict <- predict(test.cvreg, ayhanDGC.lognorm, type = 'response')
regRF.cv.pred$predict <- as.factor(regRF.cv.pred$predict)

### suffling data
labeled.data.hg.shuffled <- labeled.data.hg
labeled.data.hg.shuffled$Engramcell <- as.factor(sample(combined.meta$fos_status) )
 
test.cvreg.shuffle <- resampled.regularizedRF.crossvalidated( data = labeled.data.hg.shuffled,
                                                       under.represented.class = "Fos-",
                                                       over.represented.class = "Fos+",
                                                       trees.total = 1000,
                                                       folds = 5,
                                                       proportion.each.batch=0.8,
                                                       batches.per.fold=20)
 
 
importance.test.cvreg.shuf <- data.frame(gene = as.character( rownames(test.cvreg.shuffle$importance) ),
                                     importance_score = as.numeric(test.cvreg.shuffle$importance ) ) %>%
   arrange(desc(importance_score))
 
head(importance.test.cvreg.shuf, 50)
 
shuffled.pred <- predict(test.cvreg.shuffle, labeled.data.hg, 
                         type = 'vote', norm.votes = FALSE)
shuffled.pred <- data.frame(shuffled.pred)

classifier.integrated.shuffled$votes[,2]


roc(labeled.data.hg$Engramcell,
    shuffled.pred[,2], 
     plot=TRUE, legacy.axes=TRUE, percent=TRUE,
     xlab="False Positive Percentage", ylab="True Postive Percentage", 
     col="firebrick4", lwd=4, print.auc=TRUE)
 
 

### multi class regularized random forest


