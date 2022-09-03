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

# Loading data


setwd("/home/acampbell/PavLabEngrams/EngramCellClassifier")
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
                                     rownames(lacar2016_snNE_counts)
)#closing list 
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

# we will remove the PTZ treated cells as well before matching its genes
not.ptz <- which(lacar2016_meta$treatment != "PTZ")


#Match the gene sets by adding the gene names as rows, will strip later
DG.idx <- which(jeager2018_meta$predicted_cell_type=="DG")

# we must add 1 to the values of DG.idx and not.ptz to deal with the generow shifting the index 1
combined.counts <- jeager2018_counts[, c(DG.idx,862)] %>% 
  left_join(lacar2016_wc_counts[, c(not.ptz[not.ptz <= 82],83)], by = 'gene' , all.y = TRUE) %>%
  left_join(lacar2016_snHC_counts, by = 'gene', all.y = TRUE) %>%
  left_join(lacar2016_snNE_counts, by = 'gene', all.y = TRUE) #%>% 


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

#get index locations of adult p35 DGCs
hoch5k.GC_Adult.p35.idx <- (hochgerner5k_2018_meta$age.days.=="35") | (hochgerner5k_2018_meta$age.days.=="35*")
hoch5k.GC_Adult.p35.idx <- (hoch5k.GC_Adult.p35.idx) & (hochgerner5k_2018_meta$cluster_name == "Granule-mature")
hoch5k.GC_Adult.p35.idx <- which(hoch5k.GC_Adult.p35.idx)


#### Functions

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
                                label = c("Fos+","Fos-")
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
                       label = c("Fos+","Fos-") 
                       ){
  # returns a vector of assessments to be used to make dataframe summarizing classifiers performance
  # can be used to make df of all calssifiers trained in a single run
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
  rf.out$votes <- predict(object = classifier, newdata = labeled.data, type = 'vote', norm.votes = FALSE)
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




#crossvalidated.resampled.randomforest

#testing cross-validated
classifier <- resampled.randomForest.crossvalidated( data= labeled.data,
                                       under.represented.class = "Fos-",
                                       over.represented.class = "Fos+",
                                       trees.total = 1000,
                                       folds = 10,
                                       proportion.each.batch=0.8,
                                       batches.per.fold=20)

classifier$votes <- predict(object = classifier, newdata = labeled.data, type = 'vote', norm.votes = FALSE)

importance.df.resamptest <- data.frame(gene = as.character( rownames(classifier$importance) ),
                                       importance_score = as.numeric(classifier$importance ) ) %>%
  arrange(desc(importance_score))

head(importance.df.resamptest, 10)


roc.engramcell = roc(labeled.data$Engramcell, classifier$votes[,2], plot=TRUE, legacy.axes=TRUE, percent=TRUE,
    xlab="False Positive Percentage", ylab="True Postive Percentage", 
    col="firebrick4", lwd=4, print.auc=TRUE)


dev.off()
jpeg("ROCBinarized.jpg", width = 700, height = 700)
plot(roc.engramcell, main = "ROC of RF Classifier")
dev.off()


# Exploratory Analysis of Hochgerner
# Making some plots showing long tails of cannonical IEGs
# http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/81-ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page/
library(ggpubr)


# > which(rownames(hochgerner5k_2018_counts)=="Arc")
# [1] 1213
# > which(rownames(hochgerner5k_2018_counts)=="Fos")
# [1] 4716
# > which(rownames(hochgerner5k_2018_counts)=="Inhba")
# [1] 6470
# > which(rownames(hochgerner5k_2018_counts)=="Nptx2")
# [1] 8583

df <- t(hochgerner5k_2018_counts[c(1213, 4716, 6470, 8583), hoch5k.GC_Adult.p35.idx])
df <- apply(df, 2,as.numeric)
df <- data.frame(df)

p.fos <- ggplot(data = df, aes(x=Fos) )
p.arc <- ggplot(data = df, aes(x=Arc) )
p.inhba <- ggplot(data = df, aes(x=Inhba) )
p.nptx2 <- ggplot(data = df, aes(x=Nptx2) )

dev.off()
jpeg("hochgernerDGCs_foscounts.jpg", width = 700, height = 700)
p.fos + geom_histogram(color = "darkgreen", fill = "lightgreen") + theme_classic() +
  xlab("Reads") + ylab ("Number of Cells") + ggtitle("Fos") +
  theme(axis.text.x=element_text(size=15, face = "bold")) + 
  theme(axis.text.y=element_text(size=15, face = "bold")) +
  theme(axis.title.x=element_text(size=20, face = "bold")) + 
  theme(axis.title.y=element_text(size=20, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.5, size=22, face = "bold"))
dev.off()


dev.off()
jpeg("hochgernerDGCs_ARCcounts.jpeg", width = 700, height = 700)
p.arc + geom_histogram(color = "darkgreen", fill = "lightgreen") + theme_classic() +
  xlab("Reads") + ylab ("Number of Cells") + ggtitle("Arc") +
  theme(axis.text.x=element_text(size=15, face = "bold")) + 
  theme(axis.text.y=element_text(size=15, face = "bold")) +
  theme(axis.title.x=element_text(size=20, face = "bold")) + 
  theme(axis.title.y=element_text(size=20, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.5, size=22, face = "bold"))
dev.off()


dev.off()
jpeg("hochgernerDGCs_Inhbacounts.jpeg", width = 700, height = 700)
p.inhba + geom_histogram(color = "darkgreen", fill = "lightgreen") + theme_classic() +
  xlab("Reads") + ylab ("Number of Cells") + ggtitle("Inhba") +
  theme(axis.text.x=element_text(size=15, face = "bold")) + 
  theme(axis.text.y=element_text(size=15, face = "bold")) +
  theme(axis.title.x=element_text(size=20, face = "bold")) + 
  theme(axis.title.y=element_text(size=20, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.5, size=22, face = "bold"))
dev.off()


dev.off()
jpeg("hochgernerDGCs_Nptx2counts.jpeg", width = 700, height = 700)
p.nptx2 + geom_histogram(color = "darkgreen", fill = "lightgreen") + theme_classic() +
  xlab("Reads") + ylab ("Number of Cells") + ggtitle("Nptx2") +
  theme(axis.text.x=element_text(size=15, face = "bold")) + 
  theme(axis.text.y=element_text(size=15, face = "bold")) +
  theme(axis.title.x=element_text(size=20, face = "bold")) + 
  theme(axis.title.y=element_text(size=20, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.5, size=22, face = "bold"))
dev.off()






### Now with Hochgerner genes and running on hochgerner

# Actually running the code with cross validation

labeled.data.hochgerner2018_genes <-  log.norm(combined.counts)

hoch5k.adultDGCs.lognorm <- apply(hochgerner5k_2018_counts[, hoch5k.GC_Adult.p35.idx], 
                                  MARGIN = 1, 
                                  FUN = as.integer)

hoch5k.adultDGCs.lognorm <- log.norm( t(hoch5k.adultDGCs.lognorm) )

#gene matching
shared.genes <- intersect( colnames(hoch5k.adultDGCs.lognorm), colnames(labeled.data.hochgerner2018_genes) )
hoch5k.adultDGCs.lognorm  <- hoch5k.adultDGCs.lognorm[,colnames(hoch5k.adultDGCs.lognorm) %in% shared.genes,]
labeled.data.hochgerner2018_genes <- labeled.data.hochgerner2018_genes[,colnames(labeled.data.hochgerner2018_genes) %in% shared.genes,]
labeled.data.hochgerner2018_genes$Engramcell <- as.factor(combined.meta$fos_status)


classifier.hochgerner2018_genes <- resampled.randomForest.crossvalidated( data= labeled.data.hochgerner2018_genes,
                                                     under.represented.class = "Fos-",
                                                     over.represented.class = "Fos+",
                                                     trees.total = 1000,
                                                     folds = 10,
                                                     proportion.each.batch=0.8,
                                                     batches.per.fold=20)

# Importance
importance.df.hoch<- data.frame(gene = as.character( rownames(classifier.hochgerner2018_genes$importance) ),
                                       importance_score = as.numeric(classifier.hochgerner2018_genes$importance ) ) %>%
  arrange(desc(importance_score))

head(importance.df.hoch, 10)


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
                                 meta.data.label.column =bogus.factor)


# getting quantile thresholds
thresh.df <- as.data.frame( ninetyfive= as.numeric( quantile(on.hoch5k$Fos_pos,0.95) ),
                            ninetysevenpointfive = as.numeric( quantile(on.hoch5k$Fos_pos,0.975))
)

ninetyfive= as.numeric( quantile(on.hoch5k$Fos_pos,0.95) )
ninetysevenpointfive = as.numeric( quantile(on.hoch5k$Fos_pos,0.975))

df <- on.hoch5k[,2:3] #counts of probability
colnames(df) <- c("Fos_pos","Predicted")
p <- ggplot(data = df, aes(x=Fos_pos) )

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


#run classifer on shuffled hochgerner dataset also train on shuffled combined.counts
# then try shuffling cell ID in combined.counts and training the classifier and running it on hochgerner

## Shuffled Genes in training data


combined.shuffledgenes <-gene.shuffle(  combined.counts ) 
combined.shuffledgenes <- combined.shuffledgenes[rownames(combined.shuffledgenes) %in% shared.genes,]
combined.shuffledgenes <- gene.shuffle(combined.shuffledgenes )
combined.shuffledgenes <-  log.norm(combined.shuffledgenes)


# creating our training and validation data, we will take 30% of the
# baseline cells and and equivalent number of the engram cells to maintain balance in the
# testing data

# we do not regenerate our training and test set but
combined.shuffledgenes$Engramcell <- combined.meta$fos_status
combined.shuffledgenes$Engramcell <- as.factor(combined.shuffledgenes$Engramcell)
training_set.shuffled <- combined.shuffledgenes[which( !(combined.meta$idx %in% df.temp$idx) ), ]
validation_set.shuffled <- combined.shuffledgenes[which(combined.meta$idx %in% df.temp$idx), ]


training_set.shuffled$treatment <- combined.meta$treatment[which( !(combined.meta$idx %in% df.temp$idx) ) ]

tic()
shuffledgenes.classifier <- resample.randomForest.crossvalidated( df.in = training_set.shuffled, proportion = 0.8, 
                                                   batches = 20, trees = 1000)
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


## Shuffle genes in validation set

hoch5k.shuffledgenes <- apply(hochgerner5k_2018_counts[, hoch5k.GC_Adult.p35.idx], 
                              MARGIN = 1, 
                              FUN = as.integer
)
hoch5k.shuffledgenes <- t(hoch5k.shuffledgenes)
hoch5k.shuffledgenes  <- hoch5k.shuffledgenes[,colnames(hoch5k.shuffledgenes) %in% shared.genes,]
hoch5k.shuffledgenes <- gene.shuffle(hoch5k.shuffledgenes)
hoch5k.shuffledgenes <- log.norm( hoch5k.shuffledgenes )

hoch5k.shuffledgenes$Engramcell <- rep("Fos-", dim(hoch5k.shuffledgenes)[1])

on.hoch5kshuffled <- make.predictions.df(test.classifier,
                                         hoch5k.shuffledgenes)

#HISTOGRAM OF SHUFFLED HOCHGERN ENGRAM CELL PROBABILITY
df.shuffled <- on.hoch5kshuffled[,2:3]
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
combined.shuffledcells <- log.norm(  combined.shuffledcells )
combined.shuffledcells <- combined.shuffledcells[ sample( c(1:nrow(combined.shuffledcells)) ), ]




# creating our training and validation data, we will take 30% of the
# baseline cells and and equivalent number of the engram cells to maintain balance in the
# testing data

# we do not regenerate our training and test set but
combined.shuffledcells$Engramcell <- combined.meta$fos_status
combined.shuffledcells$Engramcell <- as.factor(combined.shuffledcells$Engramcell)
training_set.shuffledcells <- combined.shuffledcells[which( !(combined.meta$idx %in% df.temp$idx) ), ]
validation_set.shuffledcells <- combined.shuffledcells[which(combined.meta$idx %in% df.temp$idx), ]


training_set.shuffledcells$treatment <- combined.meta$treatment[which( !(combined.meta$idx %in% df.temp$idx) ) ]

tic()
shuffledcells.classifier <- resample.randomForest( df.in = training_set.shuffledcells, proportion = 0.8, 
                                                   batches = 20, trees = 1000)
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


### Aplying it to human data

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

ayhan2021_counts <- read.csv("~/test_datasets/Ayhan2021_GSE160189/GSE160189_Hippo_Counts.csv.gz")
rownames(ayhan2021_counts) <- ayhan2021_counts$gene
ayhan2021_counts[is.na(ayhan2021_counts)] <- 0
ayhan2021_counts <- ayhan2021_counts[, c(2:dim(ayhan2021_counts)[2])]


ayhan2021_meta <- read.csv("~/test_datasets/Ayhan2021_GSE160189/meta.tsv",
                           sep = '\t', header = TRUE)


#note that genes are in the first column, might as well just leave it for when we
# need to merge it with jeager anyway.
# subject id's are in the cell_id i.e. "P57_AAAGTAGGTCCAGTAT" "P57_AACCATGGTAAACACA"

# from alex's ortholog mapping

# note as per Ayhan et al., 2021 we do not want Den.Gyr3 as it is mostly from a single subject
ayhanDGC.idx <- ayhan2021_meta$Cell[(ayhan2021_meta$Cluster == "Den.Gyr2")|(ayhan2021_meta$Cluster == "Den.Gyr1")]
ayhanDGC.idx <- colnames(ayhan2021_counts) %in% ayhanDGC.idx

ayhanDGC_counts <- ayhan2021_counts[,ayhanDGC.idx]

# Orthologue matching
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


# Training model based on mouse expression using orthologos

# Normalize

labeled.data.hg <-  log.norm(combined.counts.hg)
#add label column for training
labeled.data.hg$Engramcell <- as.factor(combined.meta$fos_status)

#human data
ayhanDGCs.lognorm <- apply(ayhanDGC_counts, 
                                  MARGIN = 1, 
                                  FUN = as.integer)

ayhanDGCs.lognorm <- log.norm( t(ayhanDGCs.lognorm) )



# train classifier
classifier.hg <- resampled.randomForest.crossvalidated( data = labeled.data.hg,
                                                        under.represented.class = "Fos-",
                                                        over.represented.class = "Fos+",
                                                        trees.total = 500,
                                                        folds = 5,
                                                        proportion.each.batch=0.8,
                                                        batches.per.fold=20)

# Importance
importance.df.hg <- data.frame(gene = as.character( rownames(classifier.hg$importance) ),
                                importance_score = as.numeric(classifier.hg$importance ) ) %>%
  arrange(desc(importance_score))

head(importance.df.hg, 20)


roc.engramcell = roc(labeled.data.hg$Engramcell,
                     classifier.hg$votes[,2], 
                     plot=TRUE, legacy.axes=TRUE, percent=TRUE,
                     xlab="False Positive Percentage", ylab="True Postive Percentage", 
                     col="firebrick4", lwd=4, print.auc=TRUE)


dev.off()
jpeg("ROCBinarized_hg_to_mm_orthologs.jpg", width = 700, height = 700)
plot(roc.engramcell, main = "ROC of RF Classifier")
dev.off()

# add plot of engram cell activity like in hochgerner data
# make predictions on hochgerner2018 DGCs
bogus.factor <- labeled.data.hg$Engramcell
bogus.factor[751:1014] <- labeled.data.hg$Engramcell[1:264]
bogus.factor[1:dim(ayhanDGC_counts)[2]] <- levels(bogus.factor)[1]

ayhanDGC.predictions <- make.predictions.df(classifier.hg,
                                            ayhanDGCs.lognorm,
                                            meta.data.label.column = bogus.factor)


which(rownames(ayhanDGC_counts)=="NPTX2")

dev.off()
hist( as.numeric(ayhanDGCs.lognorm[4789,]) ) + xlab("log(INHBA)") + ylab("Counts")
dev.off()


# Histogram of Engram Probability
# http://www.sthda.com/english/wiki/ggplot2-histogram-plot-quick-start-guide-r-software-and-data-visualization

df <- ayhanDGC.predictions[,2:3] #counts of probability

ninetyfive= as.numeric( quantile(ayhanDGC.predictions$label_pos,0.95) )
ninetysevenpointfive = as.numeric( quantile(ayhanDGC.predictions$label_pos,0.975))
p <- ggplot(data = df, aes(x=label_pos) )

dev.off()
jpeg("HumanDGC_ENGRAM_Prob_distribution.jpg", width = 350, height = "350")
p + geom_histogram(color = "darkgreen", fill = "lightgreen") + theme_classic() +
  geom_vline(data=, aes( xintercept=ninetysevenpointfive, color="orange"),
             linetype="dashed") +
  geom_vline(data=, aes( xintercept=ninetyfive, color="red"),
             linetype="dashed") +
  xlab("Probability of Human DGC being an Engram Cell")+
  ylab("Counts") +
  scale_color_discrete(name = "Thresholds", labels= c("0.975", "0.95") )
dev.off()

mean(apply(ayhan2021_counts, MARGIN = 2, sum))

dev.off()
jpeg("~/PavLabEngrams/HumanDGC_EGRAMProb_distribution.jpg", width = 350, height = "350")
p + geom_histogram(color = "darkgreen", fill = "lightgreen") + theme_classic() +
  xlab("Probability of Human DGC being an Engram Cell")+
  ylab("Counts")
dev.off()

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


#  Attempt number 2, integration with seurat
# Relvant reading, good review
# 
# Shafer, M. E. (2019). Cross-species analysis of single-cell transcriptomic data.
# Frontiers in cell and developmental biology, 7, 175.
# https://www.frontiersin.org/articles/10.3389/fcell.2019.00175/full

# Affinati, A. H., Sabatini, P. V., True, C., Tomlinson, A. J., Kirigiti, M., 
# Lindsley, S. R., ... & Rupp, A. C. (2021). Cross-species analysis defines the 
# conservation of anatomically segregated VMH neuron populations. Elife, 10, e69065.
# https://elifesciences.org/articles/69065

#We need to load jeager, hochgerner and ayna and check for which genes are present

#then we can put hochgerner and ayhan  together


ayhan2021_counts <- read.csv("~/test_datasets/Ayhan2021_GSE160189/GSE160189_Hippo_Counts.csv.gz")
rownames(ayhan2021_counts) <- ayhan2021_counts$gene
ayhan2021_counts[is.na(ayhan2021_counts)] <- 0
ayhan2021_counts <- ayhan2021_counts[, c(2:dim(ayhan2021_counts)[2])]


ayhan2021_meta <- read.csv("~/test_datasets/Ayhan2021_GSE160189/meta.tsv",
                           sep = '\t', header = TRUE)

# note as per Ayhan et al., 2021 we do not want Den.Gyr3 as it is mostly from a single subject
ayhanDGC.idx <- ayhan2021_meta$Cell[(ayhan2021_meta$Cluster == "Den.Gyr2")|(ayhan2021_meta$Cluster == "Den.Gyr1")]
ayhanDGC.idx <- colnames(ayhan2021_counts) %in% ayhanDGC.idx

ayhanDGC_counts <- ayhan2021_counts[,ayhanDGC.idx]

hg_to_mm <- read.table("hg_mm_1to1_ortho_genes_DIOPT-v8.tsv", sep = '\t', header = TRUE)
# just run left join onto the appropriate column for each dataset

# filter out genes not present in all datasets
present.orthologs.idx <- (hg_to_mm$Symbol_hg %in% rownames(ayhanDGC_counts))&(hg_to_mm$Symbol_mm %in% rownames(combined.counts) )
present.orthologs.idx <- present.orthologs.idx & (hg_to_mm$Symbol_mm %in% rownames(hochgerner5k_2018_counts) )

hg_to_mm <- hg_to_mm[present.orthologs.idx,] # filter for matches, ~ 14403 present between Ayhan and Jeager/Lacar



# filtering for orthologs present in all data sets
ayhanDGC_counts$Symbol_hg <- rownames(ayhanDGC_counts) # for bringing to 
ayhanDGC_counts <- left_join(x = hg_to_mm, y = ayhanDGC_counts, by = "Symbol_hg" )
rownames(ayhanDGC_counts) <- ayhanDGC_counts$Symbol_hg
ayhanDGC_counts <- ayhanDGC_counts[,c(4:dim(ayhanDGC_counts)[2])]


combined.counts.hg <- combined.counts
combined.counts.hg$Symbol_mm <- rownames(combined.counts.hg)
combined.counts.hg <- left_join(x = hg_to_mm, y = combined.counts.hg, by = "Symbol_mm" )
rownames(combined.counts.hg) <- combined.counts.hg$Symbol_hg
combined.counts.hg <- combined.counts.hg[,c(4:(dim(combined.counts.hg)[2]) )]


hochgernerDGC_counts <- hochgerner5k_2018_counts[, hochgerner5k_2018_meta$cluster_name == "Granule-mature"]
hochgernerDGC_counts$Symbol_mm <- rownames(hochgernerDGC_counts)
hochgernerDGC_counts <- left_join(x = hg_to_mm, y = hochgernerDGC_counts, by = "Symbol_mm" )
rownames(hochgernerDGC_counts) <- hochgernerDGC_counts$Symbol_hg
hochgernerDGC_counts <- hochgernerDGC_counts[,c(4:dim(hochgernerDGC_counts)[2])]

# species label, we are going to integrate using hochgerner and 
all.cells <- cbind(ayhanDGC_counts, combined.counts.hg)
all.cells <- cbind(all.cells, hochgernerDGC_counts)

# making meta data and label to split data by
species.idx <- rep("human", dim(ayhanDGC_counts)[2]) 
species.idx <- c(species.idx, rep("mouse", dim(combined.counts.hg)[2]+dim(hochgernerDGC_counts)[2]) )

experiment <- rep("ayhan2021", dim(ayhanDGC_counts)[2])
experiment <- c( experiment, rep("jeager2018", dim(combined.counts.hg)[2]) )
experiment <- c( experiment, rep("hochgerner2018", dim(hochgernerDGC_counts)[2]) )


activity <- rep("unlabelled", dim(ayhanDGC_counts)[2]) 
activity <- c(activity, combined.meta$fos_status)
activity <- c( activity, rep("unlabelled", dim(hochgernerDGC_counts)[2]) )


all.cells.meta <- data.frame(experiment)
all.cells.meta$species <- species.idx
all.cells.meta$activity <- activity

rownames(all.cells.meta) <- colnames(all.cells)

#making a seurat object to normalize using their anchors,
# could potentially do clustering later in order to observe if cells are clustering together
integration_obj <- CreateSeuratObject(counts = all.cells, 
                                      min.cells = 0, 
                                      min.features = 0, 
                                      meta.data = all.cells.meta)

integration_obj@meta.data$species <- all.cells.meta$species

# split the dataset into a list of two seurat objects (by species)
DGC.list <- SplitObject(integration_obj, split.by = "species")

# normalize and identify variable features for each dataset independently
# may have to play with number of genes to include, so we will include more
# genes than seurat recomends becuase we may find them of use
# we will use all the genes 
DGC.list <- lapply(X = DGC.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = dim(hg_to_mm)[1])
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = DGC.list )

DGC.anchors <- FindIntegrationAnchors(object.list = DGC.list, anchor.features = features)

# this command creates an 'integrated' data assay
DGC.combined <- IntegrateData(anchorset = DGC.anchors)

DGC.combined <- ScaleData(DGC.combined , verbose = FALSE)
DGC.combined <- RunPCA(DGC.combined , npcs = 30, verbose = FALSE)
DGC.combined <- RunUMAP(DGC.combined, reduction = "pca", dims = 1:30)
DGC.combined <- FindNeighbors(DGC.combined, reduction = "pca", dims = 1:30)
DGC.combined <- FindClusters(DGC.combined , resolution = 0.5)

# this guives us 15 clusters which is kinda weird

#plotting

#p2 <- DimPlot(DGC.combined, reduction = "umap", label = TRUE, repel = TRUE)


integration_obj@meta.data$activity <- all.cells.meta$activity
integration_obj@meta.data$experiment <- all.cells.meta$experiment

#plotting
p1 <- DimPlot(DGC.combined, reduction = "umap", group.by = "species")
p2 <- DimPlot(DGC.combined, reduction = "umap", group.by = "activity",
              pt.size = 1.5)

dev.off()
jpeg("Human_and_mouse_DGCs_activity.jpg", width = 1400, height = 700)
p1 + p2
dev.off()


#plotting
p3 <- DimPlot(DGC.combined, reduction = "umap", group.by = "experiment", pt.size = 1.5)

dev.off()
jpeg("Human_and_mouse_DGCs_experiemnt.jpg", width = 1400, height = 700)
p1 + p3
dev.off()


#plotting genes
p4 <- FeaturePlot(DGC.combined, features = c("CCK", "PENK"), 
                  pt.size = 1, blend = TRUE)

dev.off()
jpeg("CCKvsPENK_DGCs_human_mouse.jpg", width = 1400, height = 700)
p4
dev.off()

#plotting genes
p5 <- FeaturePlot(DGC.combined, features = c("NPAS4"), 
                  pt.size = 1, split.by ="species")

dev.off()
jpeg("IEGs_human_mouseDGC.jpg", width = 1400, height = 700)
p5
dev.off()


DGC.integrated.list <- SplitObject(DGC.combined, split.by = "species")
humans.seurat <- DGC.integrated.list$human

derp <- as.matrix(GetAssayData(object = humans.seurat, slot = "counts")) # this is not working
# returns 0x0 matrix

# getting a dataframe back out
# https://jspaezp.github.io/sctree/reference/as.data.frame.Seurat.html
# this should help: https://satijalab.org/seurat/articles/essential_commands.html

test <- as.data.frame(DGC.integrated.list$human@assays$integrated)
test <- as.data.frame(DGC.combined)



