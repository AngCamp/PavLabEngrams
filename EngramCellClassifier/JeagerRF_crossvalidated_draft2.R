## this will be the seconds draft of the jeager data which will A run properly and B implement cross validation
# I will also include the controls here, this script can be used to generate the figures I will put into
# the final version in R markdown for the results which will just call the figures and present the code as needed


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

# needs re writing, second attempt in progress below
resample.randomForest <-function( df.in, proportion,
                                  batches, trees){
  #this function resamples from our samples and retrains new models then combines them
  # this is too prevent over fitting on cells
  trees.per.batch <- as.integer(trees/batches)
  n.cells <- trunc( sum(df.in$Engramcell=="Fos-")*proportion)
  batches <- c(1:batches)
  for( batch in batches){
    
    resample.set <- rbind(ssamp(df=df.in[df.in$Engramcell=="Fos-",], n=n.cells,
                                strata=treatment, over=0),
                          ssamp(df=df.in[df.in$Engramcell=="Fos+",], n=n.cells,
                                strata=treatment, over=0)
    )
    resample.set <- resample.set[,2:(length(resample.set))]
    
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



make.predictions.df <- function(classifier.object, test_df){
  #generate predictions for making classifier summary
  predictions <- as.data.frame(predict(classifier.object, test_df[,1:(length(test_df)-1)], type = "prob"))
  predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)] #1:2 for the number of classes
  predictions$observed <- test_df$Engramcell #this should be changed if you want to make this functions more modular
  colnames(predictions)[1:2] <- c("Fos_neg","Fos_pos")
  predictions$engramobserved <- ifelse(predictions$observed=="Fos+", 1, 0)
  predictions$inactiveobserved <- ifelse(predictions$observed=="Fos-", 1, 0)
  return(predictions)
}



# assess a single run of resampled.randomforest
assessment <- function(predictions.df){
  # returns a vector of assessments to be used to make dataframe summarizing classifiers performance
  # can be used to make df of all calssifiers trained in a single run
  TP <- sum((predictions.df$predict == "Fos+")&(predictions.df$observed == "Fos+"))
  TN <- sum((predictions.df$predict == "Fos-")&(predictions.df$observed == "Fos-"))
  FN <- sum((predictions.df$predict == "Fos-")&(predictions.df$observed == "Fos+"))
  FP <- sum((predictions.df$predict == "Fos+")&(predictions.df$observed == "Fos-"))
  
  #precision and recall as well as sumamry stats F1
  precision <- TP/(TP+FP)
  recall <- TP/(TP+FN)
  F1.score = 2 * (precision * recall) / (precision + recall)
  FPR <- FP/(TN+FP)
  FNR <- FN/(TP+FN)
  
  #getting auc
  roc.engramcell <- roc(predictions.df$engramobserved, as.numeric(predictions.df$Fos_pos) )
  AUC <- auc(roc.engramcell)
  
  return( c(F1.score, AUC, precision, recall, FPR, FNR,
            TP, FN, TN, FP) )
}

#crossvalidated.resampled.randomforest

# Actually running the code with cross validation

resamp.combined.lognorm <-  log.norm(combined.counts)

hoch5k.adultDGCs.lognorm <- apply(hochgerner5k_2018_counts[, hoch5k.GC_Adult.p35.idx], 
                                  MARGIN = 1, 
                                  FUN = as.integer
                                  )

hoch5k.adultDGCs.lognorm <- log.norm( t(hoch5k.adultDGCs.lognorm) )

#gene matching
shared.genes <- intersect( colnames(hoch5k.adultDGCs.lognorm), colnames(resamp.combined.lognorm) )
hoch5k.adultDGCs.lognorm  <- hoch5k.adultDGCs.lognorm[,colnames(hoch5k.adultDGCs.lognorm) %in% shared.genes,]
resamp.combined.lognorm <- resamp.combined.lognorm[,colnames(resamp.combined.lognorm) %in% shared.genes,]




resample.randomForest.crossvalidate <-function( data, 
                                                metadata,
                                                label.column,
                                                meta.cols.to.stratify.by,
                                                trees, 
                                                folds, 
                                                labelled.proportion = 0.8, 
                                                batches.per.fold = 20 
                                                ){
  #buidling the resampling
  folds = 3 #for debugging
  trees.per.fold = floor(trees/folds)
  
  label.column = "fos_status"
  meta.cols.to.stratify.by = c("treatment","experiment.label")
  i <- 0
  while (i < (folds-1) ) {
    
    #stratified takes a fraction of the sample,
    # as we procede step wise we need to reajust the proportion to keep the
    #number of samples taken the same, i.e. if we want 5 fold cross validation
    # firs step we take 1/5th leaving 4/5ths so the next iteration we should take
    # 1/4 of the remaining samples as 1/4*4/5=1/5 and so on and so forth,
    # the last step takes place outside of this while loop and will get the remainder if there is some
    over_n = 1/(folds-i)
    # round down otherwise stratified takes 1 more sample than it should
    over_n = floor(over_n*100)/100 
    
    temp.meta <- metadata[!(rownames(metadata) %in% used_samples), ]
    #generates folds in a stratified way, we use this index to pull from data
    temp.meta <- stratified(temp.meta, 
                            c(label.column, meta.cols.to.stratify.by), 
                            size = over_n, 
                            keep.rownames = TRUE)
    
    used_samples <- c(used_samples,temp.meta$rn) #tracking used cells (samples)
    #run resampling here with default parameters on data[,!(colnames(data) %in% temp.meta$rn)]
    # trees here is set to trees per fold
    rf.this_fold <- resample.randomForest(df.in = data[,!(colnames(data) %in% temp.meta$rn)],
                                          proportion = labelled.proportion,
                                          batches = 20, 
                                          trees = trees.per.fold))
    
    #assess add to a dataframe on data[,temp.meta$rn]
    #combine models together
    
    #update i
    i = i+1
  }# end of while loop
  
  # final fold has no need to redraw as we just take whats left
  temp.meta <- combined.meta[!(rownames(combined.meta) %in% used_samples), ]
  
  #run resampling here with default parameters on data[,!(colnames(data) %in% temp.meta$rn)]
  #assess add to a dataframe on data[,temp.meta$rn]
  #combine models together
  
  #make an s3 object with the combined models and assessment dataframe
  #return s3 model + assessment object

}

  #this function resamples from our samples and retrains new models then combines them
  # this is too prevent over fitting on cells
  
  #this needs to be moved to be within a loop occuring over folds
  trees.per.batch <- as.integer(trees/batches)
  n.cells <- trunc( sum(df.in$Engramcell=="Fos-")*proportion)
  batches <- c(1:batches)
  for( batch in batches){
    
    resample.set <- rbind(ssamp(df=df.in[df.in$Engramcell=="Fos-",], n=n.cells,
                                strata=treatment, over=0),
                          ssamp(df=df.in[df.in$Engramcell=="Fos+",], n=n.cells,
                                strata=treatment, over=0)
    )
    resample.set <- resample.set[,2:(length(resample.set))]
    
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



# creating our training and validation data, we will take 30% of the
# baseline cells and and equvilent number of the engram cells to maintain balance in the
# testing data
# NOTE ssamp() unlike sample.split(), will preserve the proportion of treatment groups so
# n = 52 is the number here it will need changing for cross validation
combined.meta$idx <- c(1:750)
df.temp <- rbind(ssamp(df=combined.meta[combined.meta$fos_status=="Fos+",], 
                       n=52, strata=treatment, over=0
),
ssamp(df=combined.meta[combined.meta$fos_status=="Fos-",], 
      n=52, strata=treatment, over=0)
)# end of rbind

resamp.combined.lognorm$Engramcell <- combined.meta$fos_status
resamp.combined.lognorm$Engramcell <- as.factor(resamp.combined.lognorm$Engramcell)



# from here we need to implement cross validation
training_set <- resamp.combined.lognorm[which( !(combined.meta$idx %in% df.temp$idx) ), ]
validation_set <- resamp.combined.lognorm[which(combined.meta$idx %in% df.temp$idx), ]


training_set$treatment <- combined.meta$treatment[which( !(combined.meta$idx %in% df.temp$idx) ) ]




tic()
test.classifier <- resample.randomForest( df.in = training_set, proportion = 0.8, 
                                          batches = 20, trees = 1000)
toc()

test.predictions <- make.predictions.df(test.classifier, validation_set)

rf.performances$resampled<- assessment(test.predictions) 

importance.df.resamptest <- data.frame(gene = as.character( rownames(test.classifier$importance) ),
                                       importance_score = as.numeric( test.classifier$importance ) ) %>%
  arrange(desc(importance_score))

head(importance.df.resamptest, 10)

#roc curve
levels(predictions.Hoch5k.lognorm$engramobserved) <- c(0,1)
#Plotting ROC...
roc.engramcell <- roc(test.predictions$engramobserved, 
                      as.numeric(test.predictions$Fos_pos) )
# there is an error here the predictions.Hoch5k.lognorm$engramobserved is showing only as 1 which cannot be true
# seomthing is wrong with the code don't know where this comes from

#roc.inactive <- roc(predictions$inactiveobserved, as.numeric(predictions$Fos_neg) )

dev.off()
jpeg("ROC_lognorm.jpg", width = 350, height = "350")
plot(roc.engramcell, col = "red", main = "ROC of RF Classifier")
dev.off()





