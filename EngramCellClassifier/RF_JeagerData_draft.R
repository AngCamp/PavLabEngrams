#Designed to run ont he pavlab server

# Load libraries
library(tidyverse)
library(dplyr)
#library(AnnotationDbi)
library(randomForest)
library(rfUtilities) # tutorial here: https://evansmurphy.wixsite.com/evansspatial/random-forest-sdm
#library(data.table)
#library(reshape2)
#library(FactoMineR)
#library(factoextra)
#library(Rtsne)
library(Seurat)
library(stringr)
library(sampler)
#library(patchwork)
#library(metap)
#library(cqn)
#library(GeoTcgaData)
library(caTools)
library(pROC)
#library(CVXR)
library(ggplot2)
library(stats)
library(Dict)
library(pheatmap)


# Loading Lacar et al., (2016)

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






## These functions will save the assesments of the classifiers
#this need changing, specifically test_set should not be a global, it is begging for issues
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


## These functions will save the assesments of the classifiers
#work in progress trying to adjust how thresholding works, which.max is not good for detecting specific threshold
# make.predictions.df <- function(classifier.object, test_df, threshold = "default"){
#   #generate predictions for making classifier summary
#   predictions <- as.data.frame(predict(classifier.object, test_df[,1:(length(test_df)-1)], type = "prob"))
# 
#   
#   if(threshold != "default"){
#     #print("if condtion works")
#     offset <- abs(0.5-threshold)
#     predictions$Fos_pos <- predictions$Fos_pos + offset
#     predictions$Fos_neg <- predictions$Fos_neg - offset
#   }else{
#     predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)] #1:2 for the number of classes
#   }
#   colnames(predictions)[1:2] <- c("Fos_neg","Fos_pos")
#   predictions$observed <- test_df$Engramcell #this should be changed if you want to make this functions more modular
#   
#   predictions$engramobserved <- ifelse(predictions$observed=="Fos+", 1, 0)
#   predictions$inactiveobserved <- ifelse(predictions$observed=="Fos-", 1, 0)
#   return(predictions)
# }


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



### DOWNSAMPLING AND LORNORM FOR PREDICTING HOCHGERNER

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



#downsampling from the count data
combined.meta$idx <- c(1:750)
df.temp <- ssamp(df=combined.meta[combined.meta$fos_status=="Fos+",], n=175, strata=treatment, over=0)
# due to representation we need to set this to n=175 to get 174 cells 

#make the new count matrix and meta data
downsamp.combinedcounts <- cbind(combined.counts[,df.temp$idx],
                                 combined.counts[,combined.meta$fos_status=="Fos-"]
)


downsamp.meta <- rbind(df.temp,
                       combined.meta[combined.meta$fos_status=="Fos-",]
) 

#logtramsform the Jeager data with all genes 
allgenes.df  <- apply(downsamp.combinedcounts, 
                                 MARGIN = 1, 
                                 FUN = logplusone
                      )

allgenes.df   <- scale( t(allgenes.df ) ) #we need it transposed so that the scaling is done per gene not cell
allgenes.df   <- data.frame( t(allgenes.df ) )


allgenes.df$Engramcell <- as.factor( downsamp.meta$fos_status )
  

### Trying with all Jeager Genes before merging with Hochgerner

split <- sample.split(allgenes.df$Engramcell, SplitRatio = 0.7)

training_set.test = subset(allgenes.df, split == TRUE)
validation_set.test = subset(allgenes.df, split == FALSE)

rf.test = randomForest(x = training_set.test[,1:(length(training_set.test)-1)],
                       y = training_set.test$Engramcell,
                       ntree = 500)


rf.test.predictions <- make.predictions.df(rf.test, validation_set.test)

rf.performances$allgenes <- assessment( rf.test.predictions )

rf.performances <- data.frame( allgenes = assessment( rf.test.predictions ) )
rownames(rf.performances) <- c("F1 Score", "AUC", "Precision", "Recall",
                               "FPR", "FNR", "True Positives", "False Negatives", 
                               "True Negatives", "False Positives")

importance.rf.test.df <- data.frame(gene = as.character( rownames(rf.test$importance) ),
                            importance_score = as.numeric( rf.test$importance ) ) %>%
  arrange(desc(importance_score))


# > head(importance.rf.test.df,10)
# gene importance_score
# 1    Inhba        1.7147604
# 2    Nptx2        1.0623350
# 3   Lingo1        1.0594203
# 4    Synpo        0.9782042
# 5      Arc        0.9184942
# 6   H2.T23        0.8468646
# 7   Shank1        0.7821768
# 8     Bdnf        0.7079833
# 9    Spry2        0.6896615
# 10 Sult2b1        0.6793700

# > head(importance.rf.test.df,10)
# gene importance_score
# 1    Inhba        1.7147604
# 2    Nptx2        1.0623350
# 3   Lingo1        1.0594203
# 4    Synpo        0.9782042
# 5      Arc        0.9184942
# 6   H2.T23        0.8468646
# 7   Shank1        0.7821768
# 8     Bdnf        0.7079833
# 9    Spry2        0.6896615
# 10 Sult2b1        0.6793700


#longnorm Hochgerner data
# > dim(hochgerner5k_2018_counts[, hoch5k.GC_Adult.p35.idx])
# [1] 14545  1014
hoch5k.adultDGCs.lognorm <- apply(hochgerner5k_2018_counts[, hoch5k.GC_Adult.p35.idx], 
                                  MARGIN = 1, FUN = as.integer
                                  )

hoch5k.adultDGCs.lognorm <- log.norm( t(hoch5k.adultDGCs.lognorm) )


#lastly we must filter for matching genes this will drop our list to 14296 genes


#doing lognormalization for the training data
downsamp.lognorm.counts <- log.norm(downsamp.combinedcounts)

#restrict to matching genes
shared.genes <- intersect( colnames(hoch5k.adultDGCs.lognorm), colnames(downsamp.lognorm.counts) )
hoch5k.adultDGCs.lognorm  <- hoch5k.adultDGCs.lognorm[,colnames(hoch5k.adultDGCs.lognorm) %in% shared.genes,]
downsamp.lognorm.counts <- downsamp.lognorm.counts[,colnames(downsamp.lognorm.counts) %in% shared.genes,]

hoch5k.adultDGCs.lognorm$Engramcell <- rep("Fos+", dim(hoch5k.adultDGCs.lognorm)[1])
hoch5k.adultDGCs.lognorm[is.na(hoch5k.adultDGCs.lognorm)] <- 0
#hoch5k.adultDGCs.lognorm$data.use <- rep("test",dim(hoch5k.adultDGCs.lognorm)[1]) 

downsamp.lognorm.counts$Engramcell <- as.factor( downsamp.meta$fos_status )
downsamp.lognorm.counts[is.na(downsamp.lognorm.counts)] <- 0
#downsamp.lognorm.counts$data.use <- rep("train/validate",dim(downsamp.lognorm.counts)[1]) 



# all.data <- bind_rows(downsamp.lognorm.counts, hoch5k.adultDGCs.lognorm)
# all.data$Engramcell <- as.factor(all.data$Engramcell)
# all.data$data.use <- as.factor(all.data$data.use)

#downsamp.lognorm.counts <- all.data[all.data$data.use == "train/validate", ]
#downsamp.lognorm.counts <- downsamp.lognorm.counts[,1:(length(downsamp.lognorm.counts)-1)]
#hoch5k.adultDGCs.lognorm <- all.data[all.data$data.use == "test", ]
#hoch5k.adultDGCs.lognorm <- hoch5k.adultDGCs.lognorm[,1:(length(hoch5k.adultDGCs.lognorm)-1)]


#split data
split <- sample.split(downsamp.lognorm.counts$Engramcell, SplitRatio = 0.7)

training_set = subset(downsamp.lognorm.counts, split == TRUE)
validation_set = subset(downsamp.lognorm.counts, split == FALSE)

#make classifier
downsamp.lognorm.classifier.forHoch5k = randomForest(x = training_set[,1:(length(training_set)-1)],
                                                     y = training_set$Engramcell,
                                                     ntree = 500)


downsamp.lognorm.predictions.forHoch5k <- make.predictions.df(downsamp.lognorm.classifier.forHoch5k, 
                                                              validation_set)

rf.performances$lognorm_forHoch5k <- assessment(downsamp.lognorm.predictions.forHoch5k) 

importance.df <- data.frame(gene = as.character( rownames(downsamp.lognorm.classifier.forHoch5k$importance) ),
                            importance_score = as.numeric( downsamp.lognorm.classifier.forHoch5k$importance ) ) %>%
  arrange(desc(importance_score))

# > rf.performances
#                  allgenes    lognorm_forHoch5k
# F1 Score         0.9278351        0.94117647
# AUC              0.9474852        0.98298817
# Precision        1.0000000        0.96000000
# Recall           0.8653846        0.92307692
# FPR              0.0000000        0.03846154
# FNR              0.1346154        0.07692308
# True Positives  45.0000000       48.00000000
# False Negatives  7.0000000        4.00000000
# True Negatives  52.0000000       50.00000000
# False Positives  0.0000000        2.00000000

# > head(importance.df,10)
# gene importance_score
# 1   Inhba        1.9997179
# 2     Arc        1.4385213
# 3   Fmnl1        1.0695626
# 4   Ptgs2        1.0350485
# 5  Lingo1        1.0343618
# 6  Epha10        1.0085155
# 7  H2-T23        0.9617053
# 8   Nptx2        0.8920429
# 9  Sorcs3        0.8731119
# 10 Mrpl13        0.8681066


#classify new data
hoch5k.adultDGCs.lognorm$Engramcell <- as.factor(rep("Fos+", dim(hoch5k.adultDGCs.lognorm)[1]))
hoch5k.adultDGCs.lognorm[is.na(hoch5k.adultDGCs.lognorm)] <- 0


predictions.Hoch5k.lognorm <- make.predictions.df(downsamp.lognorm.classifier.forHoch5k, 
                                                  hoch5k.adultDGCs.lognorm)


table(predictions.Hoch5k.lognorm$predict)
summary(predictions.Hoch5k.lognorm$Fos_pos)


count.df <- as.data.frame(t(hochgerner5k_2018_counts[,hoch5k.GC_Adult.p35.idx]) )

df <- predictions.Hoch5k.lognorm[,2:3]
df$Penk <- hoch5k.adultDGCs.lognorm$Penk
df$Arc <- hoch5k.adultDGCs.lognorm$Arc
df$Inhba <- hoch5k.adultDGCs.lognorm$Inhba
df$Lingo1 <- hoch5k.adultDGCs.lognorm$Lingo1
df$Synpo <- hoch5k.adultDGCs.lognorm$Synpo
df$Mapk4 <- hoch5k.adultDGCs.lognorm$Mapk4
df$penk_count <- as.numeric(count.df$Penk)
df$prob_bin <- as.factor(floor(df$Fos_pos*10)/10)

p <- ggplot(data = df, aes(x=prob_bin, y=penk_count) )


jpeg("Penk_vs_EngramProbability.jpg", width = 350, height = "350")
p + geom_bar(stat="identity")
dev.off()


### RESAMPLING ATTEMPT 2 without using the weird covarience matching of rf.classbalance
#https://github.com/jeffreyevans/rfUtilities/blob/master/R/rf.classBalance.R

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



# creating our training and validation data, we will take 30% of the
# baseline cells and and equvilent number of the engram cells to maintain balance in the
# testing data
combined.meta$idx <- c(1:750)
df.temp <- rbind(ssamp(df=combined.meta[combined.meta$fos_status=="Fos+",], 
                 n=52, strata=treatment, over=0
                 ),
                 ssamp(df=combined.meta[combined.meta$fos_status=="Fos-",], 
                       n=52, strata=treatment, over=0)
                 )# end of rbind

resamp.combined.lognorm$Engramcell <- combined.meta$fos_status
resamp.combined.lognorm$Engramcell <- as.factor(resamp.combined.lognorm$Engramcell)
training_set <- resamp.combined.lognorm[which( !(combined.meta$idx %in% df.temp$idx) ), ]
validation_set <- resamp.combined.lognorm[which(combined.meta$idx %in% df.temp$idx), ]


training_set$treatment <- combined.meta$treatment[which( !(combined.meta$idx %in% df.temp$idx) ) ]
#inputs
#df = training_set

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

#for loop
 #outputs
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



dev.off()
jpeg("ROC_lognorm.jpg", width = 350, height = "350")
hist(as.numeric(on.hoch5k$Fos_prob), main = "Probability of Engram Cell")
dev.off()



# head(importance.df.resamptest,10)
hoch5k.adultDGCs.lognorm$Engramcell <- rep("Fos-", dim(hoch5k.adultDGCs.lognorm)[1])

on.hoch5k <- make.predictions.df(test.classifier,
                                 hoch5k.adultDGCs.lognorm)

table(on.hoch5k$predict)
summary(on.hoch5k$Fos_pos)

dev.off()
jpeg("ROC_lognorm.jpg", width = 350, height = "350")
hist(roc.engramcell, col = "red", main = "ROC of RF Classifier")
dev.off()




count.df <- as.data.frame(t(hochgerner5k_2018_counts[,hoch5k.GC_Adult.p35.idx]) )

df <- on.hoch5k[,2:3]
df$Penk <- hoch5k.adultDGCs.lognorm$Penk
df$Arc <- hoch5k.adultDGCs.lognorm$Arc
df$Inhba <- hoch5k.adultDGCs.lognorm$Inhba
df$Lingo1 <- hoch5k.adultDGCs.lognorm$Lingo1
df$Synpo <- hoch5k.adultDGCs.lognorm$Synpo
df$Mapk4 <- hoch5k.adultDGCs.lognorm$Mapk4
df$penk_count <- as.integer(count.df$Penk)
df$prob_bin <- as.factor(floor(df$Fos_pos*20)/20)

p <- ggplot(data = df[df$penk_count>1,], aes(x=Fos_pos, y=penk_count) )

dev.off()
jpeg("Penk_vs_EngramProbability.jpg", width = 350, height = "350")
p + geom_point()
dev.off()


#the dsitribution of the probability looks very skewed as one would expect
# makes me wonder if we could find a better cut off, 

# > sum(on.hoch5k$Fos_pos>0.325)
# [1] 47

# > sum( on.hoch5k$Fos_pos>quantile(on.hoch5k$Fos_pos,0.95) )
# [1] 51

# > quantile(on.hoch5k$Fos_pos,0.95)
# 95% 
# 0.3355833 

test.predictions.newthresh <- test.predictions
#We need to show the classifier still works well at this threshold though
thresh = as.numeric( quantile(on.hoch5k$Fos_pos,0.99) )
#thresh = 0.5
test.predictions.newthresh$predict[c(test.predictions.newthresh$Fos_pos)>thresh] <- "Fos+"
test.predictions.newthresh$predict[c(test.predictions.newthresh$Fos_pos)<thresh] <- "Fos-"

rf.performances$resampled_new_thresh <- assessment( test.predictions.newthresh ) 
# way too many false positives
rf.performances
sum( as.numeric(on.hoch5k$Fos_pos) > thresh )


newpredict <- c()
for( val in c(on.hoch5k$Fos_prob) ){
  newpredict<-c(newpredict, val>thresh)
}
sum(newpredict)
sum(0.34 < on.hoch5k$Fos_prob )


## Stability of Resampling RF on Hochgerner eta l., (2018)
#install.packages("Dict")
# Dict dcumnetaion on cran https://cran.r-project.org/web/packages/Dict/readme/README.html#:~:text=Overview,(key)%20like%20a%20dictionary.

test.dict <- Dict$new(
  a = c(1,2),
  b = c(3,4,5)
)

key.vec = c("a","b","c")

test.dict[ key.vec[3] ] <- c(6,7,8,9)
test.dict["empty"] <- c()


#set up dicts to be filled in for loop
engram.dict <-Dict$new(
  a = "Start"
)

ninetyfithquantile.dict <-Dict$new(
  a = "Start"
)

ninetysevenpointfivequantile.dict <-Dict$new(
  a = "Start"
)

temp.pred.dict <- Dict$new(
  a = "start"
)

run <- as.character( c(1:10) )

tic()
for(i in c(1:10) ){
  
  #resample training and validation set
  combined.meta$idx <- c(1:750)
  df.temp <- rbind(ssamp(df=combined.meta[combined.meta$fos_status=="Fos+",], 
                         n=52, strata=treatment, over=0
  ),
  ssamp(df=combined.meta[combined.meta$fos_status=="Fos-",], 
        n=52, strata=treatment, over=0)
  )# end of rbind
  
  resamp.combined.lognorm$Engramcell <- combined.meta$fos_status
  resamp.combined.lognorm$Engramcell <- as.factor(resamp.combined.lognorm$Engramcell)
  training_set <- resamp.combined.lognorm[which( !(combined.meta$idx %in% df.temp$idx) ), ]
  validation_set <- resamp.combined.lognorm[which(combined.meta$idx %in% df.temp$idx), ]
  
  # reinstatiate classifier, validate save the results
  training_set$treatment <- combined.meta$treatment[which( !(combined.meta$idx %in% df.temp$idx) ) ]
  
  test.classifier <- resample.randomForest( df.in = training_set, proportion = 0.8, 
                                            batches = 20, trees = 600)
  temp.pred.df <- make.predictions.df(test.classifier, validation_set)
  temp.pred.dict[ run[i] ] <- temp.pred.df
  
  print(assessment(temp.pred.df))
  
  
  # label hoch cells
  on.hoch5k <- make.predictions.df(test.classifier,
                                   hoch5k.adultDGCs.lognorm)
  
  engram.dict[ run[i] ] <- which(on.hoch5k$predict=="Fos+")
  thresh <- as.numeric( quantile(on.hoch5k$Fos_pos,0.975) )
  ninetysevenpointfivequantile.dict[ run[i] ] <- which( as.numeric(on.hoch5k$Fos_pos) > thresh )
  thresh = as.numeric( quantile(on.hoch5k$Fos_pos,0.95) )
  ninetyfithquantile.dict[ run[i] ] <- which( as.numeric(on.hoch5k$Fos_pos) > thresh )
}
toc()

shared.cells <- multi.intersect(list(ninetyfithquantile.dict["1"], ninetyfithquantile.dict["2"],
                                     ninetyfithquantile.dict["3"], ninetyfithquantile.dict["4"],
                                     ninetyfithquantile.dict["5"], ninetyfithquantile.dict["6"],
                                     ninetyfithquantile.dict["7"], ninetyfithquantile.dict["8"],
                                     ninetyfithquantile.dict["9"], ninetyfithquantile.dict["10"])
                                )
# > shared.cells
# [1] 282 344 381 470 541 619 666 705 729 735 740 766 793

shared.cells <- multi.intersect(list(ninetysevenpointfivequantile.dict["1"], ninetysevenpointfivequantile.dict["2"],
                                     ninetysevenpointfivequantile.dict["3"], ninetysevenpointfivequantile.dict["4"],
                                     ninetysevenpointfivequantile.dict["5"], ninetysevenpointfivequantile.dict["6"],
                                     ninetysevenpointfivequantile.dict["7"], ninetysevenpointfivequantile.dict["8"],
                                     ninetysevenpointfivequantile.dict["9"], ninetysevenpointfivequantile.dict["10"])
                                )

# > shared.cells
# [1] 282 344 470 541 619 666 705 729 735 793

#In the engram.dict cell 282 comes up 5 times, cell 793 and cell 619 also
# got labelled with 282 twice and once on their own each on three runs no cells were labelled

### Heatmap visualization

#picking cell indeces for the heatmap
stable.cells <- c(282, 344,470, 541, 619, 666, 705, 729, 735, 793)
rand.cells <- which(!(c(1:1014) %in% stable.cells))
rand.cells <- sample( rand.cells, 10, replace = FALSE)
cells.idx <- c(stable.cells,rand.cells)

# picking gene's for the heatmap, important genes vs random ones, 10 important 10 random
head(importance.df.resamptest, 10)

top.genes <- importance.df.resamptest$gene[1:20]

rand.genes <- sample(importance.df.resamptest$gene[21:length(importance.df.resamptest$gene)],
                     20, replace = FALSE)

genes.idx <- c(top.genes, rand.genes)

#convert our gene/cells to a matrix

#convert to a matrix for pheatmap
vis.matrix <- as.matrix(hoch5k.adultDGCs.lognorm[cells.idx, genes.idx])
vis.matrix_scaled <- scale(vis.matrix)

#label the genes
gene_df <- data.frame ("Genes" = c(rep("Important Gene", 20), rep("Random Gene",20))
                       )
rownames(gene_df) <- genes.idx
# label the cells
cell_df <- data.frame ("Cells" = c(rep("Putative Engram Cell", 10), rep("Random Cell",10))
)
rownames(cell_df) <- cells.idx

# make the image file
dev.off()
jpeg("Penk_vs_EngramProbability.jpg", width = 500, height = "500")
pheatmap(t(vis.matrix), main = "Hochgerner Cells Activity State",
         cluster_rows = F, cluster_cols=F, 
         annotation_col = cell_df, annotation_row = gene_df,
         show_colnames = F, annotation_names_col = F, annotation_names_row = F)
dev.off()





for (i in run){
  print(engram.dict[i])
  print(ninetysevenpointfivequantile.dict[ run ])
  print(ninetyfithquantile.dict[ run ])
}


test.predictions <- make.predictions.df(test.classifier, validation_set)

rf.performances$resampled<- assessment(test.predictions) 

importance.df.resamptest <- data.frame(gene = as.character( rownames(test.classifier$importance) ),
                                       importance_score = as.numeric( test.classifier$importance ) ) %>%
  arrange(desc(importance_score))

head(importance.df.resamptest, 10)

on.hoch5k <- make.predictions.df(test.classifier,
                                 hoch5k.adultDGCs.lognorm)

table(on.hoch5k$predict)
summary(on.hoch5k$Fos_pos)



### EXPLORING PENK EXPRESSION IN TEST DATA TO VALIDATE PUTATIVE ENGRAM CELLS

# > table(predictions.Hoch5k.lognorm$predict)
# 
# Fos- Fos+ 
#   1004   10 

for(bin in levels(df$prob_bin)){
  print(bin)
  print(sum(df$penk_count[df$prob_bin==bin]>quantile(df$penk_count,0.75)))
  print(sum(sum(df$penk_count[df$prob_bin==bin]>quantile(df$penk_count,0.75)))/sum(df$prob_bin==bin) )
}



putative.engramcells.idx <- which(predictions.Hoch5k.lognorm$predict=="Fos+")

hoch5k.adultDGCs.lognorm$Engramcell <- "Putative Inactive"
hoch5k.adultDGCs.lognorm$Engramcell[putative.engramcells.idx] <- "Putative Engram Cell"
hoch5k.adultDGCs.lognorm$Engramcell <- as.factor(hoch5k.adultDGCs.lognorm$Engramcell)
hoch5k.adultDGCs.lognorm$prob_engram <- predictions.Hoch5k.lognorm$Fos_pos

which(rownames(hoch5k.adultDGCs.lognorm)=="Arc")


# predictions <- as.data.frame(predict(downsamp.lognorm.classifier.forHoch5k, 
#                                      hoch5k.adultDGCs.lognorm[, 1:(length(hoch5k.adultDGCs.lognorm)-1)], 
#                                      type = "prob"))

levels(predictions.Hoch5k.lognorm$engramobserved) <- c(0,1)
#Plotting ROC...
roc.engramcell <- roc(downsamp.lognorm.predictions.forHoch5k$engramobserved, 
                      as.numeric(downsamp.lognorm.predictions.forHoch5k$Fos_pos) )
# there is an error here the predictions.Hoch5k.lognorm$engramobserved is showing only as 1 which cannot be true
# seomthing is wrong with the code don't know where this comes from

#roc.inactive <- roc(predictions$inactiveobserved, as.numeric(predictions$Fos_neg) )

jpeg("ROC_lognorm.jpg", width = 350, height = "350")
plot(roc.engramcell, col = "red", main = "ROC of RF Classifier")
dev.off()


#looking at penk expression in these cells
putative.engramcells.idx <- which(predictions.Hoch5k.lognorm$predict == "Fos+")
penkcckmalat1.idx <- which(colnames(hoch5k.adultDGCs.lognorm) %in% c("Cck","Penk", "Malat1") )
                           
dfpcm <- hoch5k.adultDGCs.lognorm[putative.engramcells.idx, penkcckmalat1.idx]
colnames(dfpcm) <- c("Cck","Penk", "Malat1")
sum(dfpcm$Malat1>0)

#to plot df
df <- data.frame(Metrics = rf.performances$lognorm_forHoch5k)
rownames(df) <- rownames(rf.performances)




##### A MILLION CLASSIFIERS TRIED AND FAILED




## Training our first classifier

binarized.counts <- data.frame( lapply(combined.counts, function(x) as.character(as.integer(x>0))) ) #binarize
binarized.counts <- data.frame( t(binarized.counts), stringsAsFactors = TRUE ) #convert the strings into factors
binarized.counts$Engramcell <- as.factor(combined.meta$fos_status)

# sample slit comes from library(caTools)
#Attempt 1 regular random forest with split
split <- sample.split(binarized.counts$Engramcell, SplitRatio = 0.7)

training_set = subset(binarized.counts, split == TRUE)
test_set = subset(binarized.counts, split == FALSE)

onehot.classifier = randomForest(x = training_set[-1],
                          y = training_set$Engramcell,
                          ntree = 500)



#testing the functions above
onehotjeager.rf.predictions <- make.predictions.df(onehot.classifier, test_set)

rf.performances <- data.frame( onehot = assessment(onehotjeager.rf.predictions) )
rownames(rf.performances) <- c("F1 Score", "AUC", "Precision", "Recall",
                               "FPR", "FNR", "True Positives", "False Negatives", 
                               "True Negatives", "False Positives")



#CROSS VALIDATION
# not sure what to do with the object the cross validation returns it can't really be used to predict anything
# rfcv is from rfUtilities package
rf.cv.classifier <- rf.crossValidation(onehot.classifier, training_set[-1], 
                            normalize = FALSE, p=0.1, 
                            n=10, ntree=501)


onehotCV.rf.predictions <- make.predictions.df(rf.cv.classifier)

rf.performances$morethanoneread <- assessment(morethanoneread.rf.predictions)


###DOWN SAMPLING
#uses sampler package https://www.rdocumentation.org/packages/sampler/versions/0.2.4

combined.meta$idx <- c(1:750)
df.temp <- ssamp(df=combined.meta[combined.meta$fos_status=="Fos+",], n=175, strata=treatment, over=0)
# due to representation we need to set this to n=175 to get 174 cells 

#make the new count matrix and meta data
downsamp.combinedcounts <- cbind(combined.counts[,df.temp$idx],
                                 combined.counts[,combined.meta$fos_status=="Fos-"]
                                 )


downsamp.meta <- rbind(df.temp,
                       combined.meta[combined.meta$fos_status=="Fos-",]
                       ) 


binarized.counts <- data.frame( lapply(downsamp.combinedcounts, function(x) as.character(as.integer(x>0))) ) #binarize
binarized.counts <- data.frame( t(binarized.counts), stringsAsFactors = TRUE ) #convert the strings into factors
binarized.counts$Engramcell <- as.factor(downsamp.meta$fos_status)

# sample slit comes from library(caTools)
#Attempt 1 regular random forest with split
split <- sample.split(binarized.counts$Engramcell, SplitRatio = 0.7)

training_set = subset(binarized.counts, split == TRUE)
test_set = subset(binarized.counts, split == FALSE)

downsamp.onehot.classifier = randomForest(x = training_set[,1:(length(training_set)-1)],
                                 y = training_set$Engramcell,
                                 ntree = 500)

downsamp.onehot.predictions <- make.predictions.df(downsamp.onehot.classifier, test_set)

rf.performances$Downsampled_onehot <- assessment(downsamp.onehot.predictions)
for(i in c(1:length(rf.performances$Downsampled_onehot)) ){
  print(rownames(rf.performances)[i])
  print(rf.performances$Downsampled_onehot[i])
}

roc.engramcell <- roc(downsamp.onehot.predictions$engramobserved, 
                      as.numeric(downsamp.onehot.predictions$Fos_pos) )

# try training on raw data as well and do resampling on both,
# test on the hochgerner dataset as well

#roc.inactive <- roc(predictions$inactiveobserved, as.numeric(predictions$Fos_neg) )

jpeg("ROCBinarized.jpg", width = 350, height = "350")
plot(roc.engramcell, col = "red", main = "ROC of RF Classifier")
dev.off()
#lines(roc.inactive, col = "blue")

#RAW COUNTS DOWNSAMPLED
df.temp <- ssamp(df=combined.meta[combined.meta$fos_status=="Fos+",], n=175, strata=treatment, over=0)
# due to representation we need to set this to n=175 to get 174 cells 
# combined.counts[,df.temp$idx]
# combined.counts[,combined.meta$fos_status=="Fos-"]

#make the new count matrix and meta data
downsamp.combinedcounts <- cbind(combined.counts[,df.temp$idx],
                                 combined.counts[,combined.meta$fos_status=="Fos-"]
                                 )

#creating meta data
downsamp.meta <- rbind(df.temp,
                       combined.meta[combined.meta$fos_status=="Fos-",]
                       ) 


downsamp.raw.counts <- data.frame( t(downsamp.combinedcounts) ) #convert the strings into factors
downsamp.raw.counts$Engramcell <- as.factor(downsamp.meta$fos_status)

# sample slit comes from library(caTools)
#Attempt 1 regular random forest with split
split <- sample.split(downsamp.raw.counts$Engramcell, SplitRatio = 0.7)

training_set = subset(downsamp.raw.counts, split == TRUE)
test_set = subset(downsamp.raw.counts, split == FALSE)

downsamp.raw.classifier = randomForest(x = training_set[,1:(length(training_set)-1)],
                                          y = training_set$Engramcell,
                                          ntree = 500)


downsamp.raw.predictions <- make.predictions.df(downsamp.raw.classifier, test_set)

rf.performances$Downsampled_rawcounts <- assessment(downsamp.raw.predictions)

#to make ROC plot
roc.engramcell <- roc(downsamp.raw.predictions$engramobserved, 
                      as.numeric(downsamp.onehot.predictions$Fos_pos) )




#LOG NORMALIZED DOWNSAMPLED

combined.meta$idx <- c(1:750)
df.temp <- ssamp(df=combined.meta[combined.meta$fos_status=="Fos+",], n=175, strata=treatment, over=0)
# due to representation we need to set this to n=175 to get 174 cells 

#make the new count matrix and meta data
downsamp.combinedcounts <- cbind(combined.counts[,df.temp$idx],
                                 combined.counts[,combined.meta$fos_status=="Fos-"]
                                 )


downsamp.meta <- rbind(df.temp,
                       combined.meta[combined.meta$fos_status=="Fos-",]
                       ) 




logplusone <- function(x){
  return( log(x+1) )
}


downsamp.lognorm.counts <- apply(downsamp.combinedcounts, MARGIN = 1, 
                                       FUN = logplusone
                                       )
downsamp.lognorm.counts  <- scale(t(downsamp.lognorm.counts )) #apply is transposing the data frame for some reason
downsamp.lognorm.counts  <- data.frame( t(downsamp.lognorm.counts ) ) #convert the strings into factors
downsamp.lognorm.counts$Engramcell <- as.factor(downsamp.meta$fos_status)

#split data
split <- sample.split(downsamp.lognorm.counts$Engramcell, SplitRatio = 0.7)

training_set = subset(downsamp.lognorm.counts, split == TRUE)
test_set = subset(downsamp.lognorm.counts, split == FALSE)

#make classifier
downsamp.lognorm.classifier = randomForest(x = training_set[,1:(length(training_set)-1)],
                                       y = training_set$Engramcell,
                                       ntree = 500)


downsamp.lognorm.predictions <- make.predictions.df(downsamp.lognorm.classifier, test_set)

rf.performances$Downsampled_lognorm <- assessment(downsamp.lognorm.predictions)

#Plotting ROC...
roc.engramcell <- roc(downsamp.lognorm.predictions$engramobserved, 
                      as.numeric(downsamp.lognorm.predictions$Fos_pos) )

#roc.inactive <- roc(predictions$inactiveobserved, as.numeric(predictions$Fos_neg) )

jpeg("ROC_lognorm.jpg", width = 350, height = "350")
plot(roc.engramcell, col = "red", main = "ROC of RF Classifier")
dev.off()

#With lower thresholds
downsamp.lognorm.predictions$prob




#RANK TRANSFORM


#INCLUDING CELLS FROM HABIB ET AL. 2016



#TRYING ON HOCHGERNER dataset
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

#get shared genes
shared.genes <- intersect(rownames(downsamp.combinedcounts),
                          rownames(hochgerner5k_2018_counts)
                          )

#Binarize and transpose hochgerner adult p35 cells
hoch5k.adultDGCs.onehot <- data.frame( lapply(hochgerner5k_2018_counts[, hoch5k.GC_Adult.p35.idx], function(x) as.character(as.integer(x>0))) ) #binarize
hoch5k.adultDGCs.onehot <- data.frame( t(hoch5k.adultDGCs.onehot), stringsAsFactors = TRUE ) #convert the strings into factors
colnames(hoch5k.adultDGCs.onehot) <- rownames(hochgerner5k_2018_counts)
#lastly we must filter for matching genes this will drop our list to 14296 genes
hoch5k.adultDGCs.onehot <- hoch5k.adultDGCs.onehot[,colnames(hoch5k.adultDGCs.onehot) %in% shared.genes]


#restrict genes to testing dataset and binarizing
# We are getting downsamp.combinedcoutns and downsamp.meta from the earlier DOWNSAMPLING section, we restrict to the same cells
binarized.counts.forHoch5k <- data.frame( lapply(downsamp.combinedcounts, function(x) as.character(as.integer(x>0))) ) #binarize
binarized.counts.forHoch5k <- data.frame( t(binarized.counts.forHoch5k), stringsAsFactors = TRUE ) #convert the strings into factors
colnames(binarized.counts.forHoch5k) <- rownames(downsamp.combinedcounts)[rownames(downsamp.combinedcounts) %in% colnames(hoch5k.adultDGCs.onehot)]
#I think the line above amy be unnecessary
#again we filter for matching genes in this line
binarized.counts.forHoch5k <- binarized.counts.forHoch5k[,colnames(binarized.counts.forHoch5k) %in% shared.genes]
binarized.counts.forHoch5k$Engramcell <- downsamp.meta$fos_status

binarized.counts.forHoch5k$data.use <- rep("train/validate",dim(binarized.counts.forHoch5k)[1]) 

#merge the datasets here then try training this classifier
hoch5k.adultDGCs.onehot$Engramcell <- rep("Fos+", dim(hoch5k.adultDGCs.onehot)[1])
hoch5k.adultDGCs.onehot$data.use <- rep("test",dim(hoch5k.adultDGCs.onehot)[1]) 

all.data <- bind_rows(binarized.counts.forHoch5k, hoch5k.adultDGCs.onehot)
index <- 1:ncol(all.data)
all.data[ , index] <- lapply(all.data[ , index], as.factor)

binarized.counts.forHoch5k <- all.data[all.data$data.use == "train/validate", ]
binarized.counts.forHoch5k <- binarized.counts.forHoch5k[,1:(length(binarized.counts.forHoch5k)-1)]
hoch5k.adultDGCs.onehot <- all.data[all.data$data.use == "test", ]
hoch5k.adultDGCs.onehot <- hoch5k.adultDGCs.onehot[1:(length(hoch5k.adultDGCs.onehot)-1)]


# train classifier on restricted genes
split <- sample.split(binarized.counts.forHoch5k$Engramcell, SplitRatio = 0.7)

training_set = subset(binarized.counts.forHoch5k, split == TRUE)
validation_set = subset(binarized.counts.forHoch5k, split == FALSE)

#train
downsamp.onehot.classifier.forHoch5k = randomForest(x = training_set[,1:(length(training_set)-1)],
                                          y = training_set$Engramcell,
                                          ntree = 500)

importance.df <- data.frame(gene = as.character( rownames(downsamp.onehot.classifier.forHoch5k$importance) ),
                         importance_score = as.numeric( downsamp.onehot.classifier.forHoch5k$importance ) ) %>%
  arrange(desc(importance_score))
  

important.df <- sort(importance.df, decreasing = TRUE) # these gene's are wierd, I do not think this classifer is working well at all

#validate and measure performance
downsamp.onehot.predictions <- make.predictions.df(downsamp.onehot.classifier.forHoch5k, validation_set)

rf.performances$Onehot.downsampled.forHoch5k <- assessment(downsamp.onehot.predictions)
rf.performances #notice that as the number of genes drops our classifiers performance drops quite badly as well

#classify new data
predictions.Hoch5k <- make.predictions.df(downsamp.onehot.classifier.forHoch5k, hoch5k.adultDGCs.onehot)





### JEAGER DEGS ONLY LOGNORM 


Jeager.DEGs <- read.csv("Jeager2018_GSE98679/jeager_DEGs_grouped.csv", header = T)

colnames(Jeager.DEGs)[1] <- "GeneID"

#restrict to matching genes
shared.genes.jdegs <- multi.intersect( list(colnames(hoch5k.adultDGCs.lognorm), 
                                colnames(downsamp.lognorm.counts),
                                Jeager.DEGs$GeneID
                                )# close list
                           )#close multi.intersect

hoch5k.adultDGCs.lognorm.jdegs  <- hoch5k.adultDGCs.lognorm[,colnames(hoch5k.adultDGCs.lognorm) %in% shared.genes.jdegs]
hoch5k.adultDGCs.relu.jdegs <- data.frame(t(apply(hoch5k.adultDGCs.lognorm.jdegs, MARGIN =1, FUN = relu)))
colnames(hoch5k.adultDGCs.relu.jdegs) <- colnames(hoch5k.adultDGCs.lognorm.jdegs)
rownames(hoch5k.adultDGCs.relu.jdegs) <- rownames(hoch5k.adultDGCs.lognorm.jdegs)


downsamp.lognorm.counts.jdegs <- downsamp.lognorm.counts[,colnames(downsamp.lognorm.counts) %in% shared.genes.jdegs]


hoch5k.adultDGCs.lognorm.jdegs$Engramcell <- as.factor(rep("Fos+", dim(hoch5k.adultDGCs.lognorm.jdegs)[1]))
#hoch5k.adultDGCs.lognorm.jdegs$data.use <- rep("test",dim(hoch5k.adultDGCs.lognorm.jdegs)[1]) 

downsamp.lognorm.counts.jdegs$Engramcell <- as.factor( downsamp.meta$fos_status )
#downsamp.lognorm.counts.jdegs$data.use <- rep("train/validate",dim(downsamp.lognorm.counts.jdegs)[1]) 



#split data
split <- sample.split(downsamp.lognorm.counts.jdegs$Engramcell, SplitRatio = 0.7)

training_set = subset(downsamp.lognorm.counts.jdegs, split == TRUE)
validation_set = subset(downsamp.lognorm.counts.jdegs, split == FALSE)

#make classifier
downsamp.lognorm.classifier.forHoch5k.jdegs = randomForest(x = training_set[,1:(length(training_set)-1)],
                                                     y = training_set$Engramcell,
                                                     ntree = 500)


downsamp.lognorm.predictions.forHoch5k.jdegs <- make.predictions.df(downsamp.lognorm.classifier.forHoch5k.jdegs,
                                                              validation_set)

rf.performances$lognorm_forHoch5k_jdegs <- assessment(downsamp.lognorm.predictions.forHoch5k.jdegs)

#classify new data
hoch5k.adultDGCs.lognorm.jdegs[is.na(hoch5k.adultDGCs.lognorm.jdegs)] <- 0


predictions.Hoch5k.lognorm.jdegs <- make.predictions.df(downsamp.lognorm.classifier.forHoch5k.jdegs, 
                                                  hoch5k.adultDGCs.lognorm.jdegs)

importance.df.jdegs <- data.frame(gene = as.character( rownames(downsamp.lognorm.classifier.forHoch5k.jdegs$importance) ),
                            importance_score = as.numeric( downsamp.lognorm.classifier.forHoch5k.jdegs$importance ) ) %>%
  arrange(desc(importance_score))




### RESAMPLING, LOGNORM TRANSFORM
# so this particular method tries to ensure the classes are
# of the same covarience matrix, this is an interesting approach just from that

# normalization here: http://bioconductor.org/books/3.13/OSCA.basic/normalization.html
# https://bioconductor.org/packages/3.13/data/experiment/html/scRNAseq.html
# Consider renormalizing based on simulated resampling with proportions based on Ca2+ studies
# 
# Also this python method should be considered
# https://imbalanced-learn.org/stable/references/generated/imblearn.ensemble.BalancedRandomForestClassifier.html
# https://statistics.berkeley.edu/sites/default/files/tech-reports/666.pdf

# mean center scale then log(x+1)

logplusone <- function(x){
  return( log(x+1) )
}

#Binarize and transpose hochgerner adult p35 cells
hoch5k.adultDGCs.lognorm <- apply(hochgerner5k_2018_counts[, hoch5k.GC_Adult.p35.idx], MARGIN = 1, FUN = as.integer)
hoch5k.adultDGCs.lognorm <- apply(hoch5k.adultDGCs.lognorm,
                                  MARGIN = 1,
                                  FUN = logplusone
)
hoch5k.adultDGCs.lognorm  <- scale( t(hoch5k.adultDGCs.lognorm ) ) #apply is transposing the data frame for some reason
hoch5k.adultDGCs.lognorm  <- as.data.frame(hoch5k.adultDGCs.lognorm )
#lastly we must filter for matching genes this will drop our list to 14296 genes


#doing lognormalization for the training data
resamp.lognorm.counts <- apply(combined.counts, 
                                 MARGIN = 1, 
                                 FUN = logplusone
                               )
resamp.lognorm.counts  <- scale( t(resamp.lognorm.counts) ) #we need it transposed so that the scaling is done per gene not cell
resamp.lognorm.counts  <- data.frame( t(resamp.lognorm.counts) ) 


#restrict to matching genes
shared.genes <- intersect( colnames(hoch5k.adultDGCs.lognorm), colnames(resamp.lognorm.counts) )
hoch5k.adultDGCs.lognorm  <- hoch5k.adultDGCs.lognorm[,colnames(hoch5k.adultDGCs.lognorm) %in% shared.genes]
resamp.lognorm.counts <- resamp.lognorm.counts[,colnames(resamp.lognorm.counts) %in% shared.genes]

hoch5k.adultDGCs.lognorm$Engramcell <- rep("Fos+", dim(hoch5k.adultDGCs.lognorm)[1])
hoch5k.adultDGCs.lognorm$data.use <- rep("test",dim(hoch5k.adultDGCs.lognorm)[1]) 

resamp.lognorm.counts$Engramcell <- as.factor( combined.meta$fos_status )
resamp.lognorm.counts$data.use <- rep("train/validate", dim(resamp.lognorm.counts)[1]) 



all.data <- bind_rows(resamp.lognorm.counts, hoch5k.adultDGCs.lognorm)
all.data$Engramcell <- as.factor(all.data$Engramcell)
all.data$data.use <- as.factor(all.data$data.use)

resamp.lognorm.counts <- all.data[all.data$data.use == "train/validate", ]
resamp.lognorm.counts <- resamp.lognorm.counts[,1:(length(resamp.lognorm.counts)-1)]
hoch5k.adultDGCs.lognorm <- all.data[all.data$data.use == "test", ]
hoch5k.adultDGCs.lognorm <- hoch5k.adultDGCs.lognorm[,1:(length(hoch5k.adultDGCs.lognorm)-1)]


#split data
split <- sample.split(resamp.lognorm.counts$Engramcell, SplitRatio = 0.7)

training_set = subset(resamp.lognorm.counts, split == TRUE)
validation_set = subset(resamp.lognorm.counts, split == FALSE)

#train the calssifier
# see rfUtilities documentation here https://cran.r-project.org/web/packages/rfUtilities/rfUtilities.pdf
resampled.lognorm.classifier.forHoch5k <- rf.classBalance( ydata=training_set$Engramcell, 
                                                           xdata=training_set[,1:(length(training_set)-1)],
                                                           cbf = 1, #tests to see if there is an imbalance
                                                           sf = 1 # Majority subsampling factor, sf=1 gives  perfect balanced 
                                                           ) 

# 
downsamp.lognorm.predictions.forHoch5k.jdegs <- make.predictions.df(resampled.lognorm.classifier.forHoch5k$model,
                                                                    validation_set)

rf.performances$resampling <- assessment(downsamp.lognorm.predictions.forHoch5k.jdegs)


predictions.Hoch5k.lognorm <- predict(resampled.lognorm.classifier.forHoch5k$model,
                                      newdata = hoch5k.adultDGCs.lognorm,
                                      type = "prob")

predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)] 


resampled.importance.df <- data.frame(gene = as.character( rownames(resampled.lognorm.classifier.forHoch5k$model$importance) ),
                            importance_score = as.numeric( resampled.lognorm.classifier.forHoch5k$model$importance ) ) %>%
  arrange(desc(importance_score))


hoch5k.adultDGCs.lognorm[is.na(hoch5k.adultDGCs.lognorm)] <- 0
predictions <- as.data.frame(predict(resampled.lognorm.classifier.forHoch5k$model, 
                                     hoch5k.adultDGCs.lognorm[, 1:(length(hoch5k.adultDGCs.lognorm)-1)], 
                                     type = "prob")
                             )


downsamp.lognorm.predictions.forHoch5k.jdegs <- make.predictions.df(resampled.lognorm.classifier.forHoch5k$model,
                                                                    hoch5k.adultDGCs.lognorm[, 1:(length(hoch5k.adultDGCs.lognorm)-1)])


#Eventually when we make rpedictions...
#predict(rf.resampled.classifier$model,  hochgernerdata, type = "prob")

# model assessment from tutorial....
# # Calculate Kappa for each balanced model in ensemble 
# for(i in 1:length(cb$confusion) ) { 
#   print( accuracy(cb$confusion[[i]][,1:2])[5] ) 
# }
# 
# # Evaluate cumulative and mean confusion matrix
# accuracy( round((cb$confusion[[1]] + cb$confusion[[2]] + cb$confusion[[3]]))[,1:2] )
# accuracy( round((cb$confusion[[1]] + cb$confusion[[2]] + cb$confusion[[3]])/3)[,1:2])
# 








###  TRY RELU ACTIVATION
library(sigmoid)

relu.downsamp.combinedcounts <-

# We are getting downsamp.combinedcoutns and downsamp.meta from the earlier DOWNSAMPLING section, we restrict to the same cells
 #binarize
binarized.counts.forHoch5k <- data.frame( t(binarized.counts.forHoch5k), stringsAsFactors = TRUE ) #convert the strings into factors
colnames(binarized.counts.forHoch5k) <- rownames(downsamp.combinedcounts)[rownames(downsamp.combinedcounts) %in% colnames(hoch5k.adultDGCs.onehot)]
#again we filter for matching genes in this line
binarized.counts.forHoch5k <- binarized.counts.forHoch5k[,colnames(binarized.counts.forHoch5k) %in% shared.genes]
binarized.counts.forHoch5k$Engramcell <- downsamp.meta$fos_status



###  Try Restricting genes to jeager DEGs








# > table(predictions.Hoch5k$predict)
# 
# Fos- 
#   1014 
# > summary(predictions.Hoch5k$Inactive)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.6120  0.7100  0.7360  0.7387  0.7680  0.8860 
# > summary(predictions.Hoch5k$Engram)  
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1140  0.2320  0.2640  0.2613  0.2900  0.3880 




#VISUALIZE RANDOM FOREST
# A guide here for all thign random forest related has this 
# https://rpubs.com/markloessi/498787


predictions$observed <- test_set$Engramcell
colnames(predictions)[1:2] <- c("Fos_neg","Fos_pos")
predictions$engramobserved <- ifelse(predictions$observed=="Fos+", 1, 0)
predictions$inactiveobserved <- ifelse(predictions$observed=="Fos-", 1, 0)







#Plotting in ggplot https://www.statology.org/roc-curve-ggplot2/
predictions.rf <- as.data.frame(predict(classifier, test_set, type = "prob"))
predictions.rf$predict <- names(predictions.rf)[1:2][apply(predictions.rf[,1:2], 1, which.max)] #1:2 for the number of classes
predictions.rf$observed <- test_set$Engramcell
colnames(predictions.rf)[1:2] <- c("Fos_neg","Fos_pos")
predictions.rf$engramobserved <- ifelse(predictions.rf$observed=="Fos+", 1, 0)

predictions.rf$inactiveobserved <- ifelse(predictions.rf$observed=="Fos-", 1, 0)

roc.engramcell <- roc(predictions.rf$engramobserved, as.numeric(predictions.rf$Fos_pos) )
roc.inactive <- roc(predictions.rf$inactiveobserved, as.numeric(predictions.rf$Fos_neg) )

plot(roc.engramcell, col = "red")
lines(roc.inactive, col = "blue")



