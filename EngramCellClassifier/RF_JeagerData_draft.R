#Designed to run ont he pavlab server

# Load libraries
library(tidyverse)
library(dplyr)
#library(AnnotationDbi)
library(randomForest)
library(rfUtilities) # tutorial here: https://evansmurphy.wixsite.com/evansspatial/random-forest-sdm
#library(data.table)
#library(reshape2)
library(FactoMineR)
library(factoextra)
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

# Loading Lacar et al., (2016)
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
)
)

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
lacar2016_snNE_counts[rownames(lacar2016_snNE_counts) %in% shared.genes,]
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
#tibble::column_to_rownames(var="gene") %>% this throws an error, it would be good
#dplyr::select(-gene)

#this join is possibly inclduing the gene rows leading to mismathc number of cells later 

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


## These functions will save the assesments of the classifiers
make.predictions.df <- function(classifier.object){
  #generate predictions for making classifier summary
  predictions <- as.data.frame(predict(classifier.object, test_set, type = "prob"))
  predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)] #1:2 for the number of classes
  predictions$observed <- test_set$Engramcell
  colnames(predictions)[1:2] <- c("Fos_neg","Fos_pos")
  predictions$engramobserved <- ifelse(predictions$observed=="Fos+", 1, 0)
  predictions$inactiveobserved <- ifelse(predictions$observed=="Fos-", 1, 0)
  return(predictions)
}



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


#testing the functions above
onehotjeager.rf.predictions <- make.predictions.df(onehot.classifier)

rf.performances <- data.frame( onehot = assessment(onehotjeager.rf.predictions) )
rownames(rf.performances) <- c("F1 Score", "AUC", "Precision", "Recall",
                               "FPR", "FNR", "True Positives", "False Negatives", 
                               "True Negatives", "False Positives")

#TRYING WITH MORE THAN ONE READ AS THRESHOLD RATHER THAN 0
#
binarized.counts <- data.frame( lapply(combined.counts, function(x) as.character(as.integer(x>1))) ) #binarize
binarized.counts <- data.frame( t(binarized.counts), stringsAsFactors = TRUE ) #convert the strings into factors
binarized.counts$Engramcell <- as.factor(combined.meta$fos_status)

# sample slit comes from library(caTools)
#Attempt 1 regular random forest with split
split <- sample.split(binarized.counts$Engramcell, SplitRatio = 0.7)

training_set = subset(binarized.counts, split == TRUE)
test_set = subset(binarized.counts, split == FALSE)

morethanone.classifier = randomForest(x = training_set[-1],
                          y = training_set$Engramcell,
                          ntree = 500)

morethanoneread.rf.predictions <- make.predictions.df(classifier)

rf.performances$morethanoneread <- assessment(morethanoneread.rf.predictions)

#CROSS VALIDATION
# rfcv is from rfUtilities package
rf.cv.classifier <- rf.crossValidation(onehot.classifier, training_set[-1], 
                            normalize = FALSE, p=0.1, 
                            n=10, ntree=501)


onehotCV.rf.predictions <- make.predictions.df(rf.cv.classifier)

rf.performances$morethanoneread <- assessment(morethanoneread.rf.predictions)


###DOWN SAMPLING
#uses sampler package https://www.rdocumentation.org/packages/sampler/versions/0.2.4

combined.meta$idx <- c(1:750)
df.temp <- ssamp(df=combined.meta[combined.meta$fos_status=="Fos+",], n=176, strata=treatment, over=0)


downsamp.combinedcounts <- vbind()


#RANK TRANSFORM


#INCLUDING CELLS FROM HABIB ET AL. 2016



#FORNITO FUETZ









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



