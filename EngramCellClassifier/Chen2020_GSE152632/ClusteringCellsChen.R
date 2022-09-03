# 
# 
# 
# This is the cell clustering an cleaning for Chen et al., (2020)
# 
# Chen, M. B., Jiang, X., Quake, S. R., & Südhof, T. C. (2020).
# Persistent transcriptional programmes are associated with remote memory.
# Nature, 587(7834), 437-442.

# Engrma labelling proof of "maturation" of prefrontal engrams
# Kitamura, T., Ogawa, S. K., Roy, D. S., Okuyama, T., Morrissey, M. D., Smith, L.
# M., ... & Tonegawa, S. (2017). Engrams and circuits crucial for systems 
# consolidation of a memory. Science, 356(6333), 73-78.
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5493329/

# https://www.nature.com/articles/s41467-019-10683-2#Sec2
# replay during sleep following learning and retrival in engram cells
# compared to non engram cells from the initial novel context exploration
# Ghandour, K., Ohkawa, N., Fung, C. C. A., Asai, H., Saitoh, Y., Takekawa, T., 
# ... & Inokuchi, K. (2019). Orchestrated ensemble activities constitute a 
# hippocampal memory engram. Nature communications, 10(1), 1-14.

# 
# El-Boustani, S., Ip, J. P., Breton-Provencher, V., Knott, G. W., Okuno, H., 
# Bito, H., & Sur, M. (2018). Locally coordinated synaptic plasticity of visual
# cortex neurons in vivo. Science, 360(6395), 1349-1354.
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6366621/
# 

# https://www.nature.com/articles/s41593-022-01041-5

#Path on pavlab server
#setwd("/home/acampbell/PavLabEngrams/EngramCellClassifier")

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

#adding cell bacrcodes from coutn data to rows of metadata for seurat
rownames(chen2020_meta) <- colnames(chen2020_counts)


# 
# Briefly, genes that appeared in fewer than five cells, samples with fewer than 
# 100 genes and samples with less than 50,000 reads were excluded from the analysis. 
# Out of these cells, those with more than 30% of reads as ERCC, and more than 10% 
#   mitochondrial or 10% ribosomal were also excluded from analysis. Counts were 
#  log-normalized and then scaled where appropriate.
#

#Normal Seurat Workflow
#from here: https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
chen <- CreateSeuratObject(counts = chen2020_counts,  min.cells = 5, min.features = 100 )
chen[["percent.mt"]] <- PercentageFeatureSet(chen, pattern = "^MT-")
chen[["percent.ribo"]] <- PercentageFeatureSet(chen, pattern = "^(MRP|RP)") 
#checks for genes beginning in MRP or RP, https://www.genenames.org/data/genegroup/#!/group/1054

#selecting variable genes and normalizing
chen <- subset(chen, subset = nCount_RNA > 50000 & nFeature_RNA > 100 & percent.mt < 10 & percent.ribo < 10)
chen <- NormalizeData(chen, normalization.method = "LogNormalize", scale.factor = 10000)
chen <- FindVariableFeatures(chen, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(chen), 10)
top10

dev.off()
jpeg("varfeatChen.jpg", width = 350, height = 350)
VariableFeaturePlot(chen)
dev.off()

#scaling
all.genes <- rownames(chen)
chen <- ScaleData(chen, features = all.genes)
chen <- RunPCA(chen, features = VariableFeatures(object = chen))


dev.off()
jpeg("pcloadingsChen.jpg", width = 350, height = 350)
VizDimLoadings(chen, dims = 1:2, reduction = "pca")
dev.off()


dev.off()
jpeg("DimPlotPCA_Chen.jpg", width = 350, height = 350)
DimPlot(chen, reduction = "pca")
dev.off()


dev.off()
jpeg("HeatmapPCA_Chen.jpg", width = 350, height = 350)
DimHeatmap(chen, dims = 1, cells = 500, balanced = TRUE)
dev.off()

chen <- JackStraw(chen, num.replicate = 100)
chen <- ScoreJackStraw(chen, dims = 1:20)

dev.off()
jpeg("ElbowPlotPCA_Chen.jpg", width = 350, height = 350)
JackStrawPlot(chen, ndims = 1:20)
dev.off()

dev.off()
jpeg("ElbowPlotPCA_Chen.jpg", width = 350, height = 350)
ElbowPlot(chen)
dev.off()

chen <- RunTSNE(chen, dim.embed = 2)
chen <- AddMetaData(chen, chen2020_meta)
#clustering usign the top 20 PCs as they did in Chen
chen <- FindNeighbors(chen, reduction = "tsne", dims = 1:2)
chen <- FindClusters(chen, 
                     resolution = 0.03) #setting this parameter to 0.03 gets 6 clusters, at default 0.5 we find 21

#Plotting the clusters
dev.off()
jpeg("DimplotTSNE_Chen.jpg", width = 700, height = 700)
DimPlot(chen, reduction = "tsne")
dev.off()



#Idendtifying the big ones, Neurons, astrocytes, glial cells
# DotPlot may be good for later when comparing between groups
# https://satijalab.org/seurat/reference/dotplot

#Finding Markers
chen.markers <- FindAllMarkers(chen, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)


#checking for the markers described in chen
test <- chen.markers %>%
  group_by(cluster)

test <- data.frame(test)
test <- test[,test$p_val_adj<0.05]
layer.markers <- c("Dkkl1","Rprm","Calb2","Tesc","Tnfaip8l3","Tshz2","Lhx6")
layer.markers.idx <- which(test$gene %in% layer.markers)
test$cluster[layer.markers.idx]

#they are there and match the cluster assignments from chen
# > test$cluster[layer.markers.idx]
# [1] 0 1 2 3 4 5 6
# Levels: 0 1 2 3 4 5 6

test$p_val_adj[layer.markers.idx]
#The genes are significant even after bonferonni correction
# > test$p_val_adj[layer.markers.idx]
# [1]  0.000000e+00  0.000000e+00  0.000000e+00 6.459740e-244 1.979667e-140
# [6] 8.827732e-301  0.000000e+00




top2.chen.markers <- chen.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

top2.chen.markers <-  data.frame(top2.chen.markers)
write.csv(chen.markers, "Chen2020_ClusterMarkers.csv")

top5 <- chen.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)
top5 <-  data.frame(top5)


# Zeisel, A., Muñoz-Manchado, A. B., Codeluppi, S., Lönnerberg, P., La Manno, G., Juréus, A., ... & Linnarsson, S.
# (2015). Cell types in the mouse cortex and hippocampus revealed by single-cell RNA-seq. 
# Science, 347(6226), 1138-1142.

# Layer specific markers from above publication can be found in Chen data.  


#making list of top two markers to display on clustering graph
two.mrkrs <- c() # the list we will use to label the graph
j=0
for(i in c(1:dim(top2.chen.markers)[1]) ){
  j = j+1
  if(j==1){
    first.marker <- top2.chen.markers$gene[i]
  }# end of if statement

  if(j == 2){
    # glue the markers together and...
    this.clstr <- paste(as.character(first.marker),
          as.character(top2.chen.markers$gene[i]),
          sep = "-")
    #...add to the list as a single entry
    two.mrkrs <- c(two.mrkrs, this.clstr)
    j = 0 # reset j
  }# end of if statement
}# end of for loop

# > two.mrkrs
# [1] "Rprm-Crym"            "Dkkl1-Ptn"            "Rgs4-Slc30a3"        
# [4] "Myl4-Lypd1"           "Prkcg-Bdnf"           "3110035E14Rik-Rprm"  
# [7] "Synpr-Vip"            "Nnat-Rasl10a"         "Vip-Penk"            
# [10] "Rgs4-Nrn1"            "Sst-Crhbp"            "Ndnf-Gad1"           
# [13] "Synpr-Crh"            "3110035E14Rik-Cxcl12" "Tac2-Npy"            
# [16] "Sub1-Atpif1"          "Cnr1-Tac2"            "Rgs4-Nrep"           
# [19] "Vip-Tac2"             "Lamp5-Wfs1"           "Pcp4-Tmsb10"

#renaming cluster id's
names(two.mrkrs) <- levels(chen)
chen <- RenameIdents(chen, two.mrkrs)

#Visualize Clusters
dev.off()
jpeg("DimPlotTSNE_Chen_7clusters.jpg", width = 700, height = 700)
DimPlot(chen, reduction = "tsne", label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()


#clusters vs Snap25
dev.off()
jpeg("ClustersVsSnap25.jpg", width = 1400, height = 700)
DimPlot(chen, reduction = "tsne") + FeaturePlot(chen, reduction =  "tsne", feature = c("Snap25") )
dev.off()

#clusters vs engram label
dev.off()
jpeg("CellClustersTdplus.jpg", width = 1400, height = 700)
DimPlot(chen, reduction = "tsne", 
        label = TRUE, repel = TRUE, label.size = 5,
        pt.size = 2.5,
        shape.by = "engram_label") 
dev.off()

#clusters vs engram label
dev.off()
jpeg("EngramandConditionLabels.jpg", width = 1400, height = 700)
DimPlot(chen, reduction = "tsne", 
        pt.size = 2.5,
        shape.by = "condition_label",
        group.by = "engram_label") 
dev.off()

#Broad cell classes
dev.off()
jpeg("MajorCellClasses.jpg", width = 1400, height = 700)
FeaturePlot(chen, feature = c("Slc17a7","Vip","Pvalb","Sst") )
dev.off()

#Broad cell classes
dev.off()
jpeg("test.jpg", width = 1400, height = 700)
FeaturePlot(chen, feature = c("Snap25") )
dev.off()



#Running the classifier
library(randomForest)
library(rfUtilities)
library(tictoc)
library(caTools)
library(pROC)
library(stats)




## These functions will save the assesments of the classifiers
#this need changing, specifically test_set should not be a global, it is begging for issues

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
      assess <- assessment( pred , label = levels(data$Engramcell) )
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



### RUNNING THE SCRIPTS
# cellmarkers.col <- as.factor(chen@active.ident)
#write.csv(cellmarkers.col, "Chen2020_GSE152632/Chen2020_cellclusterlabels.csv")
markers <- read.csv("Chen2020_GSE152632/Chen2020_ClusterMarkers.csv")

#in this case we will need to ensure the training data is recieving the same proportion of cell types
#in each batch

#renaming chen clusters according to the labels in chen2020
#note its not clear how they chose these genes they are meant to be
#laminar markers from Ziesel et al., (2015)
layer.markers <- c("Dkkl1","Rprm","Calb2","Tesc","Tnfaip8l3","Tshz2","Lhx6")
names(layer.markers) <- levels(chen)
chen <- RenameIdents(chen, layer.markers)

library(sampler)
#creating chen data

log.chen <- log.norm(chen2020_counts)
log.chen$Engramcell <- as.factor(chen2020_meta$engram_label)
levels(log.chen$Engramcell) <- c("Neg", "Pos")
cellmarkers.col <- read.table("Chen2020_GSE152632/Chen2020_cellclusterlabels.csv")
#log.chen$cluster <- as.factor(cellmarkers.col[,2])
#remove TdTom-transgene, as it basically serves as a label that should not be present
log.chen <- log.chen[,!(names(log.chen) %in% "TdTom-transgene")]


#may include this later
#log.chen$experiment_group <- as.factor(chen2020_meta$source_name)

split <- sample.split(log.chen$Engramcell, SplitRatio = 0.7)

training_set.test = subset(log.chen, split == TRUE)
validation_set.test = subset(log.chen, split == FALSE)

#for loop
#outputs
tic()
classifier.mpfc <- resampled.randomForest.crossvalidated( data= labeled.data,
                                                          under.represented.class = "Td-",
                                                          over.represented.class = "Td+",
                                                          trees.total = 1000,
                                                          folds = 10,
                                                          proportion.each.batch=0.8,
                                                          batches.per.fold=20)
toc()

test.predictions <- make.predictions.df(classifier.mpfc, validation_set.test)
test.predictions$predict <- as.factor(test.predictions$predict)
test.predictions$engramobserved <- as.factor(test.predictions$engramobserved)
levels(test.predictions$engramobserved) <- levels(test.predictions$predict)
#instantiating rf performances for this script 
rf.performances <- data.frame( chen_resampled = assessment(test.predictions) )
rownames(rf.performances) <- c("F1 Score", "AUC", "Precision", "Recall",
                               "FPR", "FNR", "True Positives", "False Negatives", 
                               "True Negatives", "False Positives")

importance.df.resamptest <- data.frame(gene = as.character( rownames(test.classifier$importance) ),
                                       importance_score = as.numeric( test.classifier$importance ) ) %>%
  arrange(desc(importance_score))

head(importance.df.resamptest, 10)

#note we need to remove the TdTom-transgene from the gene set

#roc curve
levels(predictions.test$engramobserved) <- c(0,1)
#Plotting ROC...
roc.engramcell <- roc(test.predictions$engramobserved, 
                      as.numeric(test.predictions$Pos) )

# there is an error here the predictions.Hoch5k.lognorm$engramobserved is showing only as 1 which cannot be true
# seomthing is wrong with the code don't know where this comes from

#roc.inactive <- roc(predictions$inactiveobserved, as.numeric(predictions$Fos_neg) )

dev.off()
jpeg("ROC_lognormChen2020.jpg", width = 350, height = "350")
plot(roc.engramcell, col = "red", main = "ROC of Chen RF Classifier")
dev.off()

# Dropping threshold
assessment.thresholded <- function(predictions.df, thresh){
  # returns a vector of assessments to be used to make dataframe summarizing classifiers performance
  # can be used to make df of all calssifiers trained in a single run
  TP <- sum((predictions.df$Pos > thresh )&(predictions.df$observed == "Pos"))
  TN <- sum((predictions.df$Pos < thresh )&(predictions.df$observed == "Neg"))
  FN <- sum((predictions.df$Pos > thresh )&(predictions.df$observed == "Pos"))
  FP <- sum((predictions.df$Pos < thresh )&(predictions.df$observed == "Neg"))
  
  #precision and recall as well as sumamry stats F1
  precision <- TP/(TP+FP)
  recall <- TP/(TP+FN)
  F1.score = 2 * (precision * recall) / (precision + recall)
  FPR <- FP/(TN+FP)
  FNR <- FN/(TP+FN)
  
  #getting auc
  roc.engramcell <- roc(predictions.df$engramobserved, as.numeric(predictions.df$Pos) )
  AUC <- auc(roc.engramcell)
  
  return( c(F1.score, AUC, precision, recall, FPR, FNR,
            TP, FN, TN, FP) )
}

assessment.thresholded(test.predictions, 0.55)



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

#normalizing data with their cell types
cellmarkers.col <- read.csv("Chen2020_GSE152632/Chen2020_cellclusterlabels.csv")
log.chen.withintypes <- celltype.lognorm( countsdata = chen2020_counts,
                                         celltype.labels = cellmarkers.col[,2])

log.chen.withintypes$cluster <- as.factor(cellmarkers.col[,2])
log.chen.withintypes$Engramcell <- as.factor(chen2020_meta$engram_label)
levels(log.chen.withintypes$Engramcell) <- c("Neg", "Pos")
#remove TdTom-transgene, as it basically serves as a label that should not be present
# appears to have already been removed
#log.chen <- log.chen[,!(names(log.chen) %in% "TdTom-transgene")]


#may include this later
#log.chen$experiment_group <- as.factor(chen2020_meta$source_name)


#for loop
#outputs


split <- sample.split(log.chen.withintypes$Engramcell, SplitRatio = 0.7)

train = subset(log.chen.withintypes, split == TRUE)
test = subset(log.chen.withintypes, split == FALSE)

classifier.mpfc.withCT <- randomForest(x = train[,1:(length(train)-1)],
                         y = train$Engramcell,
                         ntree = 1000)



roc.engramcell.mpfc = roc(train, classifier.mpfc.withCT$votes[,2], plot=TRUE,
                          legacy.axes=TRUE, percent=TRUE,
                          xlab="False Positive Percentage", ylab="True Postive Percentage", 
                          col="firebrick4", lwd=4, print.auc=TRUE)


dev.off()
jpeg("ROCBinarized.jpg", width = 700, height = 700)
plot(roc.engramcell.mpfc, main = "ROC of mPFC RF Classifier")
dev.off()




# this stuff will take time to get working
test.predictions <- make.predictions.df(classifier.mpfc.withCT, test)
test.predictions$predict <- as.factor(test.predictions$predict)
test.predictions$engramobserved <- as.factor(test.predictions$engramobserved)
levels(test.predictions$engramobserved) <- levels(test.predictions$predict)
#instantiating rf performances for this script 
rf.performances <- data.frame( chen_resampled = assessment(test.predictions) )
rownames(rf.performances) <- c("F1 Score", "AUC", "Precision", "Recall",
                               "FPR", "FNR", "True Positives", "False Negatives", 
                               "True Negatives", "False Positives")

importance.df.resamptest <- data.frame(gene = as.character( rownames(test.classifier$importance) ),
                                       importance_score = as.numeric( test.classifier$importance ) ) %>%
  arrange(desc(importance_score))

head(importance.df.resamptest, 10)




## crossvalidation is not working for some reason
tic()
classifier.mpfc.withCT <- resampled.randomForest.crossvalidated( data= log.chen.withintypes,
                                                          under.represented.class = "Neg",
                                                          over.represented.class = "Pos",
                                                          trees.total = 1000,
                                                          folds = 5,
                                                          proportion.each.batch=0.8,
                                                          batches.per.fold=20)
toc()



classifier.mpfc.withCT$votes <- predict(object = classifier.mpfc.withCT ,
                            newdata = log.chen.withintypes,
                            type = 'vote', norm.votes = FALSE)

importance.mPFC.df.resamptest <- data.frame(gene = as.character( rownames(classifier.mpfc.withCT$importance) ),
                                       importance_score = as.numeric(classifier.mpfc.withCT$importance ) ) %>%
  arrange(desc(importance_score))

head(importance.mPFC.df.resamptest, 10)


roc.engramcell.mpfc = roc(log.chen.withintypes$Engramcell, classifier.mpfc.withCT$votes[,2], plot=TRUE,
                     legacy.axes=TRUE, percent=TRUE,
                     xlab="False Positive Percentage", ylab="True Postive Percentage", 
                     col="firebrick4", lwd=4, print.auc=TRUE)


dev.off()
jpeg("ROCBinarized.jpg", width = 700, height = 700)
plot(roc.engramcell.mpfc, main = "ROC of mPFC RF Classifier")
dev.off()





test.predictions <- make.predictions.df(classifier.mpfc, validation_set.test)
test.predictions$predict <- as.factor(test.predictions$predict)
test.predictions$engramobserved <- as.factor(test.predictions$engramobserved)
levels(test.predictions$engramobserved) <- levels(test.predictions$predict)
#instantiating rf performances for this script 
rf.performances <- data.frame( chen_resampled = assessment(test.predictions) )
rownames(rf.performances) <- c("F1 Score", "AUC", "Precision", "Recall",
                               "FPR", "FNR", "True Positives", "False Negatives", 
                               "True Negatives", "False Positives")

importance.df.resamptest <- data.frame(gene = as.character( rownames(test.classifier$importance) ),
                                       importance_score = as.numeric( test.classifier$importance ) ) %>%
  arrange(desc(importance_score))

head(importance.df.resamptest, 10)


#note we need to remove the TdTom-transgene from the gene set

#roc curve
levels(predictions.test$engramobserved) <- c(0,1)
#Plotting ROC...
roc.engramcell <- roc(test.predictions$engramobserved, 
                      as.numeric(test.predictions$Pos) )

# there is an error here the predictions.Hoch5k.lognorm$engramobserved is showing only as 1 which cannot be true
# seomthing is wrong with the code don't know where this comes from

#roc.inactive <- roc(predictions$inactiveobserved, as.numeric(predictions$Fos_neg) )

dev.off()
jpeg("ROC_lognormChen2020.jpg", width = 350, height = "350")
plot(roc.engramcell, col = "red", main = "ROC of Chen RF Classifier")
dev.off()

