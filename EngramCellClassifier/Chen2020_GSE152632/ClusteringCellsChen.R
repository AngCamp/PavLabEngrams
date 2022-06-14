# 
# 
# 
# This is the cell clustering an cleaning for Chen et al., (2020)
# 
# Chen, M. B., Jiang, X., Quake, S. R., & Südhof, T. C. (2020).
# Persistent transcriptional programmes are associated with remote memory.
# Nature, 587(7834), 437-442.

#Path on pavlab server
#setwd("/home/acampbell/PavLabEngrams/EngramCellClassifier/Chen2020_GSE152632")

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

#Load the data and Metadata
chen2020_counts <- read.csv('GSE152632_GEO_mberkchen_TRAP2_counts.csv.gz', header = TRUE)
rownames(chen2020_counts) <- chen2020_counts$X
chen2020_counts <- chen2020_counts[,2:3531]
chen2020_meta <- read.csv( 'SraRunTable.txt', header = TRUE)

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
make.predictions.df <- function(classifier.object, test_df){
  #generate predictions for making classifier summary
  predictions <- as.data.frame(predict(classifier.object, test_df[,1:(length(test_df)-1)], type = "prob"))
  predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)] #1:2 for the number of classes
  predictions$observed <- test_df$Engramcell #this should be changed if you want to make this functions more modular
  colnames(predictions)[1:2] <- c("Neg","Pos")
  predictions$engramobserved <- ifelse(predictions$observed=="Pos", 1, 0)
  predictions$inactiveobserved <- ifelse(predictions$observed=="Neg", 1, 0)
  return(predictions)
}


assessment <- function(predictions.df){
  # returns a vector of assessments to be used to make dataframe summarizing classifiers performance
  # can be used to make df of all calssifiers trained in a single run
  TP <- sum((predictions.df$predict == "Pos")&(predictions.df$observed == "Pos"))
  TN <- sum((predictions.df$predict == "Neg")&(predictions.df$observed == "Neg"))
  FN <- sum((predictions.df$predict == "Neg")&(predictions.df$observed == "Pos"))
  FP <- sum((predictions.df$predict == "Pos")&(predictions.df$observed == "Neg"))
  
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




logplusone <- function(x){
  #used in lognorm to transform the data
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


resample.randomForest <-function( df.in, proportion,
                                  batches, trees){
  #this function resamples from our samples and retrains new models then combines them
  # this is too prevent over fitting on cells
  trees.per.batch <- as.integer(trees/batches)
  n.cells <- trunc( sum(df.in$Engramcell=="Neg")*proportion)
  batches <- c(1:batches)
  for( batch in batches){
    #ssampl comes from sampler library
    resample.set <- rbind(ssamp(df=df.in[df.in$Engramcell=="Neg",], n=n.cells,
                                strata=cluster, over=0),
                          ssamp(df=df.in[df.in$Engramcell=="Pos",], n=n.cells,
                                strata=cluster, over=0)
    )
    #ssamp throws strata column at from of data frame, this removes it
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


### RUNNING THE SCRIPTS

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
log.chen$cluster <- as.factor(chen@active.ident)
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
test.classifier <- resample.randomForest( df.in = training_set.test, proportion = 0.8, 
                                          batches = 20, trees = 1000)
toc()

test.predictions <- make.predictions.df(test.classifier, validation_set.test)
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


# Cell type specific roc curve



roc.engramcell <- roc(test.predictions$engramobserved, 
                      as.numeric(test.predictions$Pos) )

dev.off()
jpeg("ROC_lognormChen2020.jpg", width = 350, height = "350")
plot(roc.engramcell, col = "red", main = "ROC of Chen RF Classifier")
dev.off()

