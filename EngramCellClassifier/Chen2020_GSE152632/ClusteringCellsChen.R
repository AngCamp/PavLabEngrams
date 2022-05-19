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
chen <- FindVariableFeatures(chen, selection.method = "vst", nfeatures = 8000)

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
chen <- ScoreJackStraw(chen, dims = 1:50)

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
chen <- FindClusters(chen, resolution = 0.5)

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
top2.chen.markers <- chen.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
top2.chen.markers

chen.markers %>%
  group_by(cluster) %>%
  slice_max(n = 1, order_by = avg_log2FC)

chen.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()


dev.off()
jpeg("DimplotTSNE_Chen.jpg", width = 1400, height = 700)
DimPlot(chen, reduction = "tsne") + FeaturePlot(chen, reduction =  "tsne", feature = c("Snap25") )
dev.off()



