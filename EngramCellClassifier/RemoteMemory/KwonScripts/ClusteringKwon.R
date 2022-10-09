# Clustering Kwon in seurat
# I could not reproduce the clusters from Chen in the Kwon data


library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(stringr)
library(purrr)

# ## Libraries
# library(randomForest)
# library(rfUtilities)
# library(stringr)
# library(sampler)
# library(caTools)
# library(pROC)
# library(stats)
# library(Dict)
# library(pheatmap)
# library(caret)
# library(data.table)

# OG Paper
# Kwon, D. Y., Xu, B., Hu, P., Zhao, Y. T., Beagan, J. A., Nofziger, J. H., ... &
# Zhou, Z. (2022). Neuronal Yin Yang1 in the prefrontal cortex regulates 
# transcriptional and behavioral responses to chronic stress in mice. Nature 
# communications, 13(1), 1-19.
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8748737/

setwd("~/test_datasets/Kwon2021_GSE145970")

kwon2021_counts <- read.csv("kwon2021_dropseq_counts.csv", header =T) # generated in cleaningkwon2021.r
genes <- kwon2021_counts$GENE # save the gene names
kwon2021_counts <- kwon2021_counts[,2:dim(kwon2021_counts)[2]]
kwon2021_counts <- apply(kwon2021_counts, 2, FUN=as.numeric) # drops the rownames
kwon2021_counts <- data.frame(kwon2021_counts)
kwon2021_counts[is.na(kwon2021_counts)] <- 0
rownames(kwon2021_counts) <- genes # put the genes back in
#> sum(is.na(kwon2021_counts)) [1] 76300000
#> > 76300000/1141360000
# [1] 0.06685007
# about 6% of the data is na's


#kwon2021_meta <- read.table("SraRunTable.txt")
# generated in ScrubletKwon2021.py
doubletscore <- read.csv("kwon2021_doubtlescores_scrublet.csv")
rownames(doubletscore) <- colnames(kwon2021_counts)

#Normal Seurat Workflow
#from here: https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
kwon2021 <- CreateSeuratObject(counts = kwon2021_counts,  meta.data = doubletscore,
                               min.cells = 10)


kwon2021[["percent.mt"]] <- PercentageFeatureSet(kwon2021, pattern = "^MT-")
#checks for genes beginning in MRP or RP, https://www.genenames.org/data/genegroup/#!/group/1054


#selecting variable genes and normalizing
kwon2021 <- subset(kwon2021, subset = nFeature_RNA > 600 & nFeature_RNA < 5000 & percent.mt < 5 & doublet_score < 0.2)
kwon2021
# by the end of filtering we are left with...
# An object of class Seurat 
# 25887 features across 30615 samples within 1 assay 
# Active assay: RNA (25887 features, 0 variable features)

# compared to their results: "As a result, 31,806 nuclei from 8 samples 
# (4 control and 4 CUS) were kept for downstream analysis"

kwon2021 <- NormalizeData(kwon2021, normalization.method = "LogNormalize", scale.factor = 10000)
kwon2021 <- FindVariableFeatures(kwon2021, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(kwon2021), 10)
top10

dev.off()
jpeg("varfeatkwon2021.jpg", width = 350, height = 350)
VariableFeaturePlot(kwon2021)
dev.off()

#scaling
all.genes <- rownames(kwon2021)
kwon2021 <- ScaleData(kwon2021, features = all.genes)
kwon2021 <- RunPCA(kwon2021, 
                   features = VariableFeatures(object = kwon2021),
                   npcs = 50) #as per the text the top 50 were selected
# looks like they didn't bother with the jackstraw plots
# 50 is the default https://satijalab.org/seurat/reference/runpca

dev.off()
jpeg("pcloadingskwon2021.jpg", width = 350, height = 350)
VizDimLoadings(kwon2021, dims = 1:2, reduction = "pca")
dev.off()


dev.off()
jpeg("DimPlotPCA_kwon2021.jpg", width = 350, height = 350)
DimPlot(kwon2021, reduction = "pca")
dev.off()


dev.off()
jpeg("HeatmapPCA_kwon2021.jpg", width = 350, height = 350)
DimHeatmap(kwon2021, dims = 1, cells = 500, balanced = TRUE)
dev.off()

kwon2021 <- JackStraw(kwon2021, num.replicate = 100)
kwon2021 <- ScoreJackStraw(kwon2021, dims = 1:20)

dev.off()
jpeg("ElbowPlotPCA_kwon2021.jpg", width = 350, height = 350)
JackStrawPlot(kwon2021, ndims = 1:20)
dev.off()

dev.off()
jpeg("ElbowPlotPCA_kwon2021.jpg", width = 350, height = 350)
ElbowPlot(kwon2021)
dev.off()

kwon2021 <- RunUMAP(kwon2021, dims=1:50, dim.embed = 2)
#clustering usign the top 20 PCs as they did in kwon2021
kwon2021 <- FindNeighbors(kwon2021, reduction = "pca", dims = 1:50)
kwon2021 <- FindClusters(kwon2021, resolution = 1)
# they tried this between 0.2 and 3 and found 38 clusters
# I find 36, they also note an extreme batch effect
# I feel at this point clustering may be best done using the integration tutorial

#Plotting the clusters
dev.off()
jpeg("DimplotTSNE_kwon2021.jpg", width = 700, height = 700)
DimPlot(kwon2021, reduction = "tsne")
dev.off()



#Idendtifying the big ones, Neurons, astrocytes, glial cells
# DotPlot may be good for later when comparing between groups
# https://satijalab.org/seurat/reference/dotplot

#Finding Markers
kwon2021.markers <- FindAllMarkers(kwon2021, only.pos = TRUE)

write_csv(kwon2021.markers, "kwon_allcell_clustering_markers.csv")

clust <- data.frame(kwon2021@meta.data$orig.ident)
clust$seurat_all_cells_clusters <- kwon2021@meta.data$seurat_clusters

write_csv(clust, "kwon_allcells_clustering_markers.csv")



#checking for the markers described in Chen et al., 2020
test <- kwon2021.markers %>%
  group_by(cluster)

test <- data.frame(test)
test <- test[,test$p_val_adj<0.05]
layer.markers <- c("Dkkl1","Rprm","Calb2","Tesc","Tnfaip8l3","Tshz2","Lhx6")
layer.markers.idx <- which(test$gene %in% layer.markers)
test$cluster[layer.markers.idx]
# > test$cluster[layer.markers.idx]
# [1] 1  18 20 23 32 34 35 35
# > layer.markers[(layer.markers %in%  test$gene)]
# [1] "Rprm"      "Tnfaip8l3" "Tshz2"

# I've got to get it down to just the neurons

DimPlot(kwon2021, reduction = "umap",  label = TRUE)

kwon2021.markers$gene[kwon2021.markers$cluster==6]

layer.markers[(layer.markers %in%  test$gene)]

#they are there and match the cluster assignments from kwon2021
# > test$cluster[layer.markers.idx]
# [1] 0 1 2 3 4 5 6
# Levels: 0 1 2 3 4 5 6

test$p_val_adj[layer.markers.idx]
#The genes are significant even after bonferonni correction
# > test$p_val_adj[layer.markers.idx]
# [1]  0.000000e+00  0.000000e+00  0.000000e+00 6.459740e-244 1.979667e-140
# [6] 8.827732e-301  0.000000e+00


top2.kwon2021.markers <- kwon2021.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

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




# clusters we want are....
#kwon2021.markers$cluster[which(kwon2021.markers$gene %in% c("Snap25","Gad2","Slc17a7") )]
# 0  1 2  3  5 7  8  10 11 13 15 16 22 24 29 31 35
FeaturePlot(kwon2021, features = "Slc17a7")# glut transporter
FeaturePlot(kwon2021, features = "Gad2") # gaba transporter
# claustrum is cluster 21
# > sum(kwon2021@meta.data$seurat_clusters==21)
# [1] 432
# > 432/30615
# [1] 0.01411073 
FeaturePlot(kwon2021, features = "Nr4a2") # claustrum excitatory cells
FeaturePlot(kwon2021, features = "Gnb4") # claustrum excitatory cells
FeaturePlot(kwon2021, features = "Ppp1r1b",  min.cutoff = 1.5) #Striatal inhibitory neurons
# clusters 6 and 14 Lmo3 and Gad2 are also markers
# > sum(kwon2021@meta.data$seurat_clusters==6)
# [1] 1345
# > sum(kwon2021@meta.data$seurat_clusters==14)
# [1] 777
# > kwon2021
# An object of class Seurat 
# 25887 features across 30615 samples within 1 assay 
# Active assay: RNA (25887 features, 2000 variable features)
# 2 dimensional reductions calculated: pca, umap
# > (1345+777)/30615
# [1] 0.06931243 so about 7 percent
FeaturePlot(kwon2021, features = "Gja1") # astrocyte marker, 
FeaturePlot(kwon2021, features = "Pdgfra") # vasculature marker, 
FeaturePlot(kwon2021, features = "Ctss") # microglia marker, 
FeaturePlot(kwon2021, features = "Enpp6") #  Oligo2
FeaturePlot(kwon2021, features = "Flt1") # endothelial cells





chen_markers <- read_csv("~/PavLabEngrams/EngramCellClassifier/Chen2020_GSE152632/Chen2020_cellclusterlabels.csv")


test <- chen_markers %>%
  group_by(cluster)

test <- data.frame(test)
chen.layer.markers <- c("Dkkl1","Rprm","Calb2","Tesc","Tnfaip8l3","Tshz2","Lhx6")
layer.markers.idx <- which(test$gene %in% chen.layer.markers)
test$cluster[layer.markers.idx]

top2.kwon2021.markers <-kwon_markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
top2.kwon2021.markers

top2.kwon2021.markers.conserved <-conserved_markers.summary %>%
  group_by(cluster_id) %>%
  slice_max(n = 2, order_by = rowmean_log2FC)
top2.kwon2021.markers.conserved

FeaturePlot(kwon2021_mpfc.combined.sct, features = "Slc17a7")# glut transporter
FeaturePlot(kwon2021_mpfc.combined.sct, features = "Gad2")
DimPlot(kwon2021_mpfc.combined.sct)





###############################
########### Second Round of Clustering to get it as close to chen as possible
###############################

mpfc_neural_clusters.idx <- which(kwon2021@meta.data$seurat_clusters %in% c(0,1,2,3,5,7,8,10,11,13,15,16,22,24,29,31,35) )
keepthese <- kwon2021@meta.data$CellID[mpfc_neural_clusters.idx]
mpfc.idx <- doubletscore$CellID %in% keepthese
kwon2021_mpfc_neurons_counts <- kwon2021_counts[,mpfc.idx]
# due to an error in the python script the cell naems are slightly different
# between the colnames() of counts and the doublescore$CellID
kwon2021_mpfc_meta <- doubletscore[mpfc.idx,]
kwon2021_mpfc_meta$allcells_seurat_cluster <- kwon2021@meta.data$seurat_clusters[which(kwon2021@meta.data$seurat_clusters %in% c(0,1,2,3,5,7,8,10,11,13,15,16,22,24,29,31,35) )]

kwon2021_mpfc_meta$condition <- as.character(lapply(kwon2021_mpfc_meta$CellID, function(x) if (grepl("Ctrl", x, fixed=TRUE)) "Control" else {x}))
kwon2021_mpfc_meta$condition <- as.character(lapply(kwon2021_mpfc_meta$condition, function(x) if (grepl("Stress", x, fixed=TRUE)) "CUS" else {x}))
kwon2021_mpfc_meta$condition <- as.factor(kwon2021_mpfc_meta$condition)

kwon2021_mpfc_meta$mouse.replicate <- as.character(lapply(kwon2021_mpfc_meta$CellID, function(x) if (grepl("Rep30", x, fixed=TRUE)) "Rep30" else {x}))
kwon2021_mpfc_meta$mouse.replicate  <- as.character(lapply(kwon2021_mpfc_meta$mouse.replicate , function(x) if (grepl("Rep31", x, fixed=TRUE)) "Rep31" else {x}))
kwon2021_mpfc_meta$mouse.replicate  <- as.character(lapply(kwon2021_mpfc_meta$mouse.replicate , function(x) if (grepl("Rep33", x, fixed=TRUE)) "Rep33" else {x}))
kwon2021_mpfc_meta$mouse.replicate  <- as.character(lapply(kwon2021_mpfc_meta$mouse.replicate , function(x) if (grepl("Rep34", x, fixed=TRUE)) "Rep34" else {x}))
kwon2021_mpfc_meta$mouse.replicate  <- as.character(lapply(kwon2021_mpfc_meta$mouse.replicate , function(x) if (grepl("Rep45", x, fixed=TRUE)) "Rep45" else {x}))
kwon2021_mpfc_meta$mouse.replicate  <- as.character(lapply(kwon2021_mpfc_meta$mouse.replicate , function(x) if (grepl("Rep46", x, fixed=TRUE)) "Rep46" else {x}))
kwon2021_mpfc_meta$mouse.replicate  <- as.character(lapply(kwon2021_mpfc_meta$mouse.replicate , function(x) if (grepl("Rep48", x, fixed=TRUE)) "Rep48" else {x}))
kwon2021_mpfc_meta$mouse.replicate  <- as.character(lapply(kwon2021_mpfc_meta$mouse.replicate , function(x) if (grepl("Rep49", x, fixed=TRUE)) "Rep49" else {x}))
kwon2021_mpfc_meta$mouse.replicate  <- as.factor(kwon2021_mpfc_meta$mouse.replicate)

rownames(kwon2021_mpfc_meta) <- colnames(kwon2021_mpfc_neurons_counts)


# follow the integration tutorial from here and do it by replicates to remove batch effects then
# cluster
# In Kwon et al., (2021) they mention a considerable batch effect so
# our clustering should be able to accomodate this

kwon2021_mpfc <- CreateSeuratObject(counts = kwon2021_mpfc_neurons_counts, 
                                    meta.data = kwon2021_mpfc_meta)

# https://satijalab.org/seurat/articles/integration_introduction.html#performing-integration-on-datasets-normalized-with-sctransform-1
# in particular lets rely on the sctransform


kwon2021_mpfc.list <- SplitObject(kwon2021_mpfc, split.by = "mouse.replicate")
kwon2021_mpfc.list <- lapply(X = kwon2021_mpfc.list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = kwon2021_mpfc.list, nfeatures = 3000)
kwon2021_mpfc.list <- PrepSCTIntegration(object.list = kwon2021_mpfc.list, 
                                         anchor.features = features)

kwon2021_mpfc.anchors <- FindIntegrationAnchors(object.list = kwon2021_mpfc.list,
                                                normalization.method = "SCT",
                                                anchor.features = features)
kwon2021_mpfc.combined.sct <- IntegrateData(anchorset = kwon2021_mpfc.anchors, 
                                            normalization.method = "SCT")

kwon2021_mpfc.combined.sct <- RunPCA(kwon2021_mpfc.combined.sct, verbose = FALSE)
# jackplots dont work on sct transformed data

kwon2021_mpfc.combined.sct <- RunUMAP(kwon2021_mpfc.combined.sct, 
                                      reduction = "pca", 
                                      dims = 1:30)

kwon2021_mpfc.combined.sct <- FindNeighbors(kwon2021_mpfc.combined.sct,
                                            reduction = "pca", dims = 1:30)

kwon2021_mpfc.combined.sct <- FindClusters(kwon2021_mpfc.combined.sct,
                                           resolution = 0.03)
# we get 7 communities out of this with the same resolution as chen
# possibly the same cell types

# we need to rewrite find conserved markers slightly to iterate it across clusters
# https://hbctraining.github.io/scRNA-seq/lessons/sc_exercises_integ_marker_identification.html

#https://www.biostars.org/p/409790/

# to use FindConservedMarkers
# install.packages('BiocManager')
# BiocManager::install('multtest')
# install.packages('metap') this one is giving an issue, non zero exit because of qqconf which wont install
# potentiall solution: https://github.com/satijalab/seurat/issues/5957
# apt-get install libfftw3-dev

# Trying while correcting for batch effects
# we want clusters that are conserved across groups
get_conserved <- function(cluster){
  FindConservedMarkers(kwon2021_mpfc.combined.sct,
                       ident.1 = cluster,
                       grouping.var = "mouse.replicate",
                       only.pos = TRUE) %>%
    rownames_to_column(var = "gene") %>%
    cbind(cluster_id = cluster, .)
}

clusters.pointzerothree <- unique( kwon2021_mpfc.combined.sct@meta.data$integrated_snn_res.0.03 )

DefaultAssay(kwon2021_mpfc.combined.sct) <- "RNA"
conserved_markers <- map_dfr(clusters.pointzerothree , get_conserved)
keep.idx <- which(colnames(conserved_markers) %in% c("cluster_id", "gene", "max_pval", "minimump_p_val") )
n <- dim(conserved_markers)[2]
conserved_markers.summary <- conserved_markers[,c(1,2,n-1,n)]

avg_log2FC.idx  <- as.logical(lapply(colnames(conserved_markers), function(x) if (grepl("avg_log2FC", x, fixed=TRUE)) TRUE else FALSE))
conserved_markers.summary$rowmean_log2FC <- rowMeans(conserved_markers[,avg_log2FC.idx])

#  Attempting to use regular marker identification
kwon_markers <- FindAllMarkers(kwon2021_mpfc.combined.sct)

test <- kwon_markers %>%
  group_by(cluster)

test <- data.frame(test)
# chen.layer.markers <- c("Dkkl1","Rprm","Calb2","Tesc","Tnfaip8l3","Tshz2","Lhx6")
# layer.markers.idx <- which(test$gene %in% chen.layer.markers)
# layer.markers.idx <- which(test$gene %in% "Lhx6")
# test$cluster[layer.markers.idx]
# test$avg_log2FC[layer.markers.idx]
# # trying without correcting for batch effects

# In the end this was fruitless, the chen cluster markers are simply not revealing
# themselves in clusters or appearing in many.

top.kwon.markers <- kwon_markers %>%
  group_by(cluster) %>%
  slice_max(n = 1, order_by = avg_log2FC)

kwon2021_mpfc_meta$mpfc_integrated_seurat_cluster <- kwon2021_mpfc.combined.sct@meta.data$seurat_clusters
kwon2021_mpfc_meta$mpfc_integrated_marker <- as.character(lapply(as.numeric(kwon2021_mpfc_meta$mpfc_integrated_seurat_cluster), function(x) top.kwon.markers$gene[x]))

write_csv( kwon2021_mpfc_meta, "~/test_datasets/Kwon2021_GSE145970/kwon_mpfc_neurons_meta.csv")


