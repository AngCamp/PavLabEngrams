# Clustering Kwon in seurat
# I could not reproduce the clusters from Chen in the Kwon data


library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(stringr)
library(purrr)
library(R.utils)

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

kwon2021_counts <- read.table("~/test_datasets/Kwon2021_GSE145970/kwon2021_dropseq_counts.csv.gz", 
                            header =T) # generated in cleaningkwon2021.r
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
kwon2021 <- FindClusters(kwon2021, resolution = 0.46) # we get 26 at this setting
# from the text:
# "Clusters were identified using the function FindCluster in Seurat with the
# resolution parameter set to 1. Cells were classified into 21–45 clusters with 
# the resolution parameter from 0.3 to 2. Clustering resolution parameters were 
# varied quantitatively based on the number of cells being clustered. After the 
# clustering results with different resolutions were compared and evaluated, we 
# chose a resolution value of 1. Using this approach we were able to assign 31,806 
# cells to 28 clusters"
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
kwon2021.markers.res0.46 <- FindAllMarkers(kwon2021, only.pos = TRUE)

top.kwon2021.markers.res0.46 <- kwon2021.markers %>%
  group_by(cluster) %>%
  slice_max(n = 1, order_by = avg_log2FC) %>% 
  data.frame()
top.kwon2021.markers.res0.46

# because we had two Tshz2 (like they did I )
top2.kwon2021.markers.res0.46 <- kwon2021.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC) %>% 
  data.frame()


kwon2021@meta.data$markers_res0.046 <- as.character(lapply(kwon2021@meta.data$RNA_snn_res.0.46, function(x) top.kwon2021.markers.res0.046$gene[as.numeric(x)]))

# rename cluster id's
clust.ident <- c("Nectin3-0", "Rorb-1","Il1rapl2-2", "Ndst4-3", "Zfpm2-4", "Rgs9-5",
                "Sst-6", "Pbx3-7", "Slc1a3-8", "Adarb2-9", "Fam19a1-10", "Plp1-11",
                "Cdh18-12", "Ntng1-13", "Tshz2-Neurod6-14", "Pdgfra-15", "Tshz2-Vwc2l-16", "Nr4a2-17",
                "Hexb-18", "Gpc5-19", "Nefm-20", "Ptgds-21", "Zbtb20-22", "Rgs5-23", 
                "9630013A20Rik-24", "Zfp804b-25")

names(clust.ident) <- levels(kwon2021)
kwon2021 <- RenameIdents(kwon2021, clust.ident)

FeaturePlot(kwon2021, label = TRUE, repel = TRUE, features = "Angpt1") 

#kwon2021.markers.roc <- FindAllMarkers(kwon2021, test.use = "roc", only.pos = TRUE)
# 
# top.kwon2021.markers.by_roc <- kwon2021.markers.roc %>%
#   group_by(cluster) %>%
#   filter(avg_log2FC>2) %>%
#   slice_max(n = 5, order_by = myAUC)
# top.kwon2021.markers.by_roc

write_csv(kwon2021.markers.res0.46, "kwon_allcells_clustering_markers.csv")

clust <- data.frame(kwon2021@meta.data$orig.ident)
clust$seurat_all_cells_clusters <- kwon2021@meta.data$seurat_clusters

write_csv(kwon2021@meta.data, "~/test_datasets/Kwon2021_GSE145970/kwon2021_clustering_meta.csv")
kwon2021 <- RenameIdents(kwon2021, top.kwon2021.markers)

kwon2021.markers <- read_csv("~/test_datasets/Kwon2021_GSE145970/kwon_allcells_clustering_markers.csv")


cluster_info <- data.frame(
  cluster = as.factor(c(0:25)), 
  marker = names(table(kwon2021@active.ident)),
  num_cells = as.numeric(table(kwon2021@active.ident))
)
cluster_info$percent <- cluster_info$num_cells/sum(cluster_info$num_cells)


#############################################################
## Finding mpfc c36ell types and corresponding clusters 
###########################################################

# Clusters we dont want
# 5, 7, 8, 11, 15, 17, 18, 19, 21, 22, 23, 24


FeaturePlot(kwon2021, label = TRUE, features = "Ppp1r1b")
FeaturePlot(kwon2021, label = TRUE, features = "Meis2")
# Rgs9 - 5 Striatal cells based on Ppp1r1b and Meis2 expression

# Plp-11 remains unidentified but its clusetering with oligodendrocytes
# does not express neuronal markers

FeaturePlot(kwon2021, label = TRUE, repel = TRUE, features = "Pdgfra")
# OPC oligodendrocyte precursors, Pdgfra-15 

FeaturePlot(kwon2021, label = TRUE, features = "Gnb4") # claustrum excitatory cells
# Nr4a2 17


FeaturePlot(kwon2021, label = TRUE, features = "Ctss") # microglia marker,
# Hexb 18

FeaturePlot(kwon2021, label = TRUE, repel =TRUE, features = "Rorb")
FeaturePlot(kwon2021, label = TRUE, repel =TRUE, features = "Apoe")
FeaturePlot(kwon2021, label = TRUE, features = "Gja1") # astrocyte marker,
# cluster Slc1a3, Gpc5,  8, 19 and 

FeaturePlot(kwon2021, label = TRUE, features = "Mgp")
FeaturePlot(kwon2021, label = TRUE, features = "Fn1")
# Ptgds-21

FeaturePlot(kwon2021, label = TRUE, features = "Flt1") # endothelial cells
# Rgs5 23

FeaturePlot(kwon2021, label = TRUE, features = "Enpp6") #  Oligo2
# 9630013A20Rik-24

FeaturePlot(kwon2021, label = TRUE, repel =TRUE, features = "Prox1") 
FeaturePlot(kwon2021, label = TRUE, repel =TRUE, features = "Zfpm2") 
# based on the small number of cells and the Prox1 and Zfpm2 expression
# these are DGCs 

FeaturePlot(kwon2021, label = TRUE, repel =TRUE, features = "Meis2")
# Rgs9 - 5, Pbx3 - 7 I suspect are Meis2 Inhibitory neurons which are white matter 
# virtually absent in Chen shoud be exculded


FeaturePlot(kwon2021, label = TRUE, repel =TRUE, features = "Rorb")
#Gpc15 19 Astrocytes

# Clusters we want....
# 0, 1, 2, 3, 4, 6, 9, 10, 12, 13, 14, 16, 20, 25
 
FeaturePlot(kwon2021, label = TRUE, repel = TRUE, features = "Enpp2")
FeaturePlot(kwon2021, label = TRUE, repel = TRUE, features = "Cux2")
# Necti3 - 0 Cux2 and Enpp2 are in this so its likely ExL

# Rorb 1, Ex_Layer_4 

FeaturePlot(kwon2021, label = TRUE, repel = TRUE, features = "Ndst4")
# 3 Ndst4 Layer2/3_Ex

FeaturePlot(kwon2021, label = TRUE, features = "Pvalb") # Inh_Pv, Inh_Sst, Inh_Ndnf
# cluster 6 Sst-6

FeaturePlot(kwon2021, label = TRUE, features = "Vip") # Inh_Vip
# Adarb2-9

#Ex_Layer5
FeaturePlot(kwon2021, label = TRUE, repel = TRUE, features = "Cpne7")
FeaturePlot(kwon2021, label = TRUE, repel = TRUE, features = "Tox")
FeaturePlot(kwon2021, label = TRUE, repel = TRUE, features = "Tshz2")
# based on tox but no Tshz2 its Layer 5 and Cpne7 Fam19a1-10

# these are Ex_L6 cells
FeaturePlot(kwon2021, label = TRUE, repel = TRUE, features = "Zfpm2")
FeaturePlot(kwon2021, label = TRUE, repel = TRUE, features = "Hs3st4")
# Zfpm2  4, Cdh18  12

# 13 Ntng1 Ex_Ntng1

# Tshz2-Vxc2l 16 one of their Tshz2 clusters

# Nefm-20 no

FeaturePlot(kwon2021, label = TRUE, repel = TRUE, features = "Dkkl1")
FeaturePlot(kwon2021, label = TRUE, repel = TRUE, features = "Rgs4")
# they express a chen marker for glutamatergic cells,
# and a marker I found 




###############################
########### Making the meta data
###############################

cortex.clusters <- c(0, 1, 2, 3, 4, 6, 9, 10, 12, 13, 14, 16, 20, 25)

mpfc_neural_clusters.idx <- which(kwon2021@meta.data$seurat_clusters %in% cortex.clusters )

# due to an error in the python script the cell naems are slightly different
# between the colnames() of counts and the doublescore$CellID
kwon2021_mpfc_neurons_meta <- doubletscore[mpfc.idx,]
kwon2021_mpfc_neurons_meta$allcells_seurat_cluster <- kwon2021@meta.data$seurat_clusters[mpfc_neural_clusters.idx ]


clust.ident <- c("Nectin3", "Rorb","Il1rapl2", "Ndst4", "Zfpm2", "Rgs9",
                 "Sst", "Pbx3", "Slc1a3", "Adarb2", "Fam19a1", "Plp1",
                 "Cdh18", "Ntng1", "Tshz2-Neurod6", "Pdgfra", "Tshz2-Vwc2l", "Nr4a2",
                 "Hexb", "Gpc5", "Nefm", "Ptgds", "Zbtb20", "Rgs5", 
                 "9630013A20Rik", "Zfp804b")

kwon2021_mpfc_neurons_meta$cluster_markers_res0.46 <- as.character(lapply(kwon2021_mpfc_neurons_meta$allcells_seurat_cluster , function(x) clust.ident[as.numeric(x)]))


kwon2021_mpfc_neurons_meta$condition <- as.character(lapply(kwon2021_mpfc_neurons_meta$CellID, function(x) if (grepl("Ctrl", x, fixed=TRUE)) "Control" else {x}))
kwon2021_mpfc_neurons_meta$condition <- as.character(lapply(kwon2021_mpfc_neurons_meta$condition, function(x) if (grepl("Stress", x, fixed=TRUE)) "CUS" else {x}))
kwon2021_mpfc_neurons_meta$condition <- as.factor(kwon2021_mpfc_neurons_meta$condition)

kwon2021_mpfc_neurons_meta$mouse.replicate <- as.character(lapply(kwon2021_mpfc_neurons_meta$CellID, function(x) if (grepl("Rep30", x, fixed=TRUE)) "Rep30" else {x}))
kwon2021_mpfc_neurons_meta$mouse.replicate  <- as.character(lapply(kwon2021_mpfc_neurons_meta$mouse.replicate , function(x) if (grepl("Rep31", x, fixed=TRUE)) "Rep31" else {x}))
kwon2021_mpfc_neurons_meta$mouse.replicate  <- as.character(lapply(kwon2021_mpfc_neurons_meta$mouse.replicate , function(x) if (grepl("Rep33", x, fixed=TRUE)) "Rep33" else {x}))
kwon2021_mpfc_neurons_meta$mouse.replicate  <- as.character(lapply(kwon2021_mpfc_neurons_meta$mouse.replicate , function(x) if (grepl("Rep34", x, fixed=TRUE)) "Rep34" else {x}))
kwon2021_mpfc_neurons_meta$mouse.replicate  <- as.character(lapply(kwon2021_mpfc_neurons_meta$mouse.replicate , function(x) if (grepl("Rep45", x, fixed=TRUE)) "Rep45" else {x}))
kwon2021_mpfc_neurons_meta$mouse.replicate  <- as.character(lapply(kwon2021_mpfc_neurons_meta$mouse.replicate , function(x) if (grepl("Rep46", x, fixed=TRUE)) "Rep46" else {x}))
kwon2021_mpfc_neurons_meta$mouse.replicate  <- as.character(lapply(kwon2021_mpfc_neurons_meta$mouse.replicate , function(x) if (grepl("Rep48", x, fixed=TRUE)) "Rep48" else {x}))
kwon2021_mpfc_neurons_meta$mouse.replicate  <- as.character(lapply(kwon2021_mpfc_neurons_meta$mouse.replicate , function(x) if (grepl("Rep49", x, fixed=TRUE)) "Rep49" else {x}))
kwon2021_mpfc_neurons_meta$mouse.replicate  <- as.factor(kwon2021_mpfc_neurons_meta$mouse.replicate)

# making neuron counts and 
keepthese <- kwon2021@meta.data$CellID[mpfc_neural_clusters.idx]
mpfc.idx <- doubletscore$CellID %in% keepthese
kwon2021_mpfc_neurons_counts <- kwon2021_counts[,mpfc.idx]
kwon2021_mpfc_neurons_counts$GENE <- rownames(kwon2021_mpfc_neurons_counts)
rownames(kwon2021_mpfc_neurons_meta) <- colnames(kwon2021_mpfc_neurons_counts)


write_csv(kwon2021_mpfc_neurons_meta, "~/test_datasets/Kwon2021_GSE145970/kwon2021_mpfc_neurons_meta.csv")
write_csv(kwon2021_mpfc_neurons_counts, "~/test_datasets/Kwon2021_GSE145970/kwon2021_mpfc_neurons_counts.csv")
gzip("~/test_datasets/Kwon2021_GSE145970/kwon2021_mpfc_neurons_counts.csv",
     destname="~/test_datasets/Kwon2021_GSE145970/kwon2021_mpfc_neurons_counts.csv.gz")


###############################
########### Second Round of Clustering to get it as close to chen as possible
###############################


# follow the integration tutorial from here and do it by replicates to remove batch effects then
# cluster
# In Kwon et al., (2021) they mention a considerable batch effect so
# our clustering should be able to accomodate this

kwon2021_mpfc <- CreateSeuratObject(counts = kwon2021_mpfc_neurons_counts, 
                                    meta.data = kwon2021_mpfc_neurons_meta)

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

kwon2021_mpfc_neurons_meta$mpfc_integrated_seurat_cluster <- kwon2021_mpfc.combined.sct@meta.data$seurat_clusters
kwon2021_mpfc_neurons_meta$mpfc_integrated_marker <- as.character(lapply(as.numeric(kwon2021_mpfc_neurons_meta$mpfc_integrated_seurat_cluster), function(x) top.kwon.markers$gene[x]))

write_csv( kwon2021_mpfc_neurons_meta, "~/test_datasets/Kwon2021_GSE145970/kwon_mpfc_neurons_meta.csv")


