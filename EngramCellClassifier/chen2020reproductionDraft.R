library(Seurat)
library(SeuratData)
library(patchwork)
library(metap)

#----Seurat explanation of algorithm in v3 here:
# Stuart, T., Butler, A., Hoffman, P., Hafemeister, C., Papalexi, E., Mauck III,
# W. M., ... & Satija, R. (2019). Comprehensive integration of single-cell data. 
# Cell, 177(7), 1888-1902.
# https://www.cell.com/cell/fulltext/S0092-8674(19)30559-8?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867419305598%3Fshowall%3Dtrue

#paper claiming to have a better method of doing batch correction

InstallData("ifnb")


# load dataset
LoadData("ifnb")

# split the dataset into a list of two seurat objects (stim and CTRL)
ifnb.list <- SplitObject(ifnb, split.by = "stim")

# normalize and identify variable features for each dataset independently
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = ifnb.list)
immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)


# this command creates an 'integrated' data assay
immune.combined <- IntegrateData(anchorset = immune.anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)

# Visualization
p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "stim")
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2

DimPlot(immune.combined, reduction = "umap", split.by = "stim")

DefaultAssay(immune.combined) <- "RNA"
nk.markers <- FindConservedMarkers(immune.combined, ident.1 = 6, grouping.var = "stim", verbose = FALSE)
head(nk.markers)

FeaturePlot(immune.combined, features = c("CD3D", "SELL", "CREM", "CD8A", "GNLY", "CD79A", "FCGR3A","CCL2", "PPBP"),
            min.cutoff = "q9")

#making name vector for numeric idents
num.names = c(0:14)
x <- c()
for(i in num.names){
  x <- c(x,as.character(num.names[i]))
}
num.names <- x
rm(x)

#marker names
mrkr.names <- c("CD14 Mono",  "CD4 Naive T", "CD4 Memory T", "CD16 Mono", "B", 
                "CD8 T", "NK", "T activated", "DC", "B Activated", "Mk", 
                "pDC", "Eryth", "Mono/Mk Doublets")

immune.combined <- RenameIdents(immune.combined, '0' = "CD14 Mono", '1' = "CD4 Naive T", '2' = "CD4 Memory T",
                                '3' = "CD16 Mono", '4' = "B", '5' = "CD8 T", '6' = "NK", '7' = "T activated", '8' = "DC",
                                '9' = "B Activated",'10' = "Mk", '11' = "pDC", '12' = "Eryth", '13' = "Mono/Mk Doublets")

immune.combined <- RenameIdents(immune.combined, 
                                old.ident.names = num.names,
                                new.ident.names = mrkr.names)
DimPlot(immune.combined, label = TRUE)



######################
#SCTtransform normalization integrated analysis

ifnb.list <- SplitObject(ifnb, split.by = "stim")
ifnb.list <- lapply(X = ifnb.list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = 3000)
ifnb.list <- PrepSCTIntegration(object.list = ifnb.list, anchor.features = features)

immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, normalization.method = "SCT",
                                         anchor.features = features)
immune.combined.sct <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT")

immune.combined.sct <- RunPCA(immune.combined.sct, verbose = FALSE)
immune.combined.sct <- RunUMAP(immune.combined.sct, reduction = "pca", dims = 1:30)

p1 <- DimPlot(immune.combined.sct, reduction = "umap", group.by = "stim")
p2 <- DimPlot(immune.combined.sct, reduction = "umap", group.by = "seurat_annotations", label = TRUE,
              repel = TRUE)
p1 + p2

#taken from above script
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)


#Attempt with Chen data

chen2020 <- CreateSeuratObject(counts = chen2020_counts,
                               meta.data = chen2020_meta)

chen2020.list <- SplitObject(chen2020, split.by = "condition_label")
chen2020.list <- lapply(X = chen2020.list, FUN = SCTransform)

# chen2020.list <- chen2020.list[c(1:7,9:20)] this was for when I was trying to correct
# for batch by mouse
features <- SelectIntegrationFeatures(object.list = chen2020.list, nfeatures = 3000)
# not clear how to pick nfeatures
#error is produced, invalid x


chen2020.list <- PrepSCTIntegration(object.list = chen2020.list, anchor.features = features)

chen2020.anchors <- FindIntegrationAnchors(object.list = chen2020.list, normalization.method = "SCT",
                                         anchor.features = features)
chen2020.combined.sct <- IntegrateData(anchorset = chen2020.anchors, normalization.method = "SCT")

chen2020.combined.sct <- RunPCA(chen2020.combined.sct, verbose = FALSE)
chen2020.combined.sct <- RunUMAP(chen2020.combined.sct, reduction = "pca", dims = 1:30)
chen2020.combined.sct <- FindNeighbors(chen2020.combined.sct, reduction = "pca", dims = 1:30)
chen2020.combined.sct <- FindClusters(chen2020.combined.sct, resolution = 0.5)

p1 <- DimPlot(chen2020.combined.sct, reduction = "umap", group.by = "condition_label")
p2 <- DimPlot(chen2020.combined.sct, reduction = "umap", group.by = "engram_label")
p3 <- DimPlot(chen2020.combined.sct, reduction = "umap", group.by = "seurat_clusters")
p1 + p2 + p3

#exploring cell clusters
DefaultAssay(chen2020.combined.sct) <- "RNA"
nk.markers <- FindConservedMarkers(chen2020.combined.sct, ident.1 = 6, grouping.var = "seurat_clusters", verbose = FALSE)
head(nk.markers)

ct.markers <- FindAllMarkers(object = chen2020.combined.sct)
setwd("C:/Users/angus/Desktop/PavLabEngrams/EngramCellClassifier")
write(data.frame(ct.markers), "chen2020clusters_condtion_label_as_batch_stcnorm_UMAP.csv")











