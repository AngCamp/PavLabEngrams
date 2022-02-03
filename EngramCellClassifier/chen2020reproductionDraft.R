library(Seurat)
library(SeuratData)
library(patchwork)
library(metap)

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

chen2020.list <- SplitObject(chen2020, split.by = "source_name")
chen2020.list <- lapply(X = chen2020.list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = chen2020.list, nfeatures = 3000)
chen2020.list <- PrepSCTIntegration(object.list = chen2020.list, anchor.features = features)

immune.anchors <- FindIntegrationAnchors(object.list = chen2020.list, normalization.method = "SCT",
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



