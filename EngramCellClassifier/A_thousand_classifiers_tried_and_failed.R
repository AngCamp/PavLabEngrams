# A thousand classifiers tried and failed

# scripts that didn't work, this is not just for classifiers but all
# normalisation methods etc etc stuff that did not work



# deal with this later
# hochgernerDGC_counts <- hochgerner5k_2018_counts[, hochgerner5k_2018_meta$cluster_name == "Granule-mature"]
# hochgernerDGC_counts$Symbol_mm <- rownames(hochgernerDGC_counts)
# hochgernerDGC_counts <- left_join(x = hg_to_mm, y = hochgernerDGC_counts, by = "Symbol_mm" )
# rownames(hochgernerDGC_counts) <- hochgernerDGC_counts$Symbol_hg
# hochgernerDGC_counts <- hochgernerDGC_counts[,c(4:dim(hochgernerDGC_counts)[2])]
# 

#cannonical correlation analysis, between hochgerner and ayhan
library(CCA) # provides cc() functions, kinda broken
library(vegan) # also performs cca
library(yacca)

#once we get the embeddings out we need to do a change of basis like so...
# V are the basis vectors, our componenets, D is some data we want to transform
# P = V*(-1)DV this gives us our projection. 
# For mouse data we can use X  components and for human we can use Y componenets
# from our CCA

# make transformations for mouse and human data
# explanation on hw to use the package here: https://www2.karlin.mff.cuni.cz/~maciak/NMST539/cvicenie11.html

hoch.lognorm <- log.norm(hochgernerDGC_counts) #made from hg orthologs list
ayhanDGC.lognorm <- log.norm(ayhanDGC_counts)

X <-  hoch.lognorm[ sample(rownames(hoch.lognorm),1000),]
Y <-  ayhanDGC.lognorm[ sample(rownames(ayhanDGC.lognorm),1000),]

X <- hoch.lognorm[, features]
Y <- ayhanDGC.lognorm[ , features ]
# from yacca package https://cran.r-project.org/web/packages/yacca/yacca.pdf
# our data is already centered,
# consider restricting to the top 2000 features from seurat
hoch_vs_ayhan.cca <- cca(X, Y, use.eigs = TRUE, max.dim = 100,
                         xcenter = FALSE, ycenter = FALSE)

hoch_vs_ayhan.cca <- cca(X, Y, use.eigs = TRUE, xcenter = FALSE, ycenter = FALSE,
                         reg.param = 1) # this reparameter will make the calculation more stable
# but give outliers and stronger effect, this may be desirable for us astcually


# attempt 2 with less cells and features
X <- hoch.lognorm[ sample( rownames(hoch.lognorm) , 1000) , features]
Y <- ayhanDGC.lognorm[ sample( rownames(ayhanDGC.lognorm) , 1000) , features ]

cca.mouse_human_DGCs <- cc(X, Y) # running this like this breaks cc 
# possibly used this: https://cran.r-project.org/web/packages/seedCCA/index.html

# try it with seurat

hoch <- CreateSeuratObject(hochgernerDGC_counts)

# I think hoch passes quality check as is
hoch[["percent.mt"]] <- PercentageFeatureSet(hoch, pattern = "^MT-")
VlnPlot(hoch, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

ayhan <- CreateSeuratObject(ayhanDGC_counts)
# Ayhan dataa looks filtered as well anything with less than 300 genes was removed
#


hoch <- NormalizeData(hoch, normalization.method = "LogNormalize", scale.factor = 10000)
ayhan <- NormalizeData(ayhan, normalization.method = "LogNormalize", scale.factor = 10000)


hoch <- FindVariableFeatures(hoch, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10.hochDGC <- head(VariableFeatures(hoch), 10)


# plot variable features with and without labels
plot1 <- VariableFeaturePlot(hoch)
plot1 <- LabelPoints(plot = plot1, points = top10.hochDGC, repel = FALSE)

plot1

all.genes.human <- rownames(hoch)
hoch  <- ScaleData(hoch, features = all.genes.human)
hoch <- RunPCA(hoch, features = VariableFeatures(object = hoch))
VizDimLoadings(hoch, dims = 1:2, reduction = "pca")

write_csv(ayhanDGC.lognorm, "humanlognorm.csv")
# humans, gonna do the clustering just to take a peak at its results
ayhan <- FindVariableFeatures(ayhan, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10.humanDGC <- head(VariableFeatures(ayhan), 10)

# plot variable features with and without labels
plot2 <- VariableFeaturePlot(ayhan)
plot2 <- LabelPoints(plot = plot2, points = top10.humanDGC, repel = TRUE)
plot2

all.genes.human <- rownames(ayhan)
ayhan  <- ScaleData(ayhan, features = all.genes.human)

# PCA results
ayhan <- RunPCA(ayhan, features = VariableFeatures(object = ayhan))

# Examine and visualize PCA results a few different ways
print(ayhan[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(hoch, dims = 1:2, reduction = "pca")

# ayhan [["percent.mt"]] <- PercentageFeatureSet(ayhan, pattern = "^MT-")
# VlnPlot(ayhan , features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# 
# plot1 <- FeatureScatter(ayhan, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(ayhan, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# plot1 + plot2

nVarFeats = 8000
hoch <- CreateSeuratObject(hochgernerDGC_counts)
#hoch <- FindVariableFeatures(hoch, selection.method = "vst", nfeatures = nVarFeats)
hoch <- NormalizeData(hoch, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes.human <- rownames(hoch)
hoch  <- ScaleData(hoch, features = all.genes.human)


ayhan <- CreateSeuratObject(ayhanDGC_counts)
ayhan <- NormalizeData(ayhan, normalization.method = "LogNormalize", scale.factor = 10000)all.genes.human <- rownames(hoch)
ayhan  <- ScaleData(ayhan, features = all.genes.human)
#ayhan <- FindVariableFeatures(ayhan, selection.method = "vst", nfeatures = nVarFeats)

# seurat gives us no way to determine significance etc
merged.cca <- RunCCA(hoch, ayhan, num.cc = 50, features = all.genes.human)
mouse.loadings.cca <- merged.cca@reductions$cca@feature.loadings

VizDimLoadings(merged.cca, dims = 1:2, reduction = "cca")

jeager <- CreateSeuratObject(combined.counts.hg[rownames(mouse.loadings.cca),])
jeager <- NormalizeData(jeager, normalization.method = "LogNormalize", scale.factor = 10000)
jeager <- ScaleData(jeager, features = all.genes.human)

#getting jeager data back out
# or conversely if we want to work with sparse matrices https://www.tutorialspoint.com/how-to-create-a-sparse-matrix-in-r
jeager.seuratnorm <- as.matrix(GetAssayData(jeager, assay="RNA", slot="data"))

# embed jeager inside the mouse coordinates
# make sure you know what these are, verify what a loading is mathematically
# may not be what 
# Relevant pub: https://scholarscompass.vcu.edu/cgi/viewcontent.cgi?article=1001&context=socialwork_pubs#:~:text=Canonical%20loadings%2C%20also%20called%20structure,and%20that%20set's%20canonical%20variate.
# This may be relevant as well:
# https://arxiv.org/pdf/1907.01693.pdf

jeager.embedded <- t(mouse.loadings.cca) %*% jeager.seuratnorm
jeager.embedded <- t( jeager.embedded )

# This neuropaper seems to have done waht you're doing here.
# Drysdale, A. T., Grosenick, L., Downar, J., Dunlop, K., Mansouri, F., Meng, Y., 
# ... & Liston, C. (2017). Resting-state connectivity biomarkers define 
# neurophysiological subtypes of depression. Nature medicine, 23(1), 28-38.
# https://www.nature.com/articles/nm.4246?TB_iframe=true&width=921.6&height=921.6

#
fospos_cca <- jeager.embedded[combined.meta$fos_status=="Fos+",]
fosneg_cca <- jeager.embedded[combined.meta$fos_status=="Fos-",]

test <- scale(jeager.embedded)

fospos_cca <- test[combined.meta$fos_status=="Fos+",]
fosneg_cca <- test[combined.meta$fos_status=="Fos-",]


test.identitymatrix <- ginv(mouse.coordindates.cca) %*% mouse.coordindates.cca
# we need the mass package to calculate the inverse as the coordinates here are a 
# retangular matrix
# https://stackoverflow.com/questions/21364060/calculate-inverse-of-a-non-square-matrix-in-r
# we can use rectangular matrices to change basis: https://math.stackexchange.com/questions/2928327/change-of-basis-of-a-linear-map-defined-by-non-square-matrix


# writting RF in caret

# it may just be easier to use seurats normalization on jeager data
# and extract the normalized data

# project jeager into CCA space
# combined.counts.hg
