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





############  Combining hochgerner and ayhan AND jeager/lacar
#  Attempt number 2, integration with seurat
# Relvant reading, good review
# 
# Shafer, M. E. (2019). Cross-species analysis of single-cell transcriptomic data.
# Frontiers in cell and developmental biology, 7, 175.
# https://www.frontiersin.org/articles/10.3389/fcell.2019.00175/full

# Affinati, A. H., Sabatini, P. V., True, C., Tomlinson, A. J., Kirigiti, M., 
# Lindsley, S. R., ... & Rupp, A. C. (2021). Cross-species analysis defines the 
# conservation of anatomically segregated VMH neuron populations. Elife, 10, e69065.
# https://elifesciences.org/articles/69065

#We need to load jeager, hochgerner, and ayna and check for which genes are present

#then we can put hochgerner and ayhan  together


ayhan2021_counts <- read.csv("~/test_datasets/Ayhan2021_GSE160189/GSE160189_Hippo_Counts.csv.gz")
rownames(ayhan2021_counts) <- ayhan2021_counts$gene
ayhan2021_counts[is.na(ayhan2021_counts)] <- 0
ayhan2021_counts <- ayhan2021_counts[, c(2:dim(ayhan2021_counts)[2])]


ayhan2021_meta <- read.csv("~/test_datasets/Ayhan2021_GSE160189/meta.tsv",
                           sep = '\t', header = TRUE)

# note as per Ayhan et al., 2021 we do not want Den.Gyr3 as it is mostly from a single subject
ayhanDGC.idx <- ayhan2021_meta$Cell[(ayhan2021_meta$Cluster == "Den.Gyr2")|(ayhan2021_meta$Cluster == "Den.Gyr1")]
ayhanDGC.idx <- colnames(ayhan2021_counts) %in% ayhanDGC.idx

ayhanDGC_counts <- ayhan2021_counts[,ayhanDGC.idx]

hg_to_mm <- read.table("hg_mm_1to1_ortho_genes_DIOPT-v8.tsv", sep = '\t', header = TRUE)
# just run left join onto the appropriate column for each dataset

# filter out genes not present in all datasets
present.orthologs.idx <- (hg_to_mm$Symbol_hg %in% rownames(ayhanDGC_counts))&(hg_to_mm$Symbol_mm %in% rownames(combined.counts) )
present.orthologs.idx <- present.orthologs.idx & (hg_to_mm$Symbol_mm %in% rownames(hochgerner5k_2018_counts) )

hg_to_mm <- hg_to_mm[present.orthologs.idx,] # filter for matches, ~ 14403 present between Ayhan and Jeager/Lacar



# filtering for orthologs present in all data sets
ayhanDGC_counts$Symbol_hg <- rownames(ayhanDGC_counts) # for bringing to 
ayhanDGC_counts <- left_join(x = hg_to_mm, y = ayhanDGC_counts, by = "Symbol_hg" )
rownames(ayhanDGC_counts) <- ayhanDGC_counts$Symbol_hg
ayhanDGC_counts <- ayhanDGC_counts[,c(4:dim(ayhanDGC_counts)[2])]


combined.counts.hg <- combined.counts
combined.counts.hg$Symbol_mm <- rownames(combined.counts.hg)
combined.counts.hg <- left_join(x = hg_to_mm, y = combined.counts.hg, by = "Symbol_mm" )
rownames(combined.counts.hg) <- combined.counts.hg$Symbol_hg
combined.counts.hg <- combined.counts.hg[,c(4:(dim(combined.counts.hg)[2]) )]


hochgernerDGC_counts <- hochgerner5k_2018_counts[, hochgerner5k_2018_meta$cluster_name == "Granule-mature"]
hochgernerDGC_counts <- data.frame( lapply(hochgernerDGC_counts,as.numeric) ) # this data is for some reason all strings
colnames(hochgernerDGC_counts) <- colnames(hochgerner5k_2018_counts)[hochgerner5k_2018_meta$cluster_name == "Granule-mature"]
hochgernerDGC_counts$Symbol_mm <- rownames(hochgerner5k_2018_counts)
hochgernerDGC_counts <- left_join(x = hg_to_mm, y = hochgernerDGC_counts, by = "Symbol_mm" )
rownames(hochgernerDGC_counts) <- hochgernerDGC_counts$Symbol_hg
hochgernerDGC_counts <- hochgernerDGC_counts[,c(4:dim(hochgernerDGC_counts)[2])]


# feature selection
# species label, we are going to integrate using hochgerner and 
all.cells <- cbind(ayhanDGC_counts, combined.counts.hg)
#all.cells <- cbind(all.cells, hochgernerDGC_counts)

# making meta data and label to split data by
species.idx <- rep("human", dim(ayhanDGC_counts)[2]) 
species.idx <- c(species.idx, rep("mouse", dim(combined.counts.hg)[2]+dim(hochgernerDGC_counts)[2]) )


experiment <- rep("ayhan2021", dim(ayhanDGC_counts)[2])
experiment <- c( experiment, rep("jeager2018", dim(combined.counts.hg)[2]) )
experiment <- c( experiment, rep("hochgerner2018", dim(hochgernerDGC_counts)[2]) )


activity <- rep("unlabelled", dim(ayhanDGC_counts)[2]) 
activity <- c(activity, combined.meta$ActivityStatus)
activity <- c(activity, rep("unlabelled", dim(hochgernerDGC_counts)[2]) )


all.cells.meta <- data.frame(experiment)
all.cells.meta$species <- species.idx
all.cells.meta$activity <- activity

rownames(all.cells.meta) <- colnames(all.cells)


# making a seurat object to normalize using their anchors,
# could potentially do clustering later in order to observe if cells are clustering together
integration_obj <- CreateSeuratObject(counts = all.cells, 
                                      min.cells = 0, 
                                      min.features = 0, 
                                      meta.data = all.cells.meta)

integration_obj@meta.data$species <- all.cells.meta$species

# split the dataset into a list of two seurat objects (by species)
DGC.list <- SplitObject(integration_obj, split.by = "species")

# normalize and identify variable features for each dataset independently
# may have to play with number of genes to include, so we will include more
# genes than seurat recommends because we may find them of use
# we will use all the genes 
DGC.list <- lapply(X = DGC.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = dim(hg_to_mm)[1])
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures( object.list = DGC.list, nfeatures = 8000 )















##########  In high iegs expressing cells.Rmd

# attempt to split data, before just using group data

# splitting data to retain samples from each, taking 25% of what ever is there 

# get reactivated data index
jeager.react.idx <- combined.meta$Reactivated=='Reactivated'

# keep cell ids for test
test.jeager.idx <- sample(combined.meta$CellID[jeager.react.idx], size = as.integer(sum(jeager.react.idx ) /4) )
test.jeager.idx <- c(test.jeager.idx, sample(combined.meta$CellID[!(jeager.react.idx)], size = as.integer(sum(!(jeager.react.idx) ) /4) ))
test.labels <- combined.meta$Reactivated[combined.meta$CellID %in% test.jeager.idx]

# do the same for chen
chen.glut.idx <- chen2020_meta$BroadCellTypes=='Excitatory'
chen.react.idx <- (chen2020_meta$Reactivated=='Reactivated')

test.chen.idx <- sample(rownames(chen2020_meta)[chen.react.idx&glut.idx], size = as.integer(sum(chen.react.idx&glut.idx) /4) )
test.chen.idx <- c(test.chen.idx, sample(rownames(chen2020_meta)[(!(chen.react.idx))&glut.idx], size = as.integer(sum((!(chen.react.idx))&glut.idx ) /4) ))
test.labels <- c(test.labels, combined.meta$Reactivated[rownames(chen2020_meta) %in% test.chen.idx])

# doing it for amygdala data
amygdalaFC2018_meta$cellID
amy.glut.idx <- amygdalaFC2018_meta$broadcellclass=="glutamatergic"
amy.react.idx <- (amygdalaFC2018_meta$Reactivated=='Reactivated')
test.amy.idx <- sample(amygdalaFC2018_meta$cellID[amy.react.idx&amy.glut.idx], size = as.integer(sum(amy.react.idx&amy.glut.idx) /4) )
test.amy.idx <- c(test.amy.idx, sample(amygdalaFC2018_meta$cellID[(!(amy.react.idx))&amy.glut.idx], size = as.integer(sum( (!(amy.react.idx) )&amy.glut.idx ) /4) ))
test.labels <- c(test.labels, amygdalaFC2018_meta$Reactivated[amygdalaFC2018_meta$cellID %in% test.amy.idx])


# training data
#rememebr to include not in test as a condition
# pick some, down sample the more common class then up sample the less comon class and 
# then train it on that, no cross validation for now
not.in.test <- !(combined.meta$CellID %in% test.jeager.idx)
# targets
train.jeager.idx <- combined.meta$CellID[jeager.react.idx&not.in.test]
# create list so we know what it's label is for classifier
train.labels <- combined.meta$Reactivated[jeager.react.idx&not.in.test] 
# controls
train.jeager.idx <- c(train.jeager.idx, combined.meta$CellID[(!jeager.react.idx)&not.in.test] )
# add to list so we know what it's label is for classifier
train.labels <- c(train.labels, combined.meta$Reactivated[(!jeager.react.idx)&not.in.test])

# chen
not.in.test <- !(rownames(chen2020_meta) %in% test.chen.idx)
# targets
train.chen.idx <- rownames(chen2020_meta)[chen.react.idx&not.in.test&chen.glut.idx]
# add to list so we know what it's label is for classifier
train.labels <- c(train.labels, chen2020_meta$ActivityStatus[chen.react.idx&not.in.test&chen.glut.idx])
# controls
train.chen.idx <- c(train.chen.idx, rownames(chen2020_meta)[(!chen.react.idx)&not.in.test&chen.glut.idx] )
# add to list so we know what it's label is for classifier
train.labels <- c(train.labels, chen2020_meta$ActivityStatus[(!chen.react.idx)&not.in.test&chen.glut.idx])

#amygdala
not.in.test <- !(amygdalaFC2018_meta$cellID %in% test.amy.idx)
#targets
train.amy.idx <- amygdalaFC2018_meta$cellID[amy.react.idx&not.in.test&amy.glut.idx]
# add to list so we know what it's label is for classifier
train.labels <- c(train.labels, amygdalaFC2018_meta$Reactivated[amy.react.idx&not.in.test&amy.glut.idx])
train.amy.idx <- c(train.amy.idx, amygdalaFC2018_meta$cellID[(!amy.react.idx)&not.in.test&amy.glut.idx] )
# add to list so we know what it's label is for classifier
train.labels <- c(train.labels, amygdalaFC2018_meta$Reactivated[(!amy.react.idx)&not.in.test&amy.glut.idx])
