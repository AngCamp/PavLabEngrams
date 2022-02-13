#
# first attempt at a classifier
setwd("C:/Users/angus/Desktop/PavLabEngrams/EngramCellClassifier")

library(tidyverse)
library(GEOquery)
library(AnnotationDbi)
library(randomForest)
library(data.table)
library(reshape2)
library(FactoMineR)
library(factoextra)
library(Rtsne)
library(dplyr)
library(Seurat)
library(stringr)
library(patchwork)
library(metap)
library(cqn)
library(GeoTcgaData)
BiocManager::install("cqn")

# Get a profile, of the engram cells just look at recent 
# 

# Lacar and Jeager --------------------------------------------------------


# Lacar et al., (2016)
lacar2016_meta <- read.csv('Lacar2016_GSE77067/SraRunTable.txt', header = TRUE)
lacar2016_snHC_counts <- read.table('Lacar2016_GSE77067/GSE77067_sn_hc_counts.txt.gz')
lacar2016_snNE_counts <- read.table('Lacar2016_GSE77067/GSE77067_sn_ne_counts.txt.gz')



# Jeager et al., (2018) meta data and counts
jeager2018_counts <- bind_cols(read.table('Jeager2018_GSE98679/GSE98679_count.txt.gz', header = TRUE, check.names = FALSE),
                               read.table('Jeager2018_GSE98679/GSE98679_v2_GSM3308862-GSM3309413_count.txt.gz', header = TRUE, check.names = FALSE))

jeager2018_meta <- read.csv('Jeager2018_GSE98679/SraRunTable.txt', header = TRUE)
jeager2018_meta = jeager2018_meta[c(1:46,599:912,47:598),] #we need to fix Jeager's files up a bit
rownames(jeager2018_meta) <- c(1:912)


# blank.to.DG <-function(x){
#   if(x==""){
#     x="DG"
#   }
# }

jeager2018_meta$predicted_cell_type <- as.character(lapply(jeager2018_meta$predicted_cell_type, function(x) if (x=="") {"DG"} else {x}))
jeager2018_meta$predicted_cell_type <- lapply(jeager2018_meta$predicted_cell_type, function(x) if (x=="") {"DG"} else {x})

#Finding engram cells
fospos <- which(grepl("_F_DG",  jeager2018_meta$source_name))
fospos <- c(fospos,361:912) # since we know all the v2 cells from time points and recall testing

neg <- which(grepl("_N_DG",  jeager2018_meta$source_name))

jeager2018_meta$fos_status <- as.factor(sapply(as.character(jeager2018_meta$source_name), function(y) if (grepl("_F_", y, fixed=TRUE)) "Fos+" else "Fos-"  ))
jeager2018_meta$fos_status[361:dim(jeager2018_meta)[1]] <- "Fos+"

which((jeager2018_meta$exposure=="HC")&(jeager2018_meta$fos_status=="Fos+"))


#filtering out cells as per instructions in Jeager et al., 2018)
under4k <- sapply(jeager2018_counts, function(y) sum(length(which(y>0))))

filtered.idx <- as.numeric(which(under4k>4000))
filtered.idx <- order(c(filtered.idx,194))

jeager2018_counts <-jeager2018_counts[,filtered.idx]
jeager2018_meta <- jeager2018_meta[filtered.idx,]


#Make seurate object
jeager2018_counts[is.na(jeager2018_counts)] = 0
jeager2018 <- CreateSeuratObject(counts = log(jeager2018_counts+1))

jeager2018 <- FindVariableFeatures(jeager2018, selection.method = "vst", nfeatures = 3530) #no idea how I chose 3530, 
# in the tutorial they chose 2000, I just chose to match the number of cells we had
#jeager2018nobatch <- ScaleData(jeager2018)

rownames(jeager2018_meta) <- colnames(jeager2018_counts)
jeager2018 <- AddMetaData(jeager2018, jeager2018_meta)




#Writting the Classifier: https://www.rdocumentation.org/packages/Seurat/versions/2.3.4/topics/ClassifyCells



#First we do integration, find anchors etc

#split the data by experiment 

# D0 not do now but I should be trying to recreate the resutls of the 
# DEG lists from the original authors and trygint to restrict the classifier to those
# genes  Suerat has a tutorial on getting the DEGs out of it. https://satijalab.org/seurat/articles/de_vignette.html



# Following this we can run the data through this function
ClassifyCells(jeagerdata, classifier, training.genes = NULL,
              training.classes = NULL, new.data = test.set, ...)


# Next steps?  What would be a good positive control?






### DG dataset to check:
# Hochgerner, H., Zeisel, A., Lönnerberg, P., & Linnarsson, S. (2018). 
# Conserved properties of dentate gyrus neurogenesis across postnatal 
# development revealed by single-cell RNA sequencing. Nature neuroscience, 
# 21(2), 290-299.
#GEO Accession
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE95753
#meta data on SRA
# https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA378266&o=acc_s%3Aa
~32000 cells, NOTE THE DATASET IS TOO LARGE DOWNLOAD IT LATER





Zhang, T. Y., Keown, C. L., Wen, X., Li, J., Vousden, D. A., Anacker, C., ... & Meaney, M. J. (2018).
Environmental enrichment increases transcriptional and epigenetic differentiation
between mouse dorsal and ventral dentate gyrus. Nature communications, 9(1), 1-11.

-methylation data but may be interestingto provide some context to our results,
worth reading



# Habib, N., Li, Y., Heidenreich, M., Swiech, L., Avraham-Davidi, I., 
# Trombetta, J. J., ... & Regev, A. (2016). Div-Seq: Single-nucleus RNA-Seq reveals
# dynamics of rare adult newborn neurons. Science, 353(6302), 925-928.
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84371

# data here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE85721

testsetpath <- "C:/Users/angus/Desktop/test_datasets"

test_habib2016_logTPM <- read.table(paste(testsetpath,"/Habib2016/GSE85721_PROCESSED_data.sNuc-Seq_Data_RSEM_log_TPM.for_GEO.txt.gz",sep=""), header = TRUE)
rownames(test_habib2016_logTPM) <- test_habib2016_logTPM[, 1]  ## set rownames
test_habib2016_logTPM <- test_habib2016_logTPM[, 2:925]   


# Hochgerner 2018 Development of DG ---------------------------------------
#calling these in the rstudieo browser: wget url 
# URLs for Hochgerner
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE104323 #
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE95752 # 
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE95315 #downloaded 5000 cells

hochgerner5k_2018_counts <- read.table(paste(testsetpath,"/Hochgerner2018/GSE95315_10X_expression_data_v2.tab.gz", sep=""))
hochgerner24k_2018_counts <- read.table(paste(testsetpath,"/Hochgerner2018/GSE104323_10X_expression_data_V2.tab.gz", sep=""))
colnames(hochgerner24k_2018_counts) <- hochgerner24k_2018_counts[1,]
rownames(hochgerner24k_2018_counts) <- hochgerner24k_2018_counts[,1]
hochgerner24k_2018_counts <- hochgerner24k_2018_counts[2:27934,2:24186]

hochgerner24k_2018_meta <- read.table(paste(testsetpath,"/Hochgerner2018/GSE104323_metadata_barcodes_24185cells.txt.gz", sep=""), 
                                      sep ="\t", header = T)

#check the publication they pooled cells and then tried to determine sex based on 
# expression see this methods section Dataset A: allocation of sampling age to cell

#this lets you check number of each cell type
table(as.character(hochgerner5k_2018_counts[3,as.character(hochgerner5k_2018_counts[2,1:5])=="35*"]))


# Testing for markers from Jeager 2018 ------------------------------------

# Summary of findings from breifly searching for the top gene's lsited as 
# persistant markers of early activation in Jeager et al., (2018) in the final results 
# section and discussion.  Jund, Tet3 Atf3 and Gadd45b are listed as the top most
# markers for this.  I also looked for arc and fos.  Jund is particularily interesting because it
# fors part of the AP1 T.F. complex which arc also forms part of and therefor upregulates.  genes invovled in long term persisten activation 
# are Blnk, Tnik, Acam Entpd4, Hes7.  
# 
# In Habib et al., (2016) is difficult to compare as it is listed in logTPM.
# few cells marked by arc. I need to check this studies average reads per cell.  600 granule cells
# were in the study and few had arc and fox co expressed, distributions of arc were not promissing or comperable
# to Jeager.  Howeve it should be noted that Jeager had a average read per cell of 1.5 million this amy be the problem for detecting Arc.
# When i looked for presence of all persitant markerss simultanouesly not cells passed threshold  but jsut looking at
# arc + one gene I found arc + Blnk 10, arc+ tnik 21,  
# 
# In Hochenberger et al., (2018) I looked at the second dataset but not by age and found a nice distribution fo fos but not for
# arc again however the average reads per cell of this study was 43000, 3x10^2 times less power.  I did find examples of
# cells which expressed arc and one of the persistant markers but increasing the number of persistant markers caused the number
# of cells passing threshold to drop dramatically.  Arc was again distributed fairly poorly, though fos was prestn at levels that 
# made sense give the reads per cell and was distrbuted the same as in the Jeager study.  One interesting thing to note was that 
# in immature granule cells both fos and arc were expressed in much greater abundance.

#Checking co-presence of the persistant markers from Jeager
engram_markers <- c("Arc","Fos","Jund", "Tet3", "Atf3", "Gadd45b")

marker.idx <- c()
for(i in engram_markers){
  marker.idx <-c(marker.idx,(which(rownames(test_habib2016_logTPM)==i)))
}

persistant_gene <- test_habib2016_logTPM[ engram_markers, habib2016.dg.idx]
persistant_gene <- data.frame(t(persistant_gene))

sum((persistant_gene$Hes7>0)&(persistant_gene$Arc>0))

#hochgerner first data set
marker.idx <- c()
for(i in engram_markers){
  marker.idx <-c(marker.idx,(which(rownames(hochgerner2018_counts)==i)))
}

hochgerner2018.dg.idx <-which(as.character(hochgerner2018_counts[3,])=="DG")

persistant_gene.hoch <- hochgerner2018_counts[as.numeric(marker.idx),as.numeric(hochgerner2018.dg.idx)]
persistant_gene.hoch <- data.frame(t(persistant_gene.hoch))


sum((persistant_gene.hoch$Fos>0)&(persistant_gene.hoch$Jund>0))


#hochgerner second larger dataset try with Immature-GC or GC-adult
hochgerner24k.dg.idx <-which(as.character(hochgerner24k_2018_meta$characteristics..cell.cluster)=="GC-adult")

marker.idx <- c()
for(i in engram_markers){
  marker.idx <-c(marker.idx,(which(rownames(hochgerner24k_2018_counts)==i)))
}

persistant_gene.hoch24k <- hochgerner24k_2018_counts[marker.idx, hochgerner24k.dg.idx]
persistant_gene.hoch24k <- data.frame(t(persistant_gene.hoch24k))


sum((persistant_gene.hoch24k$Fos>0)&(persistant_gene.hoch24k$Arc>0))

rownames(hochgerner24k_2018_counts)[str_detect(rownames(hochgerner24k_2018_counts), "Arc")]


###############reactivation

#habib in tpm remember
reactivation.markers <- c("Arc", "Fos", "Blnk", "Tnik", "Acan", "Entpd4", "Hes7")

marker.idx <- c()
for(i in reactivation.markers){
  marker.idx <-c(marker.idx,(which(rownames(hochgerner2018_counts)==i)))
}

persistant_gene <- test_habib2016_logTPM[ reactivation.markers, habib2016.dg.idx]
persistant_gene <- data.frame(t(persistant_gene))

sum((persistant_gene$Arc>0)&(persistant_gene$Hes7>0))



#Hochgerner2018 first dataset also try with "Immature-GC"
hochgerner2018.dg.idx <-which(as.character(hochgerner2018_counts[3,])=="GC-adult")

persistant_gene.hoch <- hochgerner2018_counts[marker.idx,hochgerner2018.dg.idx]
persistant_gene.hoch <- data.frame(t(persistant_gene.hoch))


sum((persistant_gene.hoch$Arc>0)&(persistant_gene.hoch$Acan>0))
sum(persistant_gene.hoch[1216,]>0)

#hochgerner second larger dataset
hochgerne24k.dg.idx <-which(as.character(hochgerner24k_2018_meta$characteristics..cell.cluster)=="GC-adult")

marker.idx <- c()
for(i in reactivation.markers){
  marker.idx <-c(marker.idx,(which(rownames(hochgerner2018_counts)==i)))
}

persistant_gene.hoch24k <- hochgerner2018_counts[marker.idx,hochgerne24k.dg.idx]
persistant_gene.hoch24k <- data.frame(t(persistant_gene.hoch24k))


sum((persistant_gene.hoch24k$Gadd45b>0)&(persistant_gene.hoch$Tet3>0))



#note that based on the colnames there appear to be some CA1 and other cells in here,
#as well as some conditions, old young etc that will need parsing, read the pub
# then get on it
condition_label <- colnames(test_habib2016_logTPM)
habib2016.dg.idx <- which(str_detect(colnames(test_habib2016_logTPM), ".DG"))
# [1] 545 pf DG cells not a bad number  of cells to test on


test_habib2016_meta <- read.csv(paste(testsetpath,"/Habib2016/SraRunTable.txt",sep=""), header =TRUE, fill =TRUE)
rownames(test_habib2016_meta) <- colnames(test_habib2016_logTPM)





# Integrating the data ----------------------------------------------------

#we must convert the jeager data to log tpm
# note that tpm may be  hihgly problematic for drawing conclusions from,
# we need to find the count data for this habib dataset

test <- countToTpm_matrix(jeager2018_counts)
test <- 
test <- merge(log(countToTpm_matrix(jeager2018_counts)), test_habib2016_logTPM, by = "row.names")



#Don't be discouraged if this does not work
data.list <- c(CreateSeuratObject(counts = log(countToTpm_matrix(jeager2018_counts[shared.genes[1]])),
                                  meta.data = jeager2018_meta),
               CreateSeuratObject(counts = test_habib2016_logTPM[shared.genes[2],habib2016.dg.idx],
                                  meta.data = test_habib2016_meta[habib2016.dg.idx,])
              )#end of data.list


jeager2018_meta$

               
data.list <- CreateSeuratObject(data.list)
data.list <- FindVariableFeatures(data.list)


