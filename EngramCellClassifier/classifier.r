#draft ofr the classifier
setwd("C:/Users/angus/Desktop/PavLabEngrams/EngramCellClassifier")

#https://www.youtube.com/watch?v=AHW0FfHDulI&ab_channel=Medeos

#Lets try to train a support vector machine classifier
#https://www.datacamp.com/community/tutorials/support-vector-machines-r

library(tidyverse)
library(GEOquery)
library(AnnotationDbi)

#if Data sets are needed pull them from the following excession numbers



library(edgeR)
#library(umap)
library(limma)
library(knitr)
#library(readxl)
library(eulerr)
library(ermineR)
library(sgof)

#library(org.Mm.eg.db)

if (!require('devtools')) install.packages('devtools'); require('devtools')
# make sure you have Rtools installed first! if not, then run:
#install.packages('installr')
#install_Rtools()
devtools::install_github('talgalili/installr')


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GEOquery")


#get Chen et al., (2020)
#important before getting the geo use this...
Sys.setenv(VROOM_CONNECTION_SIZE = 1e6)

gds_chen2020 <- getGEO("GSE152632")

chen2020_counts <- read.table('Chen2020_GSE152632/GSE152632_GEO_mberkchen_TRAP2_counts.csv.gz')
chen2020_counts <- read.csv('Chen2020_GSE152632/GSE152632_GEO_mberkchen_TRAP2_counts.csv.gz', header = TRUE)
chen2020_meta <- read.csv( 'Chen2020_GSE152632/SraRunTable.txt', header = TRUE)
#I could just look at the data using known cell markers myself 


# Jaeger, B. N., Linker, S. B., Parylak, S. L., Barron, J. J., Gallina, I. S., Saavedra, C. D., ... & Gage, F. H. (2018). 
# A novel environment-evoked transcriptional signature predicts reactivity in single dentate granule neurons.
# Nature communications, 9(1), 1-15.
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6079101/
#gds_jeager2018 <- getGEO("GSE98679")
jeager2018_counts <- read.table('Jeager2018_GSE98679/GSE98679_count.txt.gz', header = TRUE, check.names = FALSE)
jeager2018_v2_counts <- read.table('Jeager2018_GSE98679/GSE98679_v2_GSM3308862-GSM3309413_count.txt.gz', header = TRUE, check.names = FALSE)


#> which(rownames(jeager2018_counts)=="Arc")
#[1] 3021

#> which(rownames(jeager2018_counts)=="Fos")
#[1] 8233

#> which(rownames(jeager2018_counts)=="Prox1")
#[1] 35227

jeager2018_v2_tpm <- read.table('Jeager2018_GSE98679/GSE98679_v2_GSM3308862-GSM3309413_tpm.txt.gz')
jeager2018_meta <- read.csv('Jeager2018_GSE98679/SraRunTable.txt', header = TRUE)

#meta data is poorly detailed, but essentially in exposure the v2 counts are described

#note that some of the meta data, specifically the stain is in the cell barcodes in the count data


#lacar2016
lacar2016_meta <- read.csv('Lacar2016_GSE77067/SraRunTable.txt', header = TRUE)
lacar2016_wc_counts <- read.table('Lacar2016_GSE77067/GSE77067_wc_counts.txt.gz')
lacar2016_snHC_counts <- read.table('Lacar2016_GSE77067/GSE77067_sn_hc_counts.txt.gz')
lacar2016_snNE_counts <- read.table('Lacar2016_GSE77067/GSE77067_sn_ne_counts.txt.gz')


#Cho et al., 2016
cho2016_data <- read.table('Cho2016_GSE85128/GSE85128_SuppTable1_normalized_read_counts_resubmission.xls.gz')

#Rao-Ruiz
rao2019_meta=read.csv('Rao-Ruiz2019_GSE129024/SraRunTable.txt', header=TRUE)

#in the meta data the source_name columns should allow for merging, though
# the colnames in jeager2018_counts may need to be parsed

#Jeagar used Monocle for their analyses: http://cole-trapnell-lab.github.io/monocle-release/

BiocManager::install("monocle")
# we can also use seurat
# https://satijalab.org/seurat/articles/essential_commands.html

#clustering with Seurat, how to find cell types
# https://www.youtube.com/watch?v=5eYcLfs2L0g&ab_channel=LiquidBrain

GSE98679


#gives error 
Found 1 file(s)
GSE152632_series_matrix.txt.gz
trying URL 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE152nnn/GSE152632/matrix/GSE152632_series_matrix.txt.gz'
Content type 'application/x-gzip' length 91125 bytes (88 KB)
downloaded 88 KB

Error: The size of the connection buffer (131072) was not large enough
to fit a complete line:
  * Increase it by setting `Sys.setenv("VROOM_CONNECTION_SIZE")`
