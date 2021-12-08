#draft ofr the classifier
setwd("C:/Users/angus/Desktop/PavLabEngrams/EngramCellClassifier")

library(tidyverse)
library(GEOquery)
library(edgeR)
#library(umap)
library(limma)
library(knitr)
#library(readxl)
library(eulerr)
library(ermineR)
library(sgof)
library(AnnotationDbi)
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

# Jaeger, B. N., Linker, S. B., Parylak, S. L., Barron, J. J., Gallina, I. S., Saavedra, C. D., ... & Gage, F. H. (2018). 
# A novel environment-evoked transcriptional signature predicts reactivity in single dentate granule neurons.
# Nature communications, 9(1), 1-15.
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6079101/
gds_jeager2018 <- getGEO("GSE98679")
jeager2018_counts <- read.table('Jeager2018_GSE98679/GSE98679_count.txt.gz')

test <- read.table('Jeager2018_GSE98679/GSE98679_v2_GSM3308862-GSM3309413_tpm.txt.gz')


GSE98679

C:/Users/angus/Desktop/PavLabEngrams/EngramCellClassifier/GSE152632

#gives error 
Found 1 file(s)
GSE152632_series_matrix.txt.gz
trying URL 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE152nnn/GSE152632/matrix/GSE152632_series_matrix.txt.gz'
Content type 'application/x-gzip' length 91125 bytes (88 KB)
downloaded 88 KB

Error: The size of the connection buffer (131072) was not large enough
to fit a complete line:
  * Increase it by setting `Sys.setenv("VROOM_CONNECTION_SIZE")`
