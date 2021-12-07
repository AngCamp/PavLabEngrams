#draft ofr the classifier

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


#get Chen et al., (2020)
gds <- getGEO("GSE152632")

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
