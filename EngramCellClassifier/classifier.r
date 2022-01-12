#
#
#
#
#
# first attempt at a classifier
setwd("C:/Users/angus/Desktop/PavLabEngrams/EngramCellClassifier")

library(tidyverse)
library(GEOquery)
library(AnnotationDbi)
library(randomForest)

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

#splitting the data, NOT DONE YET

# Set random seed to make results reproducible:
set.seed(23)
# Calculate the size of each of the data sets:
data_set_size <- floor(nrow(iris)/2)
# Generate a random sample of "data_set_size" indexes
indexes <- sample(1:nrow(iris), size = data_set_size)
# Assign the data to the correct sets
training <- iris[indexes,]
validation1 <- iris[-indexes,]


# Perform training:
rf_classifier = randomForest(Species ~ ., data=training, ntree=100, mtry=2, importance=TRUE)











