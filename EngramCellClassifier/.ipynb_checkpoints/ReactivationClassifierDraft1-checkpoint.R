# training looking for reactivation

# caret pdf
# https://cran.r-project.org/web/packages/caret/caret.pdf

# https://machinelearningmastery.com/multi-class-imbalanced-classification/
# suggests using smote

# how caret manages class imbalances
# this si chosen through sampling option described on pg 171 of their pdf
# https://topepo.github.io/caret/subsampling-for-class-imbalances.html

#sorting labels
arclabels <- read.csv("Jeager2018_GSE98679/Jaeger2018_meta_arclabels.csv", header = TRUE)

combined.meta <- left_join(combined.meta,
                           arclabels %>% dplyr::select(Title, ArcStatus),
                           by = c("CellID" = "Title"))


jeagerdata.multilabels.hg <- labeled.data.hg[,c(1:dim(labeled.data.hg)[2])]




activation <- 

jeagerdata.multilabels.hg$Activation


ctrl <- trainControl(method = "repeatedcv", repeats = 5,
                     classProbs = TRUE,
                     summaryFunction = multiClassSummary,
                     classProbs = TRUE,
                     ## new option here:
                     sampling = "up")

model_with_down_sample <- train(Engramcell ~ ., data = imbal_train,
                                method = "RRF",
                                preProcess = c("range"),
                                verbose = FALSE,
                                trControl = ctrl)

# caret tutorial using random forest:
# https://www.machinelearningplus.com/machine-learning/caret-package/

# multi class classfication in caret
# https://stats.stackexchange.com/questions/453259/multi-class-probabilities-of-random-forest-inside-caret-model

# cv and resampling
# https://stackoverflow.com/questions/45250252/how-to-downsample-using-r-caret

# cv stakc question in caret just for good measure
# https://stackoverflow.com/questions/22909197/creating-folds-for-k-fold-cv-in-r-using-caret

# how to do it in random forest
# https://rpubs.com/jkylearmstrong/RF_Imputation_Multi_class

# models available in caret include the regularized random forest
# https://rdrr.io/cran/caret/man/models.html
# Regularized Random Forest (method = 'RRF')
# Regularized Random Forest (method = 'RRFglobal')
# they are the same and based on RRF model but they do it 
# with 