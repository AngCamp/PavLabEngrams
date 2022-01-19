library(FactoMineR)
library(factoextra)
library(ggplot2)
library(tidyr)
library(dplyr)
library(MASS)
library(reshape2)
library(cowplot)

data("diamonds")
dat <- diamonds %>% data.frame
head(dat)
dat <- dat[sample(rownames(dat), 2000),]
dat <-  dat %>% filter (x > 0, y > 0, z > 0)
for (i in 1:length(colnames(dat))){
  
  if (is.numeric(dat[, i])==TRUE)
    
    dat[, i] <- as.numeric(scale(dat[, i]))
  
  else
    
    dat[, i] <- dat[, i]
  
}


# Run the PCA
pca1 <- PCA(dat[ ,c("carat", "depth", "table", "price", "x", "y", "z", "clarity", "cut", "color")], 
            quali.sup = c(8:10), graph = FALSE)

#Visualizing the PCA using ggplot2
# extract pc scores for first two component and add to dat dataframe
dat$pc1 <- pca1$ind$coord[, 1] # indexing the first column

dat$pc2 <- pca1$ind$coord[, 2]  # indexing the second column

#We also need to extract the data for the variable 
#contributions to each of the pc axes.
pca.vars <- pca1$var$coord %>% data.frame

pca.vars$vars <- rownames(pca.vars)

pca.vars.m <- melt(pca.vars, id.vars = "vars")


# By convention, the variable contribution plot has a 
# circle around the variables that has a radius of 1. 
# Here's some code to make one.

#we dont care about this

circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

circ <- circleFun(c(0,0),2,npoints = 500)

#initial pca plot

p <- ggplot(data = dat, aes(x = pc1, y = pc2)) +
  
  geom_point()
p

# And we can customize it a bit.
# For example, coloring and shaping the points by cut.

p <- ggplot(data = dat, aes(x = pc1, y = pc2, color = cut, shape = cut)) +
  
  geom_hline(yintercept = 0, lty = 2) +
  
  geom_vline(xintercept = 0, lty = 2) +
  
  geom_point(alpha = 0.8) 


p