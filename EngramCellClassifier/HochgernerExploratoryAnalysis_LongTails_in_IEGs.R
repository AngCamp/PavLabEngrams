#  Hochgerner Exploratory Analysis
# Show long tails of Hochgerner IEGS

# Exploratory Analysis of Hochgerner
# Making some plots showing long tails of cannonical IEGs
# http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/81-ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page/
library(ggpubr)
library(ggplot2)
library(dplyr)


# > which(rownames(hochgerner5k_2018_counts)=="Arc")
# [1] 1213
# > which(rownames(hochgerner5k_2018_counts)=="Fos")
# [1] 4716
# > which(rownames(hochgerner5k_2018_counts)=="Inhba")
# [1] 6470
# > which(rownames(hochgerner5k_2018_counts)=="Nptx2")
# [1] 8583

df <- t(hochgerner5k_2018_counts[c(1213, 4716, 6470, 8583), hoch5k.GC_Adult.p35.idx])
df <- apply(df, 2,as.numeric)
df <- data.frame(df)

p.fos <- ggplot(data = df, aes(x=Fos) )
p.arc <- ggplot(data = df, aes(x=Arc) )
p.inhba <- ggplot(data = df, aes(x=Inhba) )
p.nptx2 <- ggplot(data = df, aes(x=Nptx2) )

dev.off()
jpeg("hochgernerDGCs_foscounts.jpg", width = 700, height = 700)
p.fos + geom_histogram(color = "darkgreen", fill = "lightgreen") + theme_classic() +
  xlab("Reads") + ylab ("Number of Cells") + ggtitle("Fos") +
  theme(axis.text.x=element_text(size=15, face = "bold")) + 
  theme(axis.text.y=element_text(size=15, face = "bold")) +
  theme(axis.title.x=element_text(size=20, face = "bold")) + 
  theme(axis.title.y=element_text(size=20, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.5, size=22, face = "bold"))
dev.off()


dev.off()
jpeg("hochgernerDGCs_ARCcounts.jpeg", width = 700, height = 700)
p.arc + geom_histogram(color = "darkgreen", fill = "lightgreen") + theme_classic() +
  xlab("Reads") + ylab ("Number of Cells") + ggtitle("Arc") +
  theme(axis.text.x=element_text(size=15, face = "bold")) + 
  theme(axis.text.y=element_text(size=15, face = "bold")) +
  theme(axis.title.x=element_text(size=20, face = "bold")) + 
  theme(axis.title.y=element_text(size=20, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.5, size=22, face = "bold"))
dev.off()


dev.off()
jpeg("hochgernerDGCs_Inhbacounts.jpeg", width = 700, height = 700)
p.inhba + geom_histogram(color = "darkgreen", fill = "lightgreen") + theme_classic() +
  xlab("Reads") + ylab ("Number of Cells") + ggtitle("Inhba") +
  theme(axis.text.x=element_text(size=15, face = "bold")) + 
  theme(axis.text.y=element_text(size=15, face = "bold")) +
  theme(axis.title.x=element_text(size=20, face = "bold")) + 
  theme(axis.title.y=element_text(size=20, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.5, size=22, face = "bold"))
dev.off()


dev.off()
jpeg("hochgernerDGCs_Nptx2counts.jpeg", width = 700, height = 700)
p.nptx2 + geom_histogram(color = "darkgreen", fill = "lightgreen") + theme_classic() +
  xlab("Reads") + ylab ("Number of Cells") + ggtitle("Nptx2") +
  theme(axis.text.x=element_text(size=15, face = "bold")) + 
  theme(axis.text.y=element_text(size=15, face = "bold")) +
  theme(axis.title.x=element_text(size=20, face = "bold")) + 
  theme(axis.title.y=element_text(size=20, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.5, size=22, face = "bold"))
dev.off()



