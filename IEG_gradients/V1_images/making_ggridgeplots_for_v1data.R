library(ggplots)
library(ggridges)

tasic2018v1_inhibitory.meta <-read.csv('/home/acampbell/PavLabEngrams/IEG_gradients/V1_images/tasic2018v1_inhibitory_meta.csv')

# fixing the levels of the cell clsuters
tPC1_cluster_ordering <- c('Vip Rspo4 Rxfp1 Chat','Vip Chat Htr1f','Vip Rspo1 Itga4','Vip Crispld2 Htr2c',
'Vip Gpc3 Slc18a3','Vip Col15a1 Pde1a','Vip Arhgap36 Hmcn1','Vip Ptprt Pkp2',
'Vip Igfbp6 Car10','Vip Igfbp4 Mab21l1','Vip Pygm C1ql1','Vip Lmo1 Fam159b',
'Vip Lect1 Oxtr','Vip Igfbp6 Pltp','Vip Lmo1 Myl1','Sncg Vip Nptx2','Vip Crispld2 Kcne4',
'Lamp5 Ntn1 Npy2r','Lamp5 Fam19a1 Tmem182','Lamp5 Fam19a1 Pax6','Sncg Vip Itih5',
'Lamp5 Plch2 Dock5','Sncg Gpr50','Lamp5 Krt73','Lamp5 Lsp1','Lamp5 Lhx6','Sncg Slc17a8',
'Sst Myh8 Fibin','Sst Chodl','Sst Nts','Sst Tac1 Tacr3','Sst Chrna2 Glra3','Sst Calb2 Pdlim5',
'Pvalb Vipr2','Sst Esm1','Pvalb Gabrg1','Sst Chrna2 Ptgdr','Sst Rxfp1 Prdm8','Sst Hpse Sema3c',
'Sst Tac2 Tacstd2','Sst Crh 4930553C11Rik','Sst Myh8 Etv1','Sst Rxfp1 Eya1','Sst Calb2 Necab1',
'Sst Nr2f2 Necab1','Sst Tac2 Myh4','Pvalb Th Sst','Sst Tac1 Htr1d','Sst Crhr2 Efemp1',
'Sst Hpse Cbln4','Sst Mme Fam114a1','Pvalb Akr1c18 Ntf3','Pvalb Reln Tac1','Pvalb Calb1 Sst',
'Pvalb Gpr149 Islr','Pvalb Sema3e Kank4','Pvalb Tpbg','Pvalb Reln Itm2a')

tasic2018v1_inhibitory.meta$cell_cluster <- factor(tasic2018v1_inhibitory.meta$cell_cluster,
                                                      levels= tPC1_cluster_ordering)

# IGf practice
p <-ggplot(tasic2018v1_inhibitory.meta, aes(x=Igf1, y= cell_cluster, fill = stat(x))) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
  scale_fill_viridis_c(name = "UMI LogNorm Scaled", option = "C") +
  ylab('Cell Subtype (Cluster)')+
  xlab('Igf1 Expression along PC1')
p

# Dusp1 Practice
p <-ggplot(tasic2018v1_inhibitory.meta, aes(x=Dusp1, y= cell_cluster, fill = stat(x))) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
  scale_fill_viridis_c(name = "UMI LogNorm Scaled", option = "C") +
  ylab('Cell Subtype (Cluster)')+
  xlab('Dusp1 Expression along PC1')+
  theme_ridges(center_axis_labels = TRUE)
p

tasic2018v1_inhibitory.meta$iegsum_activity

# Summed Normalized IEG expression
p <-ggplot(tasic2018v1_inhibitory.meta, aes(x=iegsum_activity, y= cell_cluster, fill=0.5 - abs(0.5 - stat(ecdf)))) +
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
  scale_fill_viridis_c(name = "Tail probability", direction = -1, option = "C") +
  ylab('Cell Subtype (Cluster)')+
  xlab('Summed Normed IEG expr.')+
  theme_ridges(center_axis_labels = TRUE)
p


# Dusp1
p <-ggplot(tasic2018v1_inhibitory.meta, aes(x=Dusp1, y= cell_cluster, fill=0.5 - abs(0.5 - stat(ecdf)))) +
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
  scale_fill_viridis_c(name = "Tail probability", direction = -1, option = "C") +
  ylab('Cell Subtype (Cluster)')+
  xlab('Dusp1 logp1(median(umi)*(x/umi_x))')+
  theme_ridges(center_axis_labels = TRUE)
p

# Igf1
p <-ggplot(tasic2018v1_inhibitory.meta, aes(x=Igf1, y= cell_cluster, fill=0.5 - abs(0.5 - stat(ecdf)))) +
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
  scale_fill_viridis_c(name = "Tail probability", direction = -1, option = "C") +
  ylab('Cell Subtype (Cluster)')+
  xlab('Igf1 logp1(median(umi)*(x/umi_x))')+
  theme_ridges(center_axis_labels = TRUE)
p

# Arc
p <-ggplot(tasic2018v1_inhibitory.meta, aes(x=Arc, y= cell_cluster, fill=0.5 - abs(0.5 - stat(ecdf)))) +
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
  scale_fill_viridis_c(name = "Tail probability", direction = -1, option = "C") +
  ylab('Cell Subtype (Cluster)')+
  xlab('Arc logp1(median(umi)*(x/umi_x))')+
  theme_ridges(center_axis_labels = TRUE)
p

# Fos
p <-ggplot(tasic2018v1_inhibitory.meta, aes(x=Fos, y= cell_cluster, fill=0.5 - abs(0.5 - stat(ecdf)))) +
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
  scale_fill_viridis_c(name = "Tail probability", direction = -1, option = "C") +
  ylab('Cell Subtype (Cluster)')+
  xlab('Fos logp1(median(umi)*(x/umi_x))')+
  theme_ridges(center_axis_labels = TRUE)
p

# Npas4
p <-ggplot(tasic2018v1_inhibitory.meta, aes(x=Npas4, y= cell_cluster, fill=0.5 - abs(0.5 - stat(ecdf)))) +
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
  scale_fill_viridis_c(name = "Tail probability", direction = -1, option = "C") +
  ylab('Cell Subtype (Cluster)')+
  xlab('Npas4 logp1(median(umi)*(x/umi_x))')+
  theme_ridges(center_axis_labels = TRUE)
p

# Nr4a1
p <-ggplot(tasic2018v1_inhibitory.meta, aes(x=Nr4a1, y= cell_cluster, fill=0.5 - abs(0.5 - stat(ecdf)))) +
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
  scale_fill_viridis_c(name = "Tail probability", direction = -1, option = "C") +
  ylab('Cell Subtype (Cluster)')+
  xlab('Nr4a1 logp1(median(umi)*(x/umi_x))')+
  theme_ridges(center_axis_labels = TRUE)
p







