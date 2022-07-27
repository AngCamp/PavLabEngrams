# this is for loading the immature granule data set from Zhou et al., (2022)
# a snRNA-seq Dentate gyrus dataset collected by SPLiT-Seq from human epileptics

# GSE185553 # URL:https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE185553
# Run,Age,Assay Type,AvgSpotLen,Bases,BioProject,BioSample,Bytes,Center Name,Consent,DATASTORE filetype,DATASTORE provider,DATASTORE region,Development_stage,Experiment,GEO_Accession (exp),Instrument,LibraryLayout,LibrarySelection,LibrarySource,Organism,Platform,ReleaseDate,Sample Name,source_name,SRA Study,tissue
# SRR16248559,4.3 Yrs,RNA-Seq,160,16194169120,PRJNA769637,SAMN22155633,7519639884,GEO,public,"fastq,sra","gs,ncbi,s3","gs.US,ncbi.public,s3.us-east-1",Child,SRX12528354,GSM5618237,NextSeq 550,PAIRED,cDNA,TRANSCRIPTOMIC,Homo sapiens,ILLUMINA,2022-04-20T00:00:00Z,GSM5618237,human hippocampus,SRP340513,human hippocampus
# SRR16248560,4.3 Yrs,RNA-Seq,160,19381193600,PRJNA769637,SAMN22155633,8915830300,GEO,public,"fastq,sra","gs,ncbi,s3","gs.US,ncbi.public,s3.us-east-1",Child,SRX12528354,GSM5618237,NextSeq 550,PAIRED,cDNA,TRANSCRIPTOMIC,Homo sapiens,ILLUMINA,2022-04-20T00:00:00Z,GSM5618237,human hippocampus,SRP340513,human hippocampus
# SRR16248561,4 Yrs/15.3 Yrs/59 Yrs/92 Yrs/0.2 Yrs,RNA-Seq,160,94917540320,PRJNA769637,SAMN22155634,44926957693,GEO,public,"fastq,sra","gs,ncbi,s3","gs.US,ncbi.public,s3.us-east-1",Child/Adolescent/Adult/Aging/Infant,SRX12528350,GSM5618238,NextSeq 550,PAIRED,cDNA,TRANSCRIPTOMIC,Homo sapiens,ILLUMINA,2022-04-20T00:00:00Z,GSM5618238,human hippocampus,SRP340513,human hippocampus
# SRR16248562,4 Yrs/15.3 Yrs/59 Yrs/92 Yrs/0.2 Yrs,RNA-Seq,160,96969022560,PRJNA769637,SAMN22155634,46347571775,GEO,public,"fastq,sra","gs,ncbi,s3","gs.US,ncbi.public,s3.us-east-1",Child/Adolescent/Adult/Aging/Infant,SRX12528350,GSM5618238,NextSeq 550,PAIRED,cDNA,TRANSCRIPTOMIC,Homo sapiens,ILLUMINA,2022-04-20T00:00:00Z,GSM5618238,human hippocampus,SRP340513,human hippocampus
# SRR16248563,0.6 Yrs/3.2 Yrs/14.7 Yrs/50 Yrs/90 Yrs,RNA-Seq,160,94180920800,PRJNA769637,SAMN22155635,44979185765,GEO,public,"fastq,sra","gs,ncbi,s3","gs.US,ncbi.public,s3.us-east-1",Infant/Child/Adolescent/Adult/Aging,SRX12528351,GSM5618239,NextSeq 550,PAIRED,cDNA,TRANSCRIPTOMIC,Homo sapiens,ILLUMINA,2022-04-20T00:00:00Z,GSM5618239,human hippocampus,SRP340513,human hippocampus
# SRR16248564,18.5 Yrs/60 Yrs/86 Yrs/16 Yrs/95 Yrs,RNA-Seq,160,93196292160,PRJNA769637,SAMN22155636,43769783401,GEO,public,"fastq,sra","gs,ncbi,s3","gs.US,ncbi.public,s3.us-east-1",Adolescent/Adult/Aging/Adolescent/Aging,SRX12528352,GSM5618240,NextSeq 550,PAIRED,cDNA,TRANSCRIPTOMIC,Homo sapiens,ILLUMINA,2022-04-20T00:00:00Z,GSM5618240,human hippocampus,SRP340513,human hippocampus
# SRR16248565,6.7 Yrs/1.2 Yrs/2.1 Yrs/88 Yrs/4 Yrs,RNA-Seq,160,85766701280,PRJNA769637,SAMN22155637,40089859390,GEO,public,"fastq,sra","gs,ncbi,s3","gs.US,ncbi.public,s3.us-east-1",Child/Infant/Infant/Aging/Child,SRX12528353,GSM5618241,NextSeq 550,PAIRED,cDNA,TRANSCRIPTOMIC,Homo sapiens,ILLUMINA,2022-04-20T00:00:00Z,GSM5618241,human hippocampus,SRP340513,human hippocampus
# 
# GSM5618237	Specimen_7 -both contains the 4.3 jsut on different plates
# GSM5618238	Specimen_6/11/16/22/39
# GSM5618239	Specimen_2/5/10/14/21
# GSM5618240	Specimen_12/17/18/40/41
# GSM5618241	Specimen_8/3/4/19/42

# opening this data: 
# for linux tar -C /GSE185553 -xvf GSE185553_RAW.tar.gz
# in R untar(tarfile = "~/test_datasets/Zhou2022/GSE185277_RAW.tar", exdir = "./GSE185277")
# count data
# ~/test_datasets/Zhou2022
# 
# GSE185277 # URL:https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE185277
# Run,Age,Assay Type,AvgSpotLen,Bases,BioProject,BioSample,Bytes,Center Name,Consent,DATASTORE filetype,DATASTORE provider,DATASTORE region,Development_stage,Experiment,GEO_Accession (exp),Instrument,LibraryLayout,LibrarySelection,LibrarySource,Organism,Platform,ReleaseDate,Sample Name,source_name,SRA Study,tissue
# SRR16206458,13 Yrs,RNA-Seq,160,6766993120,PRJNA768440,SAMN22063841,3069717289,GEO,public,"fastq,sra","gs,ncbi,s3","gs.US,ncbi.public,s3.us-east-1",Adolescent,SRX12490812,GSM5610939,NextSeq 550,PAIRED,cDNA,TRANSCRIPTOMIC,Homo sapiens,ILLUMINA,2022-04-20T00:00:00Z,GSM5610939,human hippocampus,SRP339954,human hippocampus
# SRR16206459,40 Yrs,RNA-Seq,160,81448996320,PRJNA768440,SAMN22063842,36140931633,GEO,public,"fastq,sra","gs,ncbi,s3","gs.US,ncbi.public,s3.us-east-1",Adult,SRX12490813,GSM5610940,NextSeq 550,PAIRED,cDNA,TRANSCRIPTOMIC,Homo sapiens,ILLUMINA,2022-04-20T00:00:00Z,GSM5610940,human hippocampus,SRP339954,human hippocampus
# SRR16206460,50.2 Yrs,RNA-Seq,160,27192161600,PRJNA768440,SAMN22063843,12520429981,GEO,public,"fastq,sra","gs,ncbi,s3","gs.US,ncbi.public,s3.us-east-1",Adult,SRX12490814,GSM5610941,NextSeq 550,PAIRED,cDNA,TRANSCRIPTOMIC,Homo sapiens,ILLUMINA,2022-04-20T00:00:00Z,GSM5610941,human hippocampus,SRP339954,human hippocampus
# SRR16206461,89 Yrs,RNA-Seq,160,41566176160,PRJNA768440,SAMN22063844,19457906173,GEO,public,"fastq,sra","gs,ncbi,s3","gs.US,ncbi.public,s3.us-east-1",Aging,SRX12490815,GSM5610942,NextSeq 550,PAIRED,cDNA,TRANSCRIPTOMIC,Homo sapiens,ILLUMINA,2022-04-20T00:00:00Z,GSM5610942,human hippocampus,SRP339954,human hippocampus
# SRR16194336,0.1 Yr,RNA-Seq,80,12965417920,PRJNA768440,SAMN22043370,5602673886,GEO,public,"fastq,sra","gs,ncbi,s3","gs.US,ncbi.public,s3.us-east-1",Infant,SRX12478833,GSM5609934,NextSeq 550,PAIRED,cDNA,TRANSCRIPTOMIC,Homo sapiens,ILLUMINA,2022-04-20T00:00:00Z,GSM5609934,human hippocampus,SRP339954,hippocampus
# SRR16194335,0.1 Yr,RNA-Seq,80,16272292320,PRJNA768440,SAMN22043370,6971740436,GEO,public,"fastq,sra","gs,ncbi,s3","gs.US,ncbi.public,s3.us-east-1",Infant,SRX12478833,GSM5609934,NextSeq 550,PAIRED,cDNA,TRANSCRIPTOMIC,Homo sapiens,ILLUMINA,2022-04-20T00:00:00Z,GSM5609934,human hippocampus,SRP339954,hippocampus
# 
# GSM5609934	Specimen_1
# GSM5610939	Specimen_9
# GSM5610940	Specimen_13
# GSM5610941	Specimen_15
# GSM5610942	Specimen_20


# 
# GSE198323 # URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE198323
# Run,Age,Assay Type,AvgSpotLen,Bases,BioProject,BioSample,Bytes,Center Name,Consent,DATASTORE filetype,DATASTORE provider,DATASTORE region,disease,Experiment,Instrument,Library Name,LibraryLayout,LibrarySelection,LibrarySource,Organism,Platform,ReleaseDate,Sample Name,source_name,SRA Study,tissue
# SRR18290865,88Yrs,OTHER,158,57001158192,PRJNA814627,SAMN26561906,25957392514,"DEPARTMENT OF NEUROSCIENCE, UPENN",public,"fastq,sra","gs,ncbi,s3","gs.US,ncbi.public,s3.us-east-1",Control,SRX14428578,NextSeq 550,GSM5945063,PAIRED,other,TRANSCRIPTOMIC,Homo sapiens,ILLUMINA,2022-04-19T00:00:00Z,GSM5945063,human hippocampus,SRP363403,hippocampus
# SRR18290866,85 Yrs,OTHER,160,23337002720,PRJNA814627,SAMN26561907,10963548076,"DEPARTMENT OF NEUROSCIENCE, UPENN",public,"fastq,sra","gs,ncbi,s3","gs.US,ncbi.public,s3.us-east-1",AD,SRX14428577,NextSeq 550,GSM5945062,PAIRED,other,TRANSCRIPTOMIC,Homo sapiens,ILLUMINA,2022-04-19T00:00:00Z,GSM5945062,human hippocampus,SRP363403,hippocampus
# SRR18290867,mixed sample (non-demultiplexed),OTHER,160,37557899680,PRJNA814627,SAMN26561908,18054973770,"DEPARTMENT OF NEUROSCIENCE, UPENN",public,"fastq,sra","gs,ncbi,s3","gs.US,ncbi.public,s3.us-east-1",mixed sample (non-demultiplexed),SRX14428576,NextSeq 550,GSM5945061,PAIRED,other,TRANSCRIPTOMIC,Homo sapiens,ILLUMINA,2022-04-19T00:00:00Z,GSM5945061,human hippocampus,SRP363403,hippocampus
# SRR18290869,mixed sample (non-demultiplexed),OTHER,160,49484406560,PRJNA814627,SAMN26561909,23624254078,"DEPARTMENT OF NEUROSCIENCE, UPENN",public,"fastq,sra","gs,ncbi,s3","gs.US,ncbi.public,s3.us-east-1",mixed sample (non-demultiplexed),SRX14428575,NextSeq 550,GSM5945060,PAIRED,other,TRANSCRIPTOMIC,Homo sapiens,ILLUMINA,2022-04-19T00:00:00Z,GSM5945060,hippocampus,SRP363403,hippocampus
# SRR18290871,mixed sample (non-demultiplexed),OTHER,160,47683785600,PRJNA814627,SAMN26561910,22823517250,"DEPARTMENT OF NEUROSCIENCE, UPENN",public,"fastq,sra","gs,ncbi,s3","gs.US,ncbi.public,s3.us-east-1",mixed sample (non-demultiplexed),SRX14428574,NextSeq 550,GSM5945059,PAIRED,other,TRANSCRIPTOMIC,Homo sapiens,ILLUMINA,2022-04-19T00:00:00Z,GSM5945059,human hippocampus,SRP363403,hippocampus
# SRR18290873,mixed sample (non-demultiplexed),OTHER,160,35731541440,PRJNA814627,SAMN26561911,16469146333,"DEPARTMENT OF NEUROSCIENCE, UPENN",public,"fastq,sra","gs,ncbi,s3","gs.US,ncbi.public,s3.us-east-1",mixed sample (non-demultiplexed),SRX14428573,NextSeq 550,GSM5945058,PAIRED,other,TRANSCRIPTOMIC,Homo sapiens,ILLUMINA,2022-04-19T00:00:00Z,GSM5945058,human hippocampus,SRP363403,hippocampus
# SRR18290868,mixed sample (non-demultiplexed),OTHER,160,52279646720,PRJNA814627,SAMN26561908,24533919404,"DEPARTMENT OF NEUROSCIENCE, UPENN",public,"fastq,sra","gs,ncbi,s3","gs.US,ncbi.public,s3.us-east-1",mixed sample (non-demultiplexed),SRX14428576,NextSeq 550,GSM5945061,PAIRED,other,TRANSCRIPTOMIC,Homo sapiens,ILLUMINA,2022-04-19T00:00:00Z,GSM5945061,human hippocampus,SRP363403,hippocampus
# SRR18290870,mixed sample (non-demultiplexed),OTHER,160,30625374400,PRJNA814627,SAMN26561909,14366154696,"DEPARTMENT OF NEUROSCIENCE, UPENN",public,"fastq,sra","gs,ncbi,s3","gs.US,ncbi.public,s3.us-east-1",mixed sample (non-demultiplexed),SRX14428575,NextSeq 550,GSM5945060,PAIRED,other,TRANSCRIPTOMIC,Homo sapiens,ILLUMINA,2022-04-19T00:00:00Z,GSM5945060,hippocampus,SRP363403,hippocampus
# SRR18290872,mixed sample (non-demultiplexed),OTHER,160,44775554080,PRJNA814627,SAMN26561910,20378526823,"DEPARTMENT OF NEUROSCIENCE, UPENN",public,"fastq,sra","gs,ncbi,s3","gs.US,ncbi.public,s3.us-east-1",mixed sample (non-demultiplexed),SRX14428574,NextSeq 550,GSM5945059,PAIRED,other,TRANSCRIPTOMIC,Homo sapiens,ILLUMINA,2022-04-19T00:00:00Z,GSM5945059,human hippocampus,SRP363403,hippocampus
# 
# GSM5945058	Specimen_36, _27, _26, _28
# GSM5945059	Specimen_29, _23, _32
# GSM5945060	Specimen_30, _37, _33, _31
# GSM5945061	Specimen_34, _24, _25
# GSM5945062	Specimen_35
# GSM5945063	Specimen_38
# 

setwd("/home/acampbell/PavLabEngrams/EngramCellClassifier")

library(dplyr)

# barcodes are not needed 

merged12 <- read.table("~/test_datasets/Zhou2022/GSE185277/GSM5609934_GEO_humanHippocampus_merged12.txt.gz", header = TRUE)
rownames(merged12) <- merged12$GENE
sample7plate1 <- merged12[,c(2:10001)]


sample2 <- read.table("~/test_datasets/Zhou2022/GSE185553/GSM5618239_Sample2.txt.gz", header = TRUE)
fun <- function(x){
  return( paste(x, "sample2", sep ="_") )
}
colnames(sample2) <- sapply(colnames(sample2),  fun, USE.NAMES = FALSE )
cellid <- colnames(sample2)

sample3 <- read.table("~/test_datasets/Zhou2022/GSE185553/GSM5618241_Sample3.txt.gz", header = TRUE)
fun <- function(x){
  return( paste(x, "sample3", sep ="_") )
}
colnames(sample3) <- sapply(colnames(sample3),  fun, USE.NAMES = FALSE )
cellid <- c(cellid, colnames(sample3) )

sample4 <- read.table("~/test_datasets/Zhou2022/GSE185553/GSM5618241_Sample4.txt.gz", header = TRUE)
fun <- function(x){
  return( paste(x, "sample4", sep ="_") )
}
colnames(sample4) <- sapply(colnames(sample4),  fun, USE.NAMES = FALSE )

sample5 <- read.table("~/test_datasets/Zhou2022/GSE185553/GSM5618239_Sample5.txt.gz", header = TRUE)
fun <- function(x){
  return( paste(x, "sample5", sep ="_") )
}
colnames(sample5) <- sapply(colnames(sample5),  fun, USE.NAMES = FALSE )


sample6lib1 <- read.table("~/test_datasets/Zhou2022/GSE185553/GSM5618238_Sample6_B1_lib1.txt.gz", header = TRUE)
fun <- function(x){
  return( paste(x, "sample6lib1", sep ="_") )
}
colnames(sample6lib1) <- sapply(colnames(sample6lib1),  fun, USE.NAMES = FALSE )


sample6lib2 <- read.table("~/test_datasets/Zhou2022/GSE185553/GSM5618238_Sample6-_B1_lib2.txt.gz", header = TRUE)
fun <- function(x){
  return( paste(x, "sample6lib2", sep ="_") )
}
colnames(sample6lib2) <- sapply(colnames(sample6lib2),  fun, USE.NAMES = FALSE )


sample7plate1 <- read.table("~/test_datasets/Zhou2022/GSE185553/GSM5618237_GEO_sample7_4yrs_plate1_seq12BC80_S5.deg.txt.gz", header = TRUE)
rownames(sample7plate1) <- sample7plate1$GENE
sample7plate1 <- sample7plate1[,c(2:10001)]

fun <- function(x){
  return( paste(x, "sample7plate1", sep ="_") )
}
colnames(sample7plate1) <- sapply(colnames(sample7plate1),  fun, USE.NAMES = FALSE )


sample7plate2 <- read.table("~/test_datasets/Zhou2022/GSE185553/GSM5618237_GEO_sample7_4yrs_plate2_seq12BC81_S6.deg.txt.gz", header = TRUE)
rownames(sample7plate2) <- sample7plate2$GENE
sample7plate2 <- sample7plate2[,c(2:10001)]
fun <- function(x){
  return( paste(x, "sample7plate2", sep ="_") )
}
colnames(sample7plate2) <- sapply(colnames(sample7plate2),  fun, USE.NAMES = FALSE )

sample8 <- read.table("~/test_datasets/Zhou2022/GSE185553/GSM5618241_Sample8.txt.gz", header = TRUE)
fun <- function(x){
  return( paste(x, "sample8", sep ="_") )
}
colnames(sample8) <- sapply(colnames(sample8),  fun, USE.NAMES = FALSE )


sample9 <- read.table("~/test_datasets/Zhou2022/GSE185277/GSM5610939_GEO_Sample9_seq1BC79_S4.txt.gz", header = TRUE)
rownames(sample9) <- sample9$GENE
sample9 <- sample9[,c(2:15001)]
fun <- function(x){
  return( paste(x, "sample9", sep ="_") )
}
colnames(sample9) <- sapply(colnames(sample9),  fun, USE.NAMES = FALSE )

sample10 <- read.table("~/test_datasets/Zhou2022/GSE185553/GSM5618239_Sample10.txt.gz", header = TRUE)
fun <- function(x){
  return( paste(x, "sample10", sep ="_") )
}
colnames(sample10) <- sapply(colnames(sample10),  fun, USE.NAMES = FALSE )


sample11lib1<- read.table("~/test_datasets/Zhou2022/GSE185553/GSM5618238_Sample11_B1_lib1.txt.gz", header = TRUE)                          
fun <- function(x){
  return( paste(x, "sample11lib1", sep ="_") )
}
colnames(sample11lib1) <- sapply(colnames(sample11lib1),  fun, USE.NAMES = FALSE )


sample11lib2<- read.table("~/test_datasets/Zhou2022/GSE185553/GSM5618238_Sample11_B1_lib2.txt.gz", header = TRUE)
fun <- function(x){
  return( paste(x, "sample11lib2", sep ="_") )
}
colnames(sample11lib2) <- sapply(colnames(sample11lib2),  fun, USE.NAMES = FALSE )

sample12 <- read.table("~/test_datasets/Zhou2022/GSE185553/GSM5618240_Sample12.txt.gz", header = TRUE)
fun <- function(x){
  return( paste(x, "sample12", sep ="_") )
}
colnames(sample12) <- sapply(colnames(sample12),  fun, USE.NAMES = FALSE )

sample13 <- read.table("~/test_datasets/Zhou2022/GSE185277/GSM5610940_GEO_sample13_HCT15HBW_seq12BC82_S7.deg.txt.gz",header = TRUE)
rownames(sample13) <- sample13$GENE
sample13 <- sample13[,c(2:25001)]
fun <- function(x){
  return( paste(x, "sample13", sep ="_") )
}
colnames(sample13) <- sapply(colnames(sample13),  fun, USE.NAMES = FALSE )

sample14 <- read.table("~/test_datasets/Zhou2022/GSE185553/GSM5618239_Sample14.txt.gz", header = TRUE)
fun <- function(x){
  return( paste(x, "sample14", sep ="_") )
}
colnames(sample14) <- sapply(colnames(sample14),  fun, USE.NAMES = FALSE )

sample15 <- read.table("~/test_datasets/Zhou2022/GSE185277/GSM5610941_sample15_BC80_S5.deg.txt.gz", header = TRUE)
rownames(sample15) <- sample15$GENE
sample15 <- sample15[,c(2:20001)]
fun <- function(x){
  return( paste(x, "sample15", sep ="_") )
}
colnames(sample15) <- sapply(colnames(sample15),  fun, USE.NAMES = FALSE )

sample16lib1 <- read.table("~/test_datasets/Zhou2022/GSE185553/GSM5618238_Sample16_B1_lib1.txt.gz"  , header = TRUE)                        
fun <- function(x){
  return( paste(x, "sample16lib1", sep ="_") )
}
colnames(sample16lib1) <- sapply(colnames(sample16lib1),  fun, USE.NAMES = FALSE )

sample16lib2 <- read.table("~/test_datasets/Zhou2022/GSE185553/GSM5618238_Sample16_B1_lib2.txt.gz", header = TRUE)
fun <- function(x){
  return( paste(x, "sample16lib2", sep ="_") )
}
colnames(sample16lib2) <- sapply(colnames(sample16lib2),  fun, USE.NAMES = FALSE )

sample17 <- read.table("~/test_datasets/Zhou2022/GSE185553/GSM5618240_Sample17.txt.gz", header = TRUE)
fun <- function(x){
  return( paste(x, "sample17", sep ="_") )
}
colnames(sample17) <- sapply(colnames(sample17),  fun, USE.NAMES = FALSE )

sample18 <- read.table("~/test_datasets/Zhou2022/GSE185553/GSM5618240_Sample18.txt.gz", header = TRUE)
fun <- function(x){
  return( paste(x, "sample18", sep ="_") )
}
colnames(sample18) <- sapply(colnames(sample18),  fun, USE.NAMES = FALSE )

sample19 <- read.table("~/test_datasets/Zhou2022/GSE185553/GSM5618241_Sample19.txt.gz", header = TRUE)
fun <- function(x){
  return( paste(x, "sample19", sep ="_") )
}
colnames(sample19) <- sapply(colnames(sample19),  fun, USE.NAMES = FALSE )

sample20lib1 <- read.table("~/test_datasets/Zhou2022/GSE185277/GSM5610942_Sample20_1.txt.gz", header = TRUE)
fun <- function(x){
  return( paste(x, "sample20lib1", sep ="_") )
}
colnames(sample20lib1) <- sapply(colnames(sample20lib1),  fun, USE.NAMES = FALSE )

sample20lib2 <- read.table("~/test_datasets/Zhou2022/GSE185277/GSM5610942_sample20_2.txt.gz", header = TRUE)
fun <- function(x){
  return( paste(x, "sample20lib2", sep ="_") )
}
colnames(sample20lib2) <- sapply(colnames(sample20lib2),  fun, USE.NAMES = FALSE )

sample21 <- read.table("~/test_datasets/Zhou2022/GSE185553/GSM5618239_Sample21.txt.gz", header = TRUE)
fun <- function(x){
  return( paste(x, "sample21", sep ="_") )
}
colnames(sample21) <- sapply(colnames(sample21),  fun, USE.NAMES = FALSE )

sample22lib1 <- read.table("~/test_datasets/Zhou2022/GSE185553/GSM5618238_Sample22_B1_lib1.txt.gz", header = TRUE)
fun <- function(x){
  return( paste(x, "sample22lib1", sep ="_") )
}
colnames(sample22lib1) <- sapply(colnames(sample22lib1),  fun, USE.NAMES = FALSE )

sample22lib2 <- read.table("~/test_datasets/Zhou2022/GSE185553/GSM5618238_Sample22_B1_lib2.txt.gz", header = TRUE)
fun <- function(x){
  return( paste(x, "sample22lib2", sep ="_") )
}
colnames(sample22lib2) <- sapply(colnames(sample22lib2),  fun, USE.NAMES = FALSE )

sample39lib1 <- read.table("~/test_datasets/Zhou2022/GSE185553/GSM5618238_Sample39_B1_lib1.txt.gz" , header = TRUE)
fun <- function(x){
  return( paste(x, "sample39lib1", sep ="_") )
}
colnames(sample39lib1) <- sapply(colnames(sample39lib1),  fun, USE.NAMES = FALSE )

sample39lib2 <- read.table("~/test_datasets/Zhou2022/GSE185553/GSM5618238_Sample39-lib2.txt.gz", header = TRUE)
fun <- function(x){
  return( paste(x, "sample39lib2", sep ="_") )
}
colnames(sample39lib2) <- sapply(colnames(sample39lib2),  fun, USE.NAMES = FALSE )

sample40 <- read.table("~/test_datasets/Zhou2022/GSE185553/GSM5618240_Sample40.txt.gz", header = TRUE)
fun <- function(x){
  return( paste(x, "sample40", sep ="_") )
}
colnames(sample40) <- sapply(colnames(sample40),  fun, USE.NAMES = FALSE )

sample41 <- read.table("~/test_datasets/Zhou2022/GSE185553/GSM5618240_Sample41.txt.gz", header = TRUE)
fun <- function(x){
  return( paste(x, "sample41", sep ="_") )
}
colnames(sample41) <- sapply(colnames(sample41),  fun, USE.NAMES = FALSE )

sample42 <- read.table("~/test_datasets/Zhou2022/GSE185553/GSM5618241_Sample42.txt.gz", header = TRUE)
fun <- function(x){
  return( paste(x, "sample42", sep ="_") )
}
colnames(sample42) <- sapply(colnames(sample42), fun, USE.NAMES = FALSE )


cellid <- c(colnames(sample2), colnames(sample3), colnames(sample4),
            colnames(sample5), colnames(sample6lib1), colnames(sample6lib2), colnames(sample7plate1),
            colnames(sample7plate2), colnames(sample8), colnames(sample9),
            colnames(sample10), colnames(sample11lib1), colnames(sample11lib2),
            colnames(sample12), colnames(sample13), colnames(sample14),
            colnames(sample15), colnames(sample16lib1), colnames(sample16lib2),
            colnames(sample17), colnames(sample18), colnames(sample19), 
            colnames(sample20lib1), colnames(sample20lib2), colnames(sample21),
            colnames(sample22lib1), colnames(sample39lib1), colnames(sample39lib2),
            colnames(sample40), colnames(sample41) )

multi.intersect <- function(x) Reduce(intersect, x) #takes lists of lists, c() will not work

shared.genes <- multi.intersect(list(rownames(sample2), rownames(sample3), rownames(sample4),
                                     rownames(sample5), rownames(sample6lib1), rownames(sample6lib2), rownames(sample7plate1),
                                     rownames(sample7plate2), rownames(sample8), rownames(sample9),
                                     rownames(sample10), rownames(sample11lib1), rownames(sample11lib2),
                                     rownames(sample12), rownames(sample13), rownames(sample14),
                                     rownames(sample15), rownames(sample16lib1), rownames(sample16lib2),
                                     rownames(sample17), rownames(sample18), rownames(sample19), 
                                     rownames(sample20lib1), rownames(sample20lib2), rownames(sample21),
                                     rownames(sample22lib1), rownames(sample39lib1), rownames(sample39lib2),
                                     rownames(sample40), rownames(sample41) )#closing list 
                                )#closing multi.intersect
length(shared.genes)

                                  
