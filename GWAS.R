
# Install library multtest - Run once
# source("http://www.bioconductor.org/biocLite.R")
# if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("multtest")
# biocLite("multtest")
# 
# # Install libraries gplots, LDheatmap, and genetics - Run once
# install.packages("gplots", dependencies=T)
# install.packages("LDheatmap", dependencies=T)
# install.packages("genetics", dependencies=T)

# Import libraries
library(multtest)
library("gplots")
library("LDheatmap")
library("genetics")
library("compiler") #this library is already installed in R
library(scatterplot3d)
library(BiocManager)
library(ape)
library(EMMREML)


# Install GAPIT package
source("https://raw.githubusercontent.com/jiabowang/GAPIT/refs/heads/master/gapit_functions.txt", encoding = "UTF-8")

##GWAS 
source("https://raw.githubusercontent.com/rwonneberger/R_functions/main/Ronja_functions.R")



myY<-read.table("Data/Isolate_info.txt", header=T, sep="\t", na.strings="unknown")
myY<-myY %>%
  filter(Cultivar=="Helium"|Cultivar=="Tiril", form=="N", Country=="Norway")
dim(myY)
head(myY)
myY$no<-NULL
myY$number<-NULL
myY<-myY[, c(1,4)]
myY$Cultivar<-as.character(myY$Cultivar)
myY$Cultivar[myY$Cultivar=="Helium"] <- "1"
myY$Cultivar[myY$Cultivar=="Tiril"] <- "2"
myY$Cultivar<-as.numeric(myY$Cultivar)

myY<-read.table("Data/Isolate_info.txt", header=T, sep="\t", na.strings="unknown")
myY<-myY %>%
  filter( form=="N", Country=="Norway", !is.na(rowtype))
dim(myY)
head(myY)
myY$no<-NULL
myY$number<-NULL
myY<-myY[, c(1,12)]



#Step 1: Set working directory and import data

#All isolates
myG <- read.table("Data/Hmp_Hapmap_2maffilter_2025.hmp.txt", na.strings="NA", head=F)
#myG <- read.table("BARN_50K_filtered_KNNimputed_TAGS_R20.95_Morex1.hmp.txt", na.strings="NA", head=F)

myG[1,-c(1:11)] <-gsub("^\\bX", "", myG[1,-c(1:11)] )


length(myY$Isolate)
dim(myG)[2]

#myY%<>%arrange(GEN)

# Check if both files contain the same individuals, and in the same order
myY$Isolate == myG[1,-c(1:11)] 
myY$Isolate[which(!myY$Isolate%in%myG[1,-c(1:11)])] 

# keep only individuals that are present in both files
keep_ind<-unlist(intersect(myY$Isolate, myG[1,-c(1:11)]))

# which ones are not in both files?
base::setdiff(myY$Isolate, myG[1,-c(1:11)])
base::setdiff(myG[1,-c(1:11)], myY$Isolate)

myY<-myY[myY$Isolate%in%keep_ind,]

names(myG)<-myG[1,] # Set individual names as column names to facilitate filtering
myG%<>%dplyr::select(c(1:11), all_of(keep_ind))

myG$chrom[myG$chrom == 'CM016795_1'] <- '1'
myG$chrom[myG$chrom == 'CM016796_1'] <- '2'
myG$chrom[myG$chrom == 'CM016797_1'] <- '3'
myG$chrom[myG$chrom == 'CM016798_1'] <- '4'
myG$chrom[myG$chrom == 'CM016799_1'] <- '5'
myG$chrom[myG$chrom == 'CM016800_1'] <- '6'
myG$chrom[myG$chrom == 'CM016801_1'] <- '7'
myG$chrom[myG$chrom == 'CM016802_1'] <- '8'
myG$chrom[myG$chrom == 'CM016803_1'] <- '9'
myG$chrom[myG$chrom == 'CM016804_1'] <- '10'
myG$chrom[myG$chrom == 'CM016805_1'] <- '11'
myG$chrom[myG$chrom == 'CM016806_1'] <- '12'

#Check again if both files contain the same individuals
table(myY$Isolate == myG[1,-c(1:11)])


# Run GWAS
if (!dir.exists(paste0("Results/rowtype_PC1"))){
  dir.create(paste0("Results/rowtype_PC1"))
}

setwd(paste0("Results/rowtype_PC1"))

for(i in c( "MLMM", "BLINK", "MLM", "CMLM", "FarmCPU")){
  res <- try(myGAPIT <- GAPIT(
    Y=myY,
    G=myG,
    SNP.MAF=0.01,
    model=i, 
    PCA.total=1,
    Major.allele.zero = TRUE
  ))
  if(inherits(res, "try-error"))
  {
    #error handling code, maybe just skip this iteration using
    next
  }
  #rest of iteration for case of no error
}



############ Make Manhattan plot

mlmm0<-fread("Results/Hel_Tir_N/GAPIT.Association.GWAS_Results.MLMM.Cultivar.csv")[,1:4]
blink0<-fread("Results/Hel_Tir_N/GAPIT.Association.GWAS_Results.BLINK.Cultivar.csv")[,1:4]
mlmm1<-fread("Results/Hel_Tir_N_PC1/GAPIT.Association.GWAS_Results.MLMM.Cultivar.csv")[,1:4]
blink1<-fread("Results/Hel_Tir_N_PC1/GAPIT.Association.GWAS_Results.BLINK.Cultivar.csv")[,1:4]
mlmm0$Model<-"MLMM 0 PCs"
blink0$Model<-"BLINK 0 PCs"
mlmm1$Model<-"MLMM 1 PC"
blink1$Model<-"BLINK 1 PC"


df1<-rbind(mlmm0, blink0, mlmm1, blink1)
df1%<>%filter(Chr%in%c(1,2,3,4,5,6,7,8,9,10,11,12))
df1$Chr<-as.numeric(df1$Chr)
df1$Model<-factor(df1$Model, levels=c("MLMM 0 PCs", "BLINK 0 PCs", "MLMM 1 PC", "BLINK 1 PC"))
names(df1)

thresh<--log10(0.05/nrow(mlmm0))
scaleFUN <- function(x) sprintf("%.2f", x)

tiff("Figs/GWAS_Results.tiff", width=20, height= 15, res=600, units="cm")
ggplot(df1,aes(x=as.numeric(Pos)/1000,y=-log10(P.value)))  + geom_point(size=0.7)+   theme_classic() + theme(panel.grid.minor = element_line(colour = "grey90"), panel.grid.major = element_line(colour = "grey90"), panel.border = element_rect(fill=NA))+ 
  scale_x_continuous( labels = function(x) format(x, scientific = FALSE))+
  theme(axis.text.x = element_text(angle = 90, vjust=.5)) + 
  facet_grid(Model~Chr, scales='free_x')+ xlab("Position (kbp)")+ ylab("-log10(p)")  +scale_y_continuous(labels=scaleFUN) + geom_hline(yintercept =thresh, color="red")
dev.off()
