library(adegenet)
library(hierfstat)
library(pegas)
library(ade4)
library(ggplot2)
library(vegan)
library(gplots)
library(factoextra)
library(FactoMineR)
library(mmod)
library(poppr)
library(tidyverse)


#All isolates
#Read in data.frame and save only rows with marker data in locus
Mydata <- read.table("Data/matrix_all_new_2filtermaf.txt", na.strings="NA", header=T)
dim(Mydata)
head(Mydata)[1:10,1:14]
locus<-Mydata[,-c(1:12)]



#Turn locus into Genind object
ind <- as.character(Mydata$Isolate) 
population <- as.character(Mydata$form)
Mydata1 <- df2genind(locus, ploidy = 1, ind.names = ind, pop = population, sep = "")
Mydata1
nAll(Mydata1)
div <- summary(Mydata1)
div
names(div)



#Make Genclone object 
gc<-as.genclone(Mydata1)

#Assess how much power you have to discriminate between unique individuals (shouldn't be an issue with 4000 SNPs
#genotype_curve(gc, sample=1000)
#runs forever, no need

#First overview: missing data, rare alleles and overall quality
overview<-locus_table(gc)
overview<-as.data.frame(overview)
write.table(overview, "Results/Popgen/locus_table_no_strata.txt")


#Replace "unknown" with NA
Mydata <- sapply(Mydata,function(x) {x <- gsub("unknown",NA,x)})
Mydata<-as.data.frame(Mydata)
head(Mydata)[1:10, 1:14]
class(Mydata)

#Make a new genind object
ind <- as.character(Mydata$Isolate) 
population <- as.character(Mydata$form)
Mydata1 <- df2genind(locus, ploidy = 1, ind.names = ind, pop = population, sep = "")
Mydata1
nAll(Mydata1)
div <- summary(Mydata1)
div
names(div)


#Add strata
Cul<-Mydata$Cultivar
Fyl<-Mydata$Fylke
Yea<-Mydata$year
Cou<-Mydata$Country
Form<-Mydata$form
stra<-data.frame(Form)
strata(Mydata1)<-stra
nameStrata(Mydata1)<- ~Form

#Make a new genclone
gc<-as.genclone(Mydata1)

clonecorr <- clonecorrect(gc, strata = ~Form)
clonecorr
#All individuals are MLGs

informloci(gc)
#No monomorphic loci


