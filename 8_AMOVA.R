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
library(gridExtra)

source("https://raw.githubusercontent.com/rwonneberger/R_functions/main/Ronja_functions.R")

#### Fylke

Mydata <- read.table("Data/matrix_all_new_2filtermaf.txt", na.strings="NA", header=T)

#Only keep Norwegian net form isolates with known cultivar. Only keep cultivars that have a count of at least 3
Mydata<-Mydata %>%
  filter(form=="N", Country=="Norway", Fylke!="unknown") %>%
  group_by(Fylke) %>%
  filter(n()>2)

locus<-Mydata[,-c(1:12)]

#Replace "unknown" with NA
Mydata <- sapply(Mydata,function(x) {x <- gsub("unknown",NA,x)})
Mydata<-as.data.frame(Mydata)

ind <- as.character(Mydata$Isolate) 
population <- as.character(Mydata$Fylke)
Mydata1 <- df2genind(locus, ploidy = 1, ind.names = ind, pop = population, sep = "")

#Add strata
Fyl<-Mydata$Fylke
stra<-data.frame(Fyl)
strata(Mydata1)<-stra
nameStrata(Mydata1)<- ~Fylke

MyGC<-as.genclone(Mydata1)

#By fylke
table(strata(MyGC, ~Fylke, combine = FALSE))
PTamova <- poppr.amova(Mydata1, ~Fylke,  missing="mean")
PTamova
set.seed(1999)
a <- randtest(PTamova, nrepet = 999)
a








#### Cultivar
Mydata <- read.table("Data/matrix_all_new_2filtermaf.txt", na.strings="NA", header=T)

#Only keep Norwegian net form isolates with known cultivar. Only keep cultivars that have a count of at least 3
Mydata<-Mydata %>%
  filter(form=="N", Country=="Norway", Cultivar!="unknown") %>%
  group_by(Cultivar) %>%
  filter(n()>2)

locus<-Mydata[,-c(1:12)]

#Replace "unknown" with NA
Mydata <- sapply(Mydata,function(x) {x <- gsub("unknown",NA,x)})
Mydata<-as.data.frame(Mydata)

ind <- as.character(Mydata$Isolate) 
population <- as.character(Mydata$Cultivar)
Mydata1 <- df2genind(locus, ploidy = 1, ind.names = ind, pop = population, sep = "")

#Add strata
Cul<-Mydata$Cultivar
stra<-data.frame(Cul)
strata(Mydata1)<-stra
nameStrata(Mydata1)<- ~Cultivar

MyGC<-as.genclone(Mydata1)

#By fylke
table(strata(MyGC, ~Cultivar, combine = FALSE))
PTamova <- poppr.amova(Mydata1, ~Cultivar,  missing="mean")
PTamova
set.seed(1999)
a <- randtest(PTamova, nrepet = 999)
a







#### Year
Mydata <- read.table("Data/matrix_all_new_2filtermaf.txt", na.strings="NA", header=T)

#Only keep Norwegian net form isolates with known cultivar. Only keep cultivars that have a count of at least 3
Mydata<-Mydata %>%
  filter(form=="N", Country=="Norway", year!="unknown") %>%
  group_by(year) %>%
  filter(n()>2)

locus<-Mydata[,-c(1:12)]

#Replace "unknown" with NA
Mydata <- sapply(Mydata,function(x) {x <- gsub("unknown",NA,x)})
Mydata<-as.data.frame(Mydata)

ind <- as.character(Mydata$Isolate) 
population <- as.character(Mydata$year)
Mydata1 <- df2genind(locus, ploidy = 1, ind.names = ind, pop = population, sep = "")

#Add strata
year<-Mydata$year
stra<-data.frame(year)
strata(Mydata1)<-stra
nameStrata(Mydata1)<- ~year

MyGC<-as.genclone(Mydata1)

#By fylke
table(strata(MyGC, ~year, combine = FALSE))
PTamova <- poppr.amova(Mydata1, ~year,  missing="mean")
PTamova
set.seed(1999)
a <- randtest(PTamova, nrepet = 999)
a







#### MAT
Mydata <- read.table("Data/matrix_all_new_2filtermaf.txt", na.strings="NA", header=T)

Mydata<-Mydata %>%
  filter(form=="N", Country=="Norway")

locus<-Mydata[,-c(1:12)]

#Replace "unknown" with NA
Mydata <- sapply(Mydata,function(x) {x <- gsub("unknown",NA,x)})
Mydata<-as.data.frame(Mydata)

ind <- as.character(Mydata$Isolate) 
population <- as.character(Mydata$MT)
Mydata1 <- df2genind(locus, ploidy = 1, ind.names = ind, pop = population, sep = "")

#Add strata
MT<-Mydata$MT
stra<-data.frame(MT)
strata(Mydata1)<-stra
nameStrata(Mydata1)<- ~MT

MyGC<-as.genclone(Mydata1)

#By fylke
table(strata(MyGC, ~MT, combine = FALSE))
PTamova <- poppr.amova(Mydata1, ~MT,  missing="mean")
PTamova
set.seed(1999)
a <- randtest(PTamova, nrepet = 999)
a






#### Form
Mydata <- read.table("Data/matrix_all_new_2filtermaf.txt", na.strings="NA", header=T)

Mydata<-Mydata %>%
  filter(Country=="Norway", form!="unknown")

locus<-Mydata[,-c(1:12)]

#Replace "unknown" with NA
Mydata <- sapply(Mydata,function(x) {x <- gsub("unknown",NA,x)})
Mydata<-as.data.frame(Mydata)

ind <- as.character(Mydata$Isolate) 
population <- as.character(Mydata$form)
Mydata1 <- df2genind(locus, ploidy = 1, ind.names = ind, pop = population, sep = "")

#Add strata
form<-Mydata$form
stra<-data.frame(form)
strata(Mydata1)<-stra
nameStrata(Mydata1)<- ~form

MyGC<-as.genclone(Mydata1)

#By fylke
table(strata(MyGC, ~form, combine = FALSE))
PTamova <- poppr.amova(Mydata1, ~form,  missing="mean")
PTamova
set.seed(1999)
a <- randtest(PTamova, nrepet = 999)
a









#### Year 1995 vs rest
Mydata <- read.table("Data/matrix_all_new_2filtermaf.txt", na.strings="NA", header=T)

#Only keep Norwegian net form isolates with known cultivar. Only keep cultivars that have a count of at least 3
Mydata<-Mydata %>%
  filter(form=="N", Country=="Norway", year!="unknown") %>%
  group_by(year) %>%
  filter(n()>2)

Mydata$year <- gsub("^2.*", "recent", Mydata$year)


locus<-Mydata[,-c(1:12)]

#Replace "unknown" with NA
Mydata <- sapply(Mydata,function(x) {x <- gsub("unknown",NA,x)})
Mydata<-as.data.frame(Mydata)

ind <- as.character(Mydata$Isolate) 
population <- as.character(Mydata$year)
Mydata1 <- df2genind(locus, ploidy = 1, ind.names = ind, pop = population, sep = "")

#Add strata
year<-Mydata$year
stra<-data.frame(year)
strata(Mydata1)<-stra
nameStrata(Mydata1)<- ~year

MyGC<-as.genclone(Mydata1)

#By fylke
table(strata(MyGC, ~year, combine = FALSE))
PTamova <- poppr.amova(Mydata1, ~year,  missing="mean")
PTamova
set.seed(1999)
a <- randtest(PTamova, nrepet = 999)
a





#### Only post 1995
Mydata <- read.table("Data/matrix_all_new_2filtermaf.txt", na.strings="NA", header=T)

#Only keep Norwegian net form isolates with known cultivar. Only keep cultivars that have a count of at least 3
Mydata<-Mydata %>%
  filter(form=="N", Country=="Norway", year!="unknown", year!="1995") %>%
  group_by(year) %>%
  filter(n()>2)

locus<-Mydata[,-c(1:12)]

#Replace "unknown" with NA
Mydata <- sapply(Mydata,function(x) {x <- gsub("unknown",NA,x)})
Mydata<-as.data.frame(Mydata)

ind <- as.character(Mydata$Isolate) 
population <- as.character(Mydata$year)
Mydata1 <- df2genind(locus, ploidy = 1, ind.names = ind, pop = population, sep = "")

#Add strata
year<-Mydata$year
stra<-data.frame(year)
strata(Mydata1)<-stra
nameStrata(Mydata1)<- ~year

MyGC<-as.genclone(Mydata1)

#By fylke
table(strata(MyGC, ~year, combine = FALSE))
PTamova <- poppr.amova(Mydata1, ~year,  missing="mean")
PTamova
set.seed(1999)
a <- randtest(PTamova, nrepet = 999)
a







#### Cultivar i fylke
Mydata <- read.table("Data/matrix_all_new_2filtermaf.txt", na.strings="NA", header=T)

Mydata<-Mydata %>%
  filter(form=="N", Country=="Norway", Cultivar!="unknown", Fylke!="unknown") %>%
  group_by(Cultivar, Fylke) %>%
  filter(n()>2)

locus<-Mydata[,-c(1:12)]

#Replace "unknown" with NA
Mydata <- sapply(Mydata,function(x) {x <- gsub("unknown",NA,x)})
Mydata<-as.data.frame(Mydata)

ind <- as.character(Mydata$Isolate) 
population <- as.character(Mydata$Cultivar)
Mydata1 <- df2genind(locus, ploidy = 1, ind.names = ind, pop = population, sep = "")

#Add strata
Cul<-Mydata$Cultivar
Fylke<-Mydata$Fylke
stra<-data.frame(Fylke, Cul)
strata(Mydata1)<-stra
nameStrata(Mydata1)<- ~Fylke/Cultivar

MyGC<-as.genclone(Mydata1)

#By fylke
table(strata(MyGC, ~Fylke/Cultivar, combine = FALSE))
PTamova <- poppr.amova(Mydata1, ~Fylke/Cultivar,  missing="mean")
PTamova
set.seed(1999)
a <- randtest(PTamova, nrepet = 999)
a



#### Cultivar i fylke, only Hel Tir
Mydata <- read.table("Data/matrix_all_new_2filtermaf.txt", na.strings="NA", header=T)

Mydata<-Mydata %>%
  filter(form=="N", Country=="Norway", Cultivar=="Helium"|Cultivar=="Tiril", Fylke!="unknown") %>%
  group_by(Cultivar, Fylke) %>%
  filter(n()>2)

locus<-Mydata[,-c(1:12)]

#Replace "unknown" with NA
Mydata <- sapply(Mydata,function(x) {x <- gsub("unknown",NA,x)})
Mydata<-as.data.frame(Mydata)

ind <- as.character(Mydata$Isolate) 
population <- as.character(Mydata$Cultivar)
Mydata1 <- df2genind(locus, ploidy = 1, ind.names = ind, pop = population, sep = "")

#Add strata
Cul<-Mydata$Cultivar
Fylke<-Mydata$Fylke
stra<-data.frame(Fylke, Cul)
strata(Mydata1)<-stra
nameStrata(Mydata1)<- ~Fylke/Cultivar

MyGC<-as.genclone(Mydata1)

#By fylke
table(strata(MyGC, ~Fylke/Cultivar, combine = FALSE))
PTamova <- poppr.amova(Mydata1, ~Fylke/Cultivar,  missing="mean")
PTamova
set.seed(1999)
a <- randtest(PTamova, nrepet = 999)
a





#### Cultivar, only Hel Tir
Mydata <- read.table("Data/matrix_all_new_2filtermaf.txt", na.strings="NA", header=T)

Mydata<-Mydata %>%
  filter(form=="N", Country=="Norway", Cultivar=="Helium"|Cultivar=="Tiril") %>%
  group_by(Cultivar) %>%
  filter(n()>2)

locus<-Mydata[,-c(1:12)]

#Replace "unknown" with NA
Mydata <- sapply(Mydata,function(x) {x <- gsub("unknown",NA,x)})
Mydata<-as.data.frame(Mydata)

ind <- as.character(Mydata$Isolate) 
population <- as.character(Mydata$Cultivar)
Mydata1 <- df2genind(locus, ploidy = 1, ind.names = ind, pop = population, sep = "")

#Add strata
Cul<-Mydata$Cultivar
stra<-data.frame(Cul)
strata(Mydata1)<-stra
nameStrata(Mydata1)<- ~Cultivar

MyGC<-as.genclone(Mydata1)

#By fylke
table(strata(MyGC, ~Cultivar, combine = FALSE))
PTamova <- poppr.amova(Mydata1, ~Cultivar,  missing="mean")
PTamova
set.seed(1999)
a <- randtest(PTamova, nrepet = 999)
a







#### MAT i fylke
Mydata <- read.table("Data/matrix_all_new_2filtermaf.txt", na.strings="NA", header=T)

Mydata<-Mydata %>%
  filter(form=="N", Country=="Norway", Fylke!="unknown")

locus<-Mydata[,-c(1:12)]

#Replace "unknown" with NA
Mydata <- sapply(Mydata,function(x) {x <- gsub("unknown",NA,x)})
Mydata<-as.data.frame(Mydata)

ind <- as.character(Mydata$Isolate) 
population <- as.character(Mydata$MT)
Mydata1 <- df2genind(locus, ploidy = 1, ind.names = ind, pop = population, sep = "")

#Add strata
MT<-Mydata$MT
Fyl<-Mydata$Fylke
stra<-data.frame( Fyl,MT)
strata(Mydata1)<-stra
nameStrata(Mydata1)<- ~Fyl/MT

MyGC<-as.genclone(Mydata1)

#By fylke
table(strata(MyGC, ~Fyl/MT, combine = FALSE))
PTamova <- poppr.amova(Mydata1, ~Fyl/MT,  missing="mean")
PTamova
set.seed(1999)
a <- randtest(PTamova, nrepet = 999)
a




#### MAT i Cultivar
Mydata <- read.table("Data/matrix_all_new_2filtermaf.txt", na.strings="NA", header=T)

Mydata<-Mydata %>%
  filter(form=="N", Country=="Norway", Cultivar!="unknown") %>%
  group_by(MT, Cultivar) %>%
  filter(n()>2)

locus<-Mydata[,-c(1:12)]

#Replace "unknown" with NA
Mydata <- sapply(Mydata,function(x) {x <- gsub("unknown",NA,x)})
Mydata<-as.data.frame(Mydata)

ind <- as.character(Mydata$Isolate) 
population <- as.character(Mydata$MT)
Mydata1 <- df2genind(locus, ploidy = 1, ind.names = ind, pop = population, sep = "")

#Add strata
MAT<-Mydata$MT
Cul<-Mydata$Cultivar
stra<-data.frame(Cul, MAT)
strata(Mydata1)<-stra
nameStrata(Mydata1)<- ~Cultivar/MT

MyGC<-as.genclone(Mydata1)

#By fylke
table(strata(MyGC, ~Cultivar/MT, combine = FALSE))
PTamova <- poppr.amova(Mydata1, ~Cultivar/MT,  missing="mean")
PTamova
set.seed(1999)
a <- randtest(PTamova, nrepet = 999)
a





