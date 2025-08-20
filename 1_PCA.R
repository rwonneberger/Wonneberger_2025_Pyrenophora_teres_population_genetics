library(adegenet)
library(hierfstat)
library(pegas)
library(ade4)
library(vegan)
library(gplots)
library(factoextra)
library(FactoMineR)
library(mmod)
library(poppr)
library(gridExtra)


source("https://raw.githubusercontent.com/rwonneberger/R_functions/main/Ronja_functions.R")


colors<-c('#4363d8', '#f58231', '#e6194B', '#3cb44b', '#ffe119',  '#42d4f4', '#f032e6', '#fabed4', '#469990', '#dcbeff', '#9A6324', '#fffac8', '#800000', '#aaffc3', '#000075', '#a9a9a9', '#ffffff', '#000000')

pointsize<-1
plottitlesize<-10
legendtextsize<-8
legendtitlesize<-10

theme_set(theme_minimal(base_size = 10)) 

####Fig1ABC
#All isolates
Mydata <- read.table("Data/matrix_all_new_2filtermaf.txt", na.strings="NA", header=T)
dim(Mydata)
locus<-Mydata[,-c(1:12)]

ind <- as.character(Mydata$Isolate) 
population <- as.character(Mydata$Isolate)
Mydata1 <- df2genind(locus, ploidy = 1, ind.names = ind, pop = population, sep = "")
Mydata1
nAll(Mydata1)
div <- summary(Mydata1)
div
names(div)


X <- scaleGen(Mydata1, NA.method="mean")
class(X)

pca1 <- dudi.pca(X,scannf=FALSE, center=T, nf=3)


#Plot PCA of all:by form
a<-fviz_pca_ind(pca1, geom="point",invisible="quali",  habillage=Mydata$form, pointsize=pointsize) +
  labs(title="All isolates")+
  scale_color_manual(values = c('#4363d8',  '#e6194B', "black"))+
  theme_minimal()+
  scale_shape_manual(values=c(16, 16, 16))+ theme(plot.title = element_text(size=plottitlesize), legend.text=element_text(size=legendtextsize), legend.title=element_text(size=legendtitlesize), legend.key.height = unit(0.4, 'cm'))


#Plot PCA of all: by country
b<-fviz_pca_ind(pca1, geom="point", invisible="quali", habillage=Mydata$Country, pointsize=pointsize) +
  labs(title="All isolates")+
  scale_color_manual(values = c('#a9a9a9', '#f58231', '#e6194B', '#3cb44b', '#42d4f4','#ffe119',   '#f032e6'))+
  theme_minimal()+
  scale_shape_manual(values=c(16, 16, 16, 16, 16, 16, 16))+ theme(plot.title = element_text(size=plottitlesize), legend.text=element_text(size=legendtextsize), legend.title=element_text(size=legendtitlesize),legend.key.height = unit(0.4, 'cm'))



#All isolates
Mydata <- read.table("Data/matrix_all_new_2filtermaf.txt", na.strings="NA", header=T)
dim(Mydata)

Mydata<-Mydata %>%
  filter(form=="N", Isolate!="CAWB_05_Pt_4")

dim(Mydata)
locus<-Mydata[,-c(1:12)]

ind <- as.character(Mydata$Isolate) 
population <- as.character(Mydata$Isolate)
Mydata1 <- df2genind(locus, ploidy = 1, ind.names = ind, pop = population, sep = "")
Mydata1
nAll(Mydata1)
div <- summary(Mydata1)
div
names(div)


X <- scaleGen(Mydata1, NA.method="mean")
class(X)

pca1 <- dudi.pca(X,scannf=FALSE, center=T, nf=3)


c<-fviz_pca_ind(pca1, geom="point", invisible="quali", habillage=Mydata$Country, pointsize=pointsize) +
  labs(title="All NFNB")+
  scale_color_manual(values =  c( '#f58231', '#e6194B', '#3cb44b', '#42d4f4','#ffe119',   '#f032e6'))+
  theme_minimal()+
  scale_shape_manual(values=c(16, 16, 16, 16, 16, 16, 16))+ theme(plot.title = element_text(size=plottitlesize), legend.text=element_text(size=legendtextsize), legend.title=element_text(size=legendtitlesize),legend.key.height = unit(0.4, 'cm'))





##Norwegian net form isolates only by rowtype
#All isolates
Mydata <- read.table("Data/matrix_all_new_2filtermaf.txt", na.strings="NA", header=T)
dim(Mydata)

Mydata<-Mydata %>%
  filter(Country=="Norway", form=="N", Isolate!="117LII")

dim(Mydata)
locus<-Mydata[,-c(1:12)]

ind <- as.character(Mydata$Isolate) 
population <- as.character(Mydata$Isolate)
Mydata1 <- df2genind(locus, ploidy = 1, ind.names = ind, pop = population, sep = "")
Mydata1
nAll(Mydata1)
div <- summary(Mydata1)
div
names(div)


X <- scaleGen(Mydata1, NA.method="mean")
class(X)

pca1 <- dudi.pca(X,scannf=FALSE, center=T, nf=3)



#Plot PCA of all Norwegian:by form
d<-fviz_pca_ind(pca1,  geom="point", invisible="quali",habillage=Mydata$rowtype, pointsize=pointsize) +
  labs(title="Norwegian NFNB")+
  scale_color_manual(values = c('#4363d8',  '#e6194B', "black"))+
  theme_minimal()+
  scale_shape_manual(values=c(16, 16, 16, 16, 16))+ theme(plot.title = element_text(size=plottitlesize), legend.text=element_text(size=legendtextsize), legend.title=element_text(size=legendtitlesize),legend.key.height = unit(0.4, 'cm'))







##Norwegian Helium and Tiril only
#All isolates
Mydata <- read.table("Data/matrix_all_new_2filtermaf.txt", na.strings="NA", header=T)
dim(Mydata)

Mydata<-Mydata %>%
  filter(Country=="Norway", Cultivar=="Helium" | Cultivar=="Tiril")

dim(Mydata)
locus<-Mydata[,-c(1:12)]

ind <- as.character(Mydata$Isolate) 
population <- as.character(Mydata$Isolate)
Mydata1 <- df2genind(locus, ploidy = 1, ind.names = ind, pop = population, sep = "")
Mydata1
nAll(Mydata1)
div <- summary(Mydata1)
div
names(div)


X <- scaleGen(Mydata1, NA.method="mean")
class(X)

pca1 <- dudi.pca(X,scannf=FALSE, center=T, nf=3)

#Plot PCA of Norw netform: by year
e<-fviz_pca_ind(pca1,  geom="point",  invisible="quali", habillage=Mydata$MT, pointsize=pointsize) +
  labs(title="Helium and Tiril")+
  scale_color_manual(values = c('#4363d8', '#f58231', '#e6194B', '#f032e6','#3cb44b', '#ffe119',  '#42d4f4',  '#fabed4', '#469990', '#dcbeff', '#9A6324', '#fffac8', '#800000', '#aaffc3', '#000075', '#a9a9a9', '#ffffff', '#000000'))+
  theme_minimal()+
  scale_shape_manual(values=c(16, 16, 16, 16, 16, 16, 16, 16))+ theme(plot.title = element_text(size=plottitlesize), legend.text=element_text(size=legendtextsize),legend.title=element_text(size=legendtitlesize), legend.key.height = unit(0.4, 'cm'))


e






##Norwegian net form isolates only
#All isolates
Mydata <- read.table("Data/matrix_all_new_2filtermaf.txt", na.strings="NA", header=T)
dim(Mydata)

Mydata<-Mydata %>%
  filter(Country=="Norway", form=="N")

dim(Mydata)
locus<-Mydata[,-c(1:12)]

ind <- as.character(Mydata$Isolate) 
population <- as.character(Mydata$Isolate)
Mydata1 <- df2genind(locus, ploidy = 1, ind.names = ind, pop = population, sep = "")
Mydata1
nAll(Mydata1)
div <- summary(Mydata1)
div
names(div)


X <- scaleGen(Mydata1, NA.method="mean")
class(X)

pca1 <- dudi.pca(X,scannf=FALSE, center=T, nf=3)



#Plot PCA of Norw netform: by fylke
f<-fviz_pca_ind(pca1,  geom="point",  invisible="quali", habillage=Mydata$Fylke, pointsize=pointsize) +
  labs(title="Norwegian NFNB")+
  scale_color_manual(values = c('#4363d8', '#f58231', '#e6194B', '#f032e6','#3cb44b', '#ffe119',  '#42d4f4',  '#fabed4', '#469990', '#dcbeff', '#9A6324', '#fffac8', '#800000', '#aaffc3', '#000075', '#a9a9a9', '#ffffff', '#000000'))+
  theme_minimal()+
  scale_shape_manual(values=c(16, 16, 16, 16, 16, 16, 16, 16, 16, 16))+ theme(plot.title = element_text(size=plottitlesize), legend.text=element_text(size=legendtextsize), legend.title=element_text(size=legendtitlesize),legend.key.height = unit(0.4, 'cm'))



##Norwegian net form isolates only
#All isolates
Mydata <- read.table("Data/matrix_all_new_2filtermaf.txt", na.strings="NA", header=T)
dim(Mydata)

Mydata<-Mydata %>%
  filter(Country=="Norway", form=="N", Isolate!="117LII")

dim(Mydata)
locus<-Mydata[,-c(1:12)]

ind <- as.character(Mydata$Isolate) 
population <- as.character(Mydata$Isolate)
Mydata1 <- df2genind(locus, ploidy = 1, ind.names = ind, pop = population, sep = "")
Mydata1
nAll(Mydata1)
div <- summary(Mydata1)
div
names(div)


X <- scaleGen(Mydata1, NA.method="mean")
class(X)

pca1 <- dudi.pca(X,scannf=FALSE, center=T, nf=3)



#Plot PCA of all Norwegian:by form
g<-fviz_pca_ind(pca1,  geom="point", invisible="quali",habillage=Mydata$year, pointsize=pointsize) +
  labs(title="Norwegian NFNB")+
  scale_color_manual(values = c("#e6194B",   "#f58231", "#3cb44b", "#4363d8", "#42d4f4"))+
  theme_minimal()+
  scale_shape_manual(values=c(16, 16, 16, 16, 16))+ theme(plot.title = element_text(size=plottitlesize), legend.text=element_text(size=legendtextsize), legend.title=element_text(size=legendtitlesize),legend.key.height = unit(0.4, 'cm'))



tiff("Figs/Fig1all.tiff", unit="cm", res=600, width=17.4, height=23.4)
ggarrange(a, b, c, d, e, f, g, ncol=2, nrow=4, labels=c("a", "b", "c", "d", "e", "f", "g"))
dev.off()


