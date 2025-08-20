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
library(picante)
library(pals)





color.plot.phylo.custom<-function(
  # this function takes a phylogeny and a dataframe with the trait of interest and the taxa columns explicitly passed
  # The number of colors defaults to the number of factors (this cannot be chnaged safely) for factors and 12 for continous traits
  # A user supplied vector of colors can be passed to col.names, be sure this equals num.breaks
  # leg.title provides a title for the legend
  phylo, df, trait, taxa.names,
  num.breaks = ifelse(is.factor(df[,trait]),
                      length(levels(df[,trait])), 12),
  col.names = rainbow(ifelse(length(num.breaks) > 1, length(num.breaks) - 1, num.breaks)),
  cut.labs = NULL,
  leg.title = NULL,
  main = trait,
  leg.cex = 0.4,
  tip.labs = NULL,
  ...
)
{
  # Get the initial par settings so they can be restored
  init.par <- par(mar = c(0, 0, 1, 0))
  # some data input error checking, all taxa in tree and df
  # no missing data values
  stopifnot( trait %in% names(df), taxa.names %in% names(df),
             class(df) == "data.frame", class(phylo) == "phylo")
  len.tips <- length(phylo$tip.label)
  len.taxa <- length(df[,taxa.names])
  if (len.tips != len.taxa |
      sum(phylo$tip.label %in% df[,taxa.names]) != len.taxa) {
    stop("ERROR. Missing taxa in tree or data frame; # tips: ",
         len.tips, "# taxa: ", len.taxa, "# tips in df: ",
         sum(phylo$tip.label %in% df[,taxa.names]))
  }
  # ensure that the order of the data frame matches the tips
  order <- match(phylo$tip.label, df[,taxa.names])
  ordered.trait <- df[trait][order,]
  # cut up the trait and assign a list of colors
  if(is.factor(ordered.trait)){
    levs <- levels(ordered.trait)
    tip.color <- rep("black", times = len.taxa)
    tip.color <- col.names[match(ordered.trait, levs)]
  }else{
    tip.color = as.character(cut(
      ordered.trait,
      breaks = num.breaks,
      labels = col.names
    ))
    # levs gets used in legend, make one for continous data
    levs <- levels(cut(ordered.trait, breaks = num.breaks))
  }
  if(!is.null(tip.labs)) {
    phylo$tip.label <- df[tip.labs][order,]
  }
  plot.phylo(
    phylo,
    ## y.lim = c(0,80),
    cex = 0.5,
    tip.color = tip.color,
    main = main,
    ...
  )
  title(line = 0)
  
  if(is.null(cut.labs)) cut.labs <- levs
  legend("bottomleft", cut.labs,
         fill = col.names,
         inset = 0.0, title = leg.title,
         cex = leg.cex
  )
  # restore settings ?
  on.exit(par(init.par))
}

## Change inset = 0.0 to move legend


## Dendro all Norw N iso by Cul

        Mydata <- read.table("Data/matrix_all_new_2filtermaf.txt", na.strings="NA", header=T)

#Replace "unknown" with NA
Mydata <- sapply(Mydata,function(x) {x <- gsub("unknown",NA,x)})
Mydata<-as.data.frame(Mydata)
names(Mydata)[names(Mydata) == "Isolate"] <- "taxa"

Mydata<-Mydata %>%
  filter(form=="N", Country=="Norway") 

locus<-Mydata[,-c(1:12)]

ind <- as.character(Mydata$taxa) 
population <- as.character(Mydata$taxa)
Mydata1 <- df2genind(locus, ploidy = 1, ind.names = ind, pop = population, sep = "")

tree<-aboot(Mydata1, cutoff=50, quiet=T)

data<-read.table("Data/Isolate_info.txt", na.strings="NA", header=T)
names(data)[names(data) == "Isolate"] <- "taxa"
data$no<-NULL
data<-data %>%
  filter(form=="N", Country=="Norway") 
data <- sapply(data,function(x) {x <- gsub("unknown",NA,x)})
data<-as.data.frame(data)
rownames(data)<-data$taxa

match.phylo.data(tree, data)

mycols=kelly( n=16)
mycols<-mycols[-c(1,2)]
mycols<-mycols[c(14, 1:13)]



tiff("Genetic_distances_FST/dendro_Norw_N_cultivar_circular.tif", res=600, width=3500, height=3500)
color.plot.phylo.custom(tree, data, type="fan", "Cultivar", "taxa",main="Dendrogram of all Norwegian net form isolates \n (colored by barley variety)", cex.main=0.5, col=mycols)
dev.off()



#### Dendro all isolates by form
setwd("//storage-al.slu.se/home$/rawo0001/My Documents/Projects/Netblotch/Manuscript/Analysis2025")
Mydata <- read.table("Data/matrix_all_new_2filtermaf.txt", na.strings="NA", header=T)

#Replace "unknown" with NA
Mydata <- sapply(Mydata,function(x) {x <- gsub("unknown",NA,x)})
Mydata<-as.data.frame(Mydata)
names(Mydata)[names(Mydata) == "Isolate"] <- "taxa"



locus<-Mydata[,-c(1:12)]

ind <- as.character(Mydata$taxa) 
population <- as.character(Mydata$taxa)
Mydata1 <- df2genind(locus, ploidy = 1, ind.names = ind, pop = population, sep = "")

tree<-aboot(Mydata1, cutoff=50, quiet=T)

data<-read.table("Data/Isolate_info.txt", na.strings="NA", header=T)
names(data)[names(data) == "Isolate"] <- "taxa"
data$no<-NULL
data <- sapply(data,function(x) {x <- gsub("unknown",NA,x)})
data<-as.data.frame(data)
rownames(data)<-data$taxa
match.phylo.data(tree, data)



mycolors=c("red", "blue", "black")
tiff("Results/Genetic_distances_FST/dendro_all_form._circular.tiff", res=600, width=4500, height=4500)
color.plot.phylo.custom(tree, data, type="fan", "form", "taxa",main="Dendrogram of all Norwegian net form isolates \n (colored by form)", cex.main=0.5, col=mycolors)
dev.off()
