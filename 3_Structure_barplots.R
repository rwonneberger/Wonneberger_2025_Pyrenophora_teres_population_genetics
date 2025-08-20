#devtools::install_github('royfrancis/pophelper')
library(pophelper)
library(data.table)
library(ggplot2)
library(gridExtra)

mergeQ2 <- function(qlist) {
  
  is.qlist(qlist)
  if(diff(range(as.integer(tabulateQ(qlist)$ind)))!=0) stop("mergeQ: Number of individuals differ between runs.")
  
  # Computes mean cell-wise across dataframes
  # x A list of numeric dataframes
  # 
  mergy <- function(x) {
    return(list(Reduce(`+`, x)/length(x)))
  }
  
  # if all runs have same K, merge as is
  if(diff(range(as.integer(tabulateQ(qlist)$k)))==0) {
    labels <- summariseQ(tabulateQ(qlist))$k
    x <- mergy(qlist)
    names(x) <- labels
  }else{
    # if runs have different K, split them and merge within sublists
    qlist <- sortQ2(qlist)
    labels <- summariseQ(tabulateQ(qlist,sorttable=FALSE))$k
    x <- unlist(lapply(splitQ(qlist),mergy),recursive=FALSE)
    names(x) <- labels
  }
  
  return(as.qlist(x))
}

sortQ2 <- function(qlist,by="k",decreasing=FALSE,debug=FALSE) {
  
  is.qlist(qlist)
  if(length(by)==0) stop("sortQ: Argument 'by' must not be length zero.")
  if(!is.character(by)) stop("sortQ: Argument 'by' must be a character.")
  if(!is.logical(decreasing)) stop("sortQ: Argument 'decreasing' must be a logical datatype.")
  
  fun1 <- function(x) as.matrix(unlist(attributes(x)))
  a <- lapply(qlist,fun1)
  if(debug) print(a)
  if(any(!sapply(a,function(x) any(grepl(paste0(by,collapse="|"),rownames(x)))))) {
    stop(paste0("One or more of the attributes provided in by (",by,") is missing in one or more runs. If 'ind' or 'k' is missing, use 'as.qlist()' to add them."))
  }
  
  # get df of attributes
  b <- as.data.frame(t(as.data.frame(lapply(a,function(x,y) x[y,],by),stringAsFactors=FALSE)),stringsAsFactors=FALSE)
  fun2 <- function(x) if(all(!is.na(as.numeric(as.character(x))))) {return(as.numeric(as.character(x)))}else{return(x)}
  b <- as.data.frame(sapply(b,fun2),stringAsFactors=FALSE)
  
  if(debug) {print(str(b)); print(b)}
  
  # order
  ord <- do.call(order,b[,by,drop=FALSE])
  if(decreasing) ord <- rev(ord)
  # sort qlist
  return(qlist[ord])
}


source("https://raw.githubusercontent.com/rwonneberger/R_functions/main/Ronja_functions.R")


meta<-fread("Isolate_names.txt")

sfiles <- list.files(path="Teres/2500050000/Results/",pattern="2500050000_run_*_")
# basic usage
slist <- sortQ2(readQ(sfiles))
if(length(unique(sapply(slist,nrow)))==1) slist <- lapply(slist,"rownames<-",meta$Isolate)
slist_1 <- alignK(slist)
slist_2 <- mergeQ2(slist_1)

slist_3<-slist_2[c(2:14)]
class(slist)
slist_3<-as.qlist(slist_3)




meta <- read.delim("meta_for_pophelper.txt", header=T,stringsAsFactors=F, na.strings = "un")
head(meta)
head(metasub)
metasub <- meta[c(2,11, 4,7, 8, 10)]

#Make multiline plot for K=4 and 5


tiff("Figs/multi_k5.tiff", units="cm", res=300, height=32, width=25)

p <- plotQMultiline(slist_3[4],exportplot=F,returnplot=T,spl=43, lpp=8,grplab=metasub,ordergrp=T, useindlab=T, grplabsize=4, indlabsize = 4)
grid.arrange(p$plot[[1]][[1]])

dev.off()



#make two single plots for k=2 and K=3

slist_3<-slist_2[c(2:3)]
class(slist)
slist_3<-as.qlist(slist_3)



meta <- read.delim("meta_for_pophelper.txt", header=T,stringsAsFactors=F, na.strings = "un")
head(meta)
head(metasub)
metasub <- meta[c(2,11)]

plotQ(slist_3,imgoutput="join",sharedindlab=F,grplab=metasub,
      selgrp="form",ordergrp=T, sortind="Cluster1", showlegend=T, 
      
      legendkeysize = 9, 
      legendtextsize = 6,
      indlabheight=0.07,
      grplabsize=2.1,
      splabsize=6,
      basesize=9, 

      indlabspacer=-1,divgrp = c("form", "Cultivar"), divtype=1,splab = c("k=2", "k=3"),grplabangle=0,  height=2.5, width=17.4,units="cm", dpi=600,exportpath="Figs",
      barbordercolour="black",barbordersize=0,outputfilename="K2_3",imgtype="tiff")

