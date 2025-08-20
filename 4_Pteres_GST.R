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



################### Helium vs Tiril

#All isolates
#Read in data.frame and save only rows with marker data in locus
Mydata <- read.table("Data/matrix_all_new_2filtermaf.txt", na.strings="NA", header=T)

Mydata<-Mydata %>%
  filter(Cultivar=="Helium"|Cultivar=="Tiril", form=="N", Country=="Norway")


locus<-Mydata[,-c(1:12)]


#Replace "unknown" with NA
Mydata <- sapply(Mydata,function(x) {x <- gsub("unknown",NA,x)})
Mydata<-as.data.frame(Mydata)


#Turn locus into Genind object
ind <- as.character(Mydata$Isolate) 
population <- as.character(Mydata$Cultivar)
Mydata1 <- df2genind(locus, ploidy = 1, ind.names = ind, pop = population, sep = "")
Mydata1
nAll(Mydata1)
div <- summary(Mydata1)
div


#Calculate Hedrick's GST between subgroups
gst<-Gst_Nei(Mydata1)
#### Global GST between Norw N Hel and Tir 0.3084996


#Turn output into df and format
gst<-as.data.frame(gst)
gst$rs<-rownames(gst)
rownames(gst)<-NULL

#Add chrom and pos for visualization
Hap<-read.table("Data/Hapmap_2maffilter.txt", header=T)
gst1<-merge(gst, Hap, by="rs")
gst2<-gst1[,c(1, 2, 5, 6)]
gst2$chrom[gst2$chrom == 'CM016795_1'] <- '1'
gst2$chrom[gst2$chrom == 'CM016796_1'] <- '2'
gst2$chrom[gst2$chrom == 'CM016797_1'] <- '3'
gst2$chrom[gst2$chrom == 'CM016798_1'] <- '4'
gst2$chrom[gst2$chrom == 'CM016799_1'] <- '5'
gst2$chrom[gst2$chrom == 'CM016800_1'] <- '6'
gst2$chrom[gst2$chrom == 'CM016801_1'] <- '7'
gst2$chrom[gst2$chrom == 'CM016802_1'] <- '8'
gst2$chrom[gst2$chrom == 'CM016803_1'] <- '9'
gst2$chrom[gst2$chrom == 'CM016804_1'] <- '10'
gst2$chrom[gst2$chrom == 'CM016805_1'] <- '11'
gst2$chrom[gst2$chrom == 'CM016806_1'] <- '12'

gst3<-gst2%>%filter(!chrom%in%c("NPOS01000013_1", "NPOS01000019_1", "NPOS01000020_1", "NPOS01000029_1"))
gst3$chrom<-as.numeric(gst3$chrom)

gst3 = gst3 %>% dplyr::select(rs, chrom, pos, per.locus)
colnames(gst3)[4]<-"HelTir"


gst3%<>%filter(HelTir !="NaN")
gst3$GST_roll<-zoo::rollmean(gst3$HelTir,10,fill = list(NA, NULL, NA))


scaleFUN <- function(x) sprintf("%.1f", x)




tiff("Figs/GST_HelTir.tiff", units="cm", res=600, height=6, width=17.4)
ggplot(gst3,aes(x=as.numeric(pos)/1000,y=HelTir))  + geom_point(size=0.7)+   theme_classic() + theme(panel.grid.minor = element_line(colour = "grey90"), panel.grid.major = element_line(colour = "grey90"), panel.border = element_rect(fill=NA))+ 
  scale_x_continuous( labels = function(x) format(x, scientific = FALSE))+
  theme(axis.text.x = element_text(angle = 90, vjust=.5)) + 
  facet_grid(~chrom, scales='free_x')+ xlab("Position (kbp)")+ ylab("GST")+ scale_y_continuous(labels=scaleFUN)+
  scale_y_continuous(labels=scaleFUN) + theme(strip.text.x = element_text(size =5.5), axis.text=element_text(size=6), axis.title=element_text(size=8))

dev.off()
