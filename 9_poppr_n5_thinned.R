
############################## subpop n5  
source("https://raw.githubusercontent.com/rwonneberger/R_functions/main/Ronja_functions.R")


snp<-read.table("Data/matrix_all_new_2filtermaf_thinned.txt", sep="\t", header=T)
snp<-snp[-c(2:12)]
info<-read.table("Results/Popgen/Poppr_Norw_all_sub_n5.txt", sep="\t", header=T)
mat<-merge(info, snp, by="Isolate")
mat$no<-NULL
mat$number<-NULL
dim(mat)
head(mat)
str(mat)
#write.table(mat, "Popgen/Poppr_matrix_Norw_all_sub_n5.txt", sep="\t", row.names=F)

#Mydata <- read.table("Popgen/Poppr_matrix_Norw_all_sub_n5.txt", na.strings="NA", header=T)
Mydata<-mat
dim(Mydata)
head(Mydata)[1:10, 1:14]
locus<-Mydata[,-c(1:13)]


ind <- as.character(Mydata$Isolate) 
population <- as.character(Mydata$Sub)
Mydata1 <- df2genind(locus, ploidy = 1, ind.names = ind, pop = population, sep = "")
Mydata1<-missingno(Mydata1, type="ignore")

a1<-poppr(Mydata1, sample=999)
  write.table(a1, "Results/Popgen/Poppr_Norw_all_sub_n5_thinned.csv", sep=",", row.names = F)


