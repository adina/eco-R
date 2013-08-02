setwd("~/Documents/Projects/KBase/Data_for_venn/")

library(ggplot2)
library(reshape)
library(Vennerable)
library(plyr)

data<-read.delim("merged_fy.txt", header=T, row.names=1)
head(data)
length(data[1,])
head(data[, 1:9])
head(data[, 10:27])

#derep by summing
LM <- rowSums(data[10:13])
MI <- rowSums(data[14:17])
MM <- rowSums(data[18:20])
SM <- rowSums(data[21:23])
WS<-rowSums(data[24:27])

new<-cbind(m[, 1:9], LM, MI, MM, SM, WS)

##reclass all factor strings as characters
##i <- sapply(m, is.factor)
##m[i] <- lapply(m[i], as.character)

contig<-cbind(new[,1], new[10:14])
head(contig)

##consolidate duplicated rows by summing
contig.con<-ddply(contig, "new[,1]", numcolwise(sum))

con<-cbind(contig.con[, 2:6])
head(con)
rownames(con)<-contig.con[, 1]
head(con)
m<-as.matrix(con)
head(m)
for (i in 1:nrow(m)){
 for(j in 1:ncol(m)) {
 m[i, j]<-ifelse(m[i,j]>0, rownames(m)[i], m[i, j])
 }
 }

LM<-as.vector(m[, 1][m[, 1]!=0])
MI<-as.vector(m[, 2][m[, 2]!=0])
MM<-as.vector(m[, 3][m[, 3]!=0])
SM<-as.vector(m[, 4][m[, 4]!=0])
WS<-as.vector(m[, 5][m[, 5]!=0])

lv_contig<-vector(mode="list")
lv_contig$LM<-LM
lv_contig$MM<-MM
lv_contig$MI<-MI
lv_contig$SM<-SM
lv_contig$WS<-WS

V_contig<-Venn(lv_contig)
plot(V_contig, doWeights=FALSE, type="ellipses")


#############################################
##Cazy_fam
##consolidate duplicated rows by summing
cazy1<-cbind(new[,3], new[10:14])
head(cazy1)
length(cazy1[,1])
##consolidate duplicated rows by summing
cazy1.con<-ddply(cazy1, "new[,3]", numcolwise(sum))
head(cazy1.con)
length(cazy1.con[, 1])

con<-cbind(cazy1.con[, 2:6])
head(con)
rownames(con)<-cazy1.con[, 1]
head(con)
m<-as.matrix(con)
head(m)
for (i in 1:nrow(m)){
 for(j in 1:ncol(m)) {
 m[i, j]<-ifelse(m[i,j]>0, rownames(m)[i], m[i, j])
 }
 }

LM<-as.vector(m[, 1][m[, 1]!=0])
MI<-as.vector(m[, 2][m[, 2]!=0])
MM<-as.vector(m[, 3][m[, 3]!=0])
SM<-as.vector(m[, 4][m[, 4]!=0])
WS<-as.vector(m[, 5][m[, 5]!=0])

lv_contig<-vector(mode="list")
lv_contig$LM<-LM
lv_contig$MM<-MM
lv_contig$MI<-MI
lv_contig$SM<-SM
lv_contig$WS<-WS

V_contig<-Venn(lv_contig)
plot(V_contig, doWeights=FALSE, type="ellipses")