setwd("~/Documents/Projects/KBase/Data_for_venn/")

library(ggplot2)
library(reshape)
library(Vennerable)
library(plyr)

abundance_data<-read.delim("merged.abundance.ann.txt", header=T, row.names=1, strip.white=TRUE)
head(abundance_data)
length(abundance_data[1,])
head(abundance_data[, 2:19])

#filter samples. get rid of contigs of which each sample contains less than min_count. 
min_count = 5
in_all_samples <- subset(abundance_data, abundance_data$PF_LM_H08 >=min_count |abundance_data$PF_LM_H14 >=min_count | abundance_data$PF_LM_H16 >=min_count | abundance_data$PF_LM_H03 >=min_count | abundance_data$PF_MI_H01 >=min_count |abundance_data$PF_MI_H06 >=min_count | abundance_data$PF_MI_H12 >=min_count |abundance_data$PF_MI_H13 >=min_count | abundance_data$PF_MM_H17 >=min_count |abundance_data$PF_MM_H19 >=min_count | abundance_data$PF_MM_H20 >=min_count |abundance_data$PF_SM_H02 >=min_count | abundance_data$PF_SM_H10 >=min_count |abundance_data$PF_SM_H11 >=min_count | abundance_data$PF_WS_H04 >=min_count |abundance_data$PF_WS_H07 >=min_count | abundance_data$PF_WS_H09 >=min_count |abundance_data$PF_WS_H15 >=min_count)
length(in_all_samples[,1])
length(abundance_data[,1])

head(in_all_samples[, 2:19])


#derep by summing
LM <- rowSums(in_all_samples[2:5])
MI <- rowSums(in_all_samples[6:9])
MM <- rowSums(in_all_samples[10:12])
SM <- rowSums(in_all_samples[13:15])
WS<-rowSums(in_all_samples[16:19])

new<-cbind(in_all_samples[, 1], in_all_samples[, 20:27], LM, MI, MM, SM, WS)
head(new)
length(new[1,])

##reclass all factor strings as characters
##i <- sapply(m, is.factor)
##m[i] <- lapply(m[i], as.character)

contig<-cbind(new[,1], new[10:14])
head(contig)

##consolidate duplicated rows by summing
contig.con<-ddply(contig, "new[,1]", numcolwise(sum))
length(contig[, 1])
length(contig.con[, 1])

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

#############################################
##Cazy_fam2
##consolidate duplicated rows by summing
cazy2<-cbind(new[,4], new[10:14])
head(cazy2)
length(cazy2[,1])
##consolidate duplicated rows by summing
cazy2.con<-ddply(cazy2, "new[,4]", numcolwise(sum))
head(cazy2.con)
length(cazy2.con[, 1])

con<-cbind(cazy2.con[, 2:6])
head(con)
rownames(con)<-cazy2.con[, 1]
head(con)
m<-as.matrix(con)
head(m)

#############################################
##t3
##consolidate duplicated rows by summing
t3<-cbind(new[,9], new[10:14])
head(t3)
length(t3[,1])
##consolidate duplicated rows by summing
t3.con<-ddply(t3, "new[,9]", numcolwise(sum))
head(t3.con)
length(t3.con[, 1])

con<-cbind(t3.con[, 2:6])
head(con)
rownames(con)<-t3.con[, 1]
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

##################
## To find out which cazy family groups or taxa gorups that are unique to each fraction, do following grouping and sort. starting from matrix m. 
head(m)
unique.lm<-cbind(m[, 3], m[, 2], m[, 4])
m.new<-cbind(m, rowSums(unique.lm))
head(m.new)
sort.unique.lm<-m.new[order(m.new[, 6]), ]
head(sort.unique.lm)
