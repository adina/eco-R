#inputting data
data<-read.csv(file.choose())
#I reorganized the data first to have samples as rows and contigs as columns; I brute forced it in excel
#check data
head(data[,1:5])
dim(data)
library(vegan)

#I am reorganizing the data for easier plotting later
data$SoilFrac<-factor(data$SoilFrac, levels=c("Micro","SM","MM","LM","WS"))

#I am transforming the data for relative abundances in each row; this is needed to test community composition

contig.trans<-decostand(data[,4:365],"total")

#perform an NMDS for visualization
contig.nmds<-metaMDS(contig.trans,k=3, autotransform=FALSE)

#plot it in ggplot using this function I made in NMDS_for_community.R (just run the whole thing): put the colors you want in the order that I reordered soil fractions before.  The format for the function is ggplot.NMDS("the nmds you want to plot","groups to color","colors you want")

colors<-c("black","blue","purple","red","grey")
ggplot.NMDS(contig.nmds, data$SoilFrac, colors)

#now to run the adonis

adonis(contig.trans~data$SoilFrac, permutations=9999)### notsignficant 

####Now performing analysis at Cazy_fam2, everything is summed at that cazy level within samples
cazy.2.data<-read.csv("contigs_summed_at_cazy_fam_2.csv")
head(cazy.2.data)
dim(cazy.2.data)
cazy.2.data$SoilFrac<-factor(cazy.2.data$SoilFrac, levels=c("Micro","SM","MM","LM","WS"))
cazy.2.trans<-decostand(cazy.2.data[,3:8],"total")
cazy.2.nmds<-metaMDS(cazy.2.trans,autotransform=FALSE)
colors<-c("black","blue","purple","red","grey")
ggplot.NMDS(cazy.2.nmds, cazy.2.data$SoilFrac, colors)
adonis(cazy.2.trans~cazy.2.data$SoilFrac, permutations=9999)  ###not significant


####Now performing analysis at Cazy_fam, everything is summed at that cazy level within samples
cazy.data<-read.csv("contigs_summed_at_cazy_fam.csv")
head(cazy.data)
dim(cazy.data)
cazy.data$SoilFrac<-factor(cazy.data$SoilFrac, levels=c("Micro","SM","MM","LM","WS"))
cazy.trans<-decostand(cazy.data[,3:73],"total")
cazy.nmds<-metaMDS(cazy.trans,k=3,autotransform=FALSE) #went to more k's, usually do this to keep stress under .1
colors<-c("black","blue","purple","red","grey")
ggplot.NMDS(cazy.nmds, cazy.data$SoilFrac, colors)
adonis(cazy.trans~cazy.data$SoilFrac, permutations=9999)  ###not significant


