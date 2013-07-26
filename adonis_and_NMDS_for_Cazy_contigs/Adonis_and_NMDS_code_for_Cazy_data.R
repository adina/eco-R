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

adonis(contig.trans~data$SoilFrac, permutations=9999)


