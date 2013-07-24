data<-read.csv(file.choose())
names(data)

library(vegan)

#transform abundances for relative comparisons on columns that do not contain metadata. In this case I am using a cazy dataset with the first three columns containing sample name, crop, and soil fraction

enz.trans<-cbind(data[,1:3],decostand(data[,4:73],"total"))


#the function is a y~x1+x2 format like this: community matrix~Trt1*Trt2
adonis(enz.trans[,4:73]~enz.trans$Crop*enz.trans$SoilFrac)

for nmds you can do it this way

enz.mds<-metaMDS(enz.trans[,4:73])
plot(enz.mds)