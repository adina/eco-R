##data set is merged metadata, abundance of each contig, cazy info for each contig, and taxonomy of each contig
data<-read.csv(file.choose())

head(data)
library(ggplot2)
library(plyr)

#set up for errorbars
limits<-aes(ymax=MEAN+SE, ymin=MEAN-SE)
names(data)

#summary stats for figure by some level of taxonomy and enzyme
stats<-ddply(data, .(t3,Cazy_fam2), summarise,
MEAN=mean(Abundance, na.rm=T),
SE=sd(Abundance, na.rm=T)/sqrt(length(Abundance))
)

head(stats)
#weird blank rows were removed
stats<-stats[-c(1:3),]

ggplot(stats, aes(t3, MEAN, colour=Cazy_fam2))+geom_pointrange(limits)+coord_flip()+theme_bw()+theme(aspect.ratio=1)+scale_y_log10()

