

M1.eg=M1.genes[20:31,c(8,11)]
#dimnames(M1.eg)[[2]]=c("LPS Treatment","C1q Treatment","Social Stress")
dimnames(M1.eg)[[2]]=c("Microglia Distinctive","SCZ profile")
dimnames(M1.eg)[[1]]=c("Gene A","Gene B","Gene C","Gene D","Gene E","Gene F","Gene G","Gene H","Gene I","Gene J","Gene K","Gene L")
library(gplots)
heatmap <- heatmap(M1.eg, Rowv=TRUE, Colv=NA, col =redgreen(50),cexCol=1.2, scale="column", margins=c(9,10))
dev.off()


xx=rnorm(1000)
yy=xx*-2+rnorm(100)

dev.off()
plot(xx,yy,xlab="Microglial Activation for 'i'",ylab="Antipsychotic Response Score")


library(ggplot2)
TreatedMicroglia <- data.frame(length = rnorm(1000, 3, 4))
Controls <- data.frame(length = rnorm(1000, 15, 5))

#Now, combine your two dataframes into one.  First make a new column in each.
TreatedMicroglia$group <- 'Treated Microglia'
Controls$group<- 'Controls'

#and combine into your new data frame vegLengths
Data <- rbind(TreatedMicroglia,Controls)

#now make your lovely plot
ggplot(Data, aes(length, fill = group)) + geom_density(alpha = 0.2)



