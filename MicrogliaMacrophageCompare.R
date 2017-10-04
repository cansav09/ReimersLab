## Comparisons of Cell Types and Treatments in Healy et al data

input="/Users/cansav091/Desktop/Neurogenomics Data/Healy M0 M1 M2Lists/HealyDataGEOfile"
output="/Users/cansav091/Desktop/Neurogenomics Data/Healy M0 M1 M2Lists/MicrogliaMacrophageCompare"

dir()
#################### Plots #############################################

setwd(output)
all.cells=read.csv("Healy Data Genes Two Way ANOVA.csv")

jpeg("FC vs ANOVA -log10pval Healy Data Two Way ANOVA.jpeg",width=1500,height=1500)
xx=cor(geo.avg,-log10(pval))
plot(abs(geo.fc[,1]),-log10(pval),pch="",xlab="Microglia/Macrophage FC",cex.lab=2 )
textxy(abs(geo.fc[,1]),-log10(pval),gene.abund,offset=0,cex=1.3)
abline(h=1.30,col="red")
abline(v=1.30,col="blue")
dev.off()


macrophage.activ=read.csv("Healy Data Genes M0 vs M1 Macrophage.csv" )
microglia.activ=read.csv("Healy Data Genes M0 vs M1 Microglia.csv" )

xx=order(macrophage.activ$X)
yy=order(microglia.activ$X)

jpeg("Activation p values Microglia vs Macrophage.jpeg",width=1500,height=1500)
plot(-log10(macrophage.activ$p.value[xx]),-log10(microglia.activ$p.value[yy]),pch="",xlab="-log10(pval) Macrophage",ylab="-log10(pval) Microglia",cex.lab=2 )
textxy(-log10(macrophage.activ$p.value[xx]),-log10(microglia.activ$p.value[yy]),sort(gene.abund),offset=0,cex=1.3)
abline(h=1.30,col="red")
abline(v=1.30,col="red")
dev.off()


macrophage.activ=read.csv("Healy Data Fold Changes M0 vs M1 Macrophage.csv" )
microglia.activ=read.csv("Healy Data Fold Changes M0 vs M1 Microglia.csv" )

jpeg("Activation FC Microglia vs Macrophage.jpeg",width=1500,height=1500)
plot(log(macrophage.activ$macrophages.M1.macrophages.M0),log(microglia.activ$microglia.M1.microglia.M0),pch="",xlab="FC Macrophage M1/M0",ylab="Microglia M1/M0",cex.lab=2 )
textxy(log(macrophage.activ$macrophages.M1.macrophages.M0),log(microglia.activ$microglia.M1.microglia.M0),microglia.activ$X,offset=0,cex=1.3)
abline(h=1.30,col="red")
abline(v=1.30,col="red")
dev.off()

write.csv(cor(geo.avg),file="Correlations across Macrophages and Microglia and Treatments.csv")

jpeg("M0 vs M1 Microglia Avg.jpeg",width=1500,height=1500)
plot(geo.avg[,1],geo.avg[,3],pch="",xlab=dimnames(geo.avg)[[2]][1],ylab=dimnames(geo.avg)[[2]][3],cex.lab=2 )
textxy(geo.avg[,1],geo.avg[,3],gene.abund,offset=0,cex=1.3)
dev.off()

jpeg("M1 Macrophage vs M1 Microglia Avg.jpeg",width=1500,height=1500)
plot(geo.avg[,3],geo.avg[,8],pch="",xlab=dimnames(geo.avg)[[2]][3],ylab=dimnames(geo.avg)[[2]][8],cex.lab=2 )
textxy(geo.avg[,3],geo.avg[,8],gene.abund,offset=0,cex=1.3)
dev.off()

jpeg("M0 Macrophage vs M0 Microglia Avg.jpeg",width=1500,height=1500)
plot(geo.avg[,3],geo.avg[,8],pch="",xlab=dimnames(geo.avg)[[2]][3],ylab=dimnames(geo.avg)[[2]][8],cex.lab=2 )
textxy(geo.avg[,3],geo.avg[,8],gene.abund,offset=0,cex=1.3)
dev.off()

read.csv

jpeg(paste0("Scatterplot Matrix for ",gene.family[jj],"Individual Probes.jpeg"))
scatterplotMatrix(,reg.line=FALSE,var.labels=gene.family.names[tmp])
dev.off()



