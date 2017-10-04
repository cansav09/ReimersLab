### Comparing Zhang data to lps controls 

zhang=readRDS("mouse.counts.rds")

tmp=match(toupper(dimnames(zhang)[[1]]),toupper(rownames(microglia.genes)))

sum(!is.na(tmp)) ## N=3687/4428 genes that match

matched.fc.lps.z=microglia.genes[match(toupper(dimnames(zhang)[[1]]),toupper(rownames(microglia.genes))),1]
matched.fc.lps.z=matched.fc.lps[!is.na(matched.fc.lps)]

tmp=rownames(microglia.genes)[tmp]
matched.microglia.genes.zhang=tmp[!is.na(tmp)]
tmp=match(matched.microglia.genes.zhang,dimnames(lps.avg)[[1]])
vehicle.microglia.avg=lps.avg[tmp,5]

matched.genes.lps.z=lps.avg[!is.na(matched.fc.lps.z)]

matched.genes.zhang1=zhang[!is.na(tmp),11]
matched.genes.zhang2=zhang[!is.na(tmp),12]

matched.genes.zhangavg=apply(cbind(matched.genes.zhang1,matched.genes.zhang2),1,mean)

### correlate Zhang average with the vehicle microglia average
plot(log2(matched.genes.zhangavg),log2(vehicle.microglia.avg), xlab="log2 Signal Zhang RNA-Seq",ylab="log2 Signal Srinvasan RNA-Seq",ylim=c(-5,10))
abline(0,0)
cor(log2(matched.genes.zhangavg),log2(vehicle.microglia.avg)) ## R=.25


### Correlating each sample individually
plot(log2(matched.genes.zhang1),log2(vehicle.microglia.avg), xlab="log2 Signal Zhang RNA-Seq",ylab="log2 Signal Srinvasan RNA-Seq",pch=".")
abline(log2(matched.genes.zhang1),log2(vehicle.microglia.avg))
cor(log2(matched.genes.zhang1),log2(vehicle.microglia.avg)) ## R=.25

plot(log2(matched.genes.zhang2),log2(vehicle.microglia.avg), xlab="log2 Signal Zhang RNA-Seq",ylab="log2 Signal Srinvasan RNA-Seq")
abline(log2(matched.genes.zhang2),log2(vehicle.microglia.avg))
cor(log2(matched.genes.zhang2),log2(vehicle.microglia.avg)) ## R=.25





