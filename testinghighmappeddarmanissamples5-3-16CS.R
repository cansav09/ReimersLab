#### Compare the higher proportion of mapped read samples to zhang and see how that looks 

### For looking at high mapped 
neuron.fc=high.mapped.rpkm.avg[,7]/apply(high.mapped.rpkm.avg[,c(1,6,8)],1,mean)
### Fold change only using the oligodensdrocytess,astrocytes, and neurons as comparisons
##Sort genes from largest to smallest fc
neuron.signal=high.mapped.rpkm.avg[,7][order(neuron.fc,decreasing=TRUE)]
neuron.genes=dimnames(high.mapped.rpkm.avg)[[1]][order(neuron.fc,decreasing=TRUE)]
nonneuron.signal=apply(high.mapped.rpkm.avg[,c(1,6,8)],1,mean)[order(neuron.fc,decreasing=TRUE)]
## adjust for cell proportion 
neuron.fc=neuron.fc[order(neuron.fc,decreasing=TRUE)]
neuron.fc[which(neuron.fc=="Inf")]=8000
####################################################################################
### For looking at low mapped

neuron.fc=low.mapped.rpkm.avg[,6]/apply(low.mapped.rpkm.avg[,c(1,5,7)],1,mean)
### Fold change only using the oligodensdrocytess,astrocytes, and neurons as comparisons
##Sort genes from largest to smallest fc
neuron.signal=low.mapped.rpkm.avg[,6][order(neuron.fc,decreasing=TRUE)]
neuron.genes=dimnames(low.mapped.rpkm.avg)[[1]][order(neuron.fc,decreasing=TRUE)]
nonneuron.signal=apply(low.mapped.rpkm.avg[,c(1,5,7)],1,mean)[order(neuron.fc,decreasing=TRUE)]
## adjust for cell proportion 
neuron.fc=neuron.fc[order(neuron.fc,decreasing=TRUE)]
neuron.fc[which(neuron.fc=="Inf")]=8000
####################################################################################


darmanis=cbind(neuron.genes,neuron.signal,nonneuron.signal,neuron.fc)
zhang=cbind(toupper(neuron.genes.zhang),neuron.signal.zhang,nonneuron.signal.zhang,neuron.fc.zhang)

darmanis=darmanis[which(darmanis[,2]>100),]
zhang=zhang[which(zhang[,2]>100),]


### Let's get rid of genes in the zhang list that don't appear to have any human equivalent
library(org.Hs.eg.db)
EGIDs <- keys(org.Hs.eg.db)
genesymbol.key=unlist(mget(EGIDs, org.Hs.egSYMBOL))# N=56340 genes in this homo sapiens database

###go back and get the full list from Zhang 
#Get rid of zhang genes without human matches 604 genes not matching (by using symbols anyway)
zhang=zhang[which(!is.na(match(zhang[,1],genesymbol.key))),] 

## Now lets make a top 500 genes and cross reference them to zhang
overlap=zhang[match(darmanis[,1],zhang[,1]),1]
overlap=overlap[!is.na(overlap)] # N=359 of the genes that are in the zhang top 500 and are in the darmanis with an rpkm>0

nn=zhang[match(overlap,zhang[,1]),1]
mm=darmanis[match(overlap,darmanis[,1]),1]
which((nn==mm)=="FALSE")### Check gene IDs match


## Plotting the fold change across these two datasets
nn=log2(as.numeric(zhang[match(overlap,zhang[,1]),4]))
mm=log2(as.numeric(darmanis[match(overlap,darmanis[,1]),4]))
cor(mm,nn)

jpeg("HighProportionMappedZhangvsDarmanisFoldChangeNeuronsLOG25-4-16CS.jpeg")
plot(mm,nn,xlab="Zhang Fold Change",ylab="Darmanis Fold Change",sub="R=0.61")
dev.off()

##Plotting abundances across the two datasets 
nn=log2(as.numeric(zhang[match(overlap,zhang[,1]),2]))
mm=log2(as.numeric(darmanis[match(overlap,darmanis[,1]),2]))
cor(mm,nn)

jpeg("HighProportionMappedZhangvsDarmanisSignalNeuronsLOG25-4-16CS.jpeg")
plot(mm,nn,xlab="Zhang Signal",ylab="Darmanis Signal",sub="R=0.22")
dev.off()