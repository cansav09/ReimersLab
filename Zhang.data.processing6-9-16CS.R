## Goal: Look at how the other cell types correlate
# Step 1) Load the Zhang original data
# Step 2) Calculate the fold changes for the neurons from Zhang
# Step 3) Calculate the fold changes for then neurons in the Darmanis dataset



# Step 1) Load the raw data from GSE file
setwd('~/Desktop/Neurogenomics Data/CellProportionsProject/GSE52564_RAWcsv/')
fnames=dir()
zhang.data = matrix(nr=22462,nc=18)
tmp= read.csv(fnames[1])
rownames(zhang.data)<-tmp[,1]
for ( i in 1:18) {
  tmp= read.csv(fnames[i])
  zhang.data[,i] <- tmp[,2] 
}
setwd('..')
colnames(zhang.data)<-sub('_.*C','.C',sub('.csv','',fnames))
saveRDS(zhang.data, file='mouse.counts.RDS')

## DO all the gene names match?? yup. 
zhang.genes=data.frame(1:22462)
for ( i in 1:18) {
  tmp= read.csv(fnames[i])
  zhang.genes[,i] <- tmp[,1] 
}
setwd('..')


# Step 2) Calculate the averages for each cell type

CellType=substr(dimnames(zhang.data)[[2]],12,nchar(dimnames(zhang.data)[[2]])-1)
CellType[13]="Endothelial"

zhang.fpkm.avg=matrix(nrow=22462,ncol=7)
for (i in 1:22462){
  nn=tapply(zhang.data[i,],CellType,FUN=mean)
  zhang.fpkm.avg[i,]=nn
}
colnames(zhang.fpkm.avg)=sort(unique(CellType))
dimnames(zhang.fpkm.avg)[[1]]=dimnames(zhang.data)[[1]]

library(calibrate)
zhang.gene.avg=apply(zhang.fpkm.avg,1,mean)
plot(zhang.gene.avg)
textxy(1:22462,zhang.gene.avg,dimnames(zhang.data)[[1]])
dev.off()

### CD45 gene in both samples (marker for microglia)
zhang.fpkm.avg[which(dimnames(zhang.fpkm.avg)[[1]]=="Ptprc"),]
human.rpkm.avg[which(dimnames(human.rpkm.avg)[[1]]=="PTPRC"),]

zhang.fpkm.abun=zhang.fpkm.avg[which(zhang.fpkm.avg[,5]>1),]## N=64029 genes

#6) Get rid of genes that don't have a human match.
library(org.Hs.eg.db)
EGIDs <- keys(org.Hs.eg.db)
genesymbol.key=unlist(mget(EGIDs, org.Hs.egSYMBOL))# N=56340 genes in this homo sapiens database

zhang.fpkm.abun=zhang.fpkm.abun[which(!is.na(match(toupper(dimnames(zhang.fpkm.abun)[[1]]),genesymbol.key))),]
## N=55482 genes that have a human match

#7) Get fold changes for microglia
##calculate fold change by the average of the averages across the cell types
microglia.fc.zhang=zhang.fpkm.abun[,4]/apply(zhang.fpkm.abun[,c(1,6,5)],1,mean)


### Fold change only using the oligodensdrocytess,astrocytes, and microglias as comparisons
##Sort genes from largest to smallest fc
microglia.signal.zhang=zhang.fpkm.abun[,4][order(microglia.fc.zhang,decreasing=TRUE)]
microglia.genes.zhang=rownames(zhang.fpkm.abun)[order(microglia.fc.zhang,decreasing=TRUE)]
nonmicroglia.signal.zhang=apply(zhang.fpkm.abun[,c(1,5,6)],1,mean)[order(microglia.fc.zhang,decreasing=TRUE)]
microglia.fc.zhang=microglia.fc.zhang[order(microglia.fc.zhang,decreasing=TRUE)]

###Make a microglia distinctive gene list file
write.table(cbind(microglia.signal.zhang[1:500],nonmicroglia.signal[1:500],microglia.fc.zhang[1:500]),"MicrogliaDistinct Based on Zhang6-6-16CS",col.names=c("microglia signal","nonmicroglia signal","fold change"))

#6) Get fold changes for neurons
neuron.fc.zhang=zhang.fpkm.abun[,6]/apply(zhang.fpkm.abun[,c(1,4,5)],1,mean)
### Fold change only using the oligodensdrocytess,astrocytes, and neurons as comparisons
##Sort genes from largest to smallest fc
neuron.signal.zhang=zhang.fpkm.abun[,6][order(neuron.fc.zhang,decreasing=TRUE)]
neuron.genes.zhang=rownames(zhang.fpkm.abun)[order(neuron.fc.zhang,decreasing=TRUE)]
nonneuron.signal.zhang=apply(zhang.fpkm.abun[,c(1,4,5)],1,mean)[order(neuron.fc.zhang,decreasing=TRUE)]
## adjust for cell proportion 
neuron.fc.zhang=neuron.fc.zhang[order(neuron.fc.zhang,decreasing=TRUE)]

## Get fold changes for zhang data
astrocytes.fc.zhang=zhang.fpkm.abun[,1]/apply(zhang.fpkm.abun[,c(6,4,5)],1,mean)
### Fold change only using the oligodensdrocytess,astrocytes, and astrocytess as comparisons
##Sort genes from largest to smallest fc
astrocytes.signal.zhang=zhang.fpkm.abun[,1][order(astrocytes.fc.zhang,decreasing=TRUE)]
astrocytes.genes.zhang=rownames(zhang.fpkm.abun)[order(astrocytes.fc.zhang,decreasing=TRUE)]
nonastrocytes.signal.zhang=apply(zhang.fpkm.abun[,c(6,4,5)],1,mean)[order(astrocytes.fc.zhang,decreasing=TRUE)]
## adjust for cell proportion 
astrocytes.fc.zhang=astrocytes.fc.zhang[order(astrocytes.fc.zhang,decreasing=TRUE)]
