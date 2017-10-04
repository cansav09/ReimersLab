#GOAL:Load the previously processed single cell type data from Darmanis et al. and use it to create a list of genes most differentailly regulated in human microglia 
#Step 1: Load Data and libraries needed
#Step 2: Annotate the data with which samples are which cells. 
#Step 3: Sort by the fold changes between microglia and other cell types
#Step 4: Get rid of genes that don't have enough counts to be considered trustworthy
#Step 5: Create lists to use for somparison to Zhang et al mouse microglia gene list
#Step 6: Check for the differences between MiSeq N=138 and HiSeq N=328 
##(This ends up being irrelevant because basically all the samples we care about are done in HiSeq except for four (3 astrocytes and 1 neuron sample))

## Step 1) Load original GEO data files
setwd('~/Desktop/Neurogenomics Data/CellProportionsProject/GSE67835_RAW/')
fnames=dir()
darmanis = matrix(nr=22088,nc=466)
tmp= read.table(fnames[1])
rownames(darmanis)<-tmp[,1]
for ( ii in 1:466) {
  tmp= read.table(fnames[ii])
  if (length(tmp[[1]]) != 22088) {warning('sample',ii,'has',);next}
  darmanis[,ii] <- tmp[,2] 
}

setwd('..')
dimnames(darmanis)[[2]]<-sub('_.*C','.C',sub('.csv.gz','',fnames))
saveRDS( darmanis, file='human.counts.RDS')

library(Biobase)
library(GEOquery)
library(Rsubread)
library(limma)
library(edgeR)

data.genes=dimnames(darmanis)[[1]]####N=22084 genes, last four rows are summary data

#2) Load the information that matches the cell type to the sample ID and match the cell type to the dataset 
## For whatever reason the GEO files came into two separate files
CellType=as.vector(t(read.table("SampleType.txt",sep="\t")))
CellID=as.vector(t(read.table("SampleIDs.txt",sep="\t")))
CellType=c(CellType,as.vector(t(read.table("SampleType2.txt",sep="\t"))))
CellID=c(CellID,as.vector(t(read.table("SampleIDs2.txt",sep="\t"))))

##combining the cell type to their sample IDs
Type2ID=rbind(substr(CellType,12,nchar(CellType)),CellID)
Type2ID=rbind(Type2ID, c(rep("Hiseq",328),rep("MiSeq",138)))
##matching the CellType ID's to the samples in the human.counts dataset
tmp=match(substr(dimnames(darmanis)[[2]],1,10),CellID)##These are the coordinates of the Cell Type that actually correspond to our dataset 

###Putting the cell type information in the same order as the data. 
Type2ID=Type2ID[,tmp]

####checking that things are in the right order
match(substr(dimnames(darmanis)[[2]],1,10),Type2ID[2,])==1:466

###replace sample IDs with cell types for the column names 
dimnames(darmanis)[[2]]=Type2ID[1,]
###figuring out how many of each cell type we have
summary(as.factor(dimnames(darmanis)[[2]]))

#astrocytes       endothelial   fetal_quiescent fetal_replicating     hybrid 
#62                20               110                25                46 
#microglia        neurons       oligodendrocytes        OPC 
#16               131                38                18 

#3) Convert each log10(CPM) from the original dataset into RPKMs
## Formula:  RPKM=10^(log10CPMs)/genelength/1000
## First need to match the gene names to their transcript lengths

##Match the annotation "NM_000 etc" IDs with the ids that correspond to the probes (referred to as "nm")
gene2NM=readRDS("genes.RDS")
### How many genes make this up? 
length(unique(gene2NM$Gene.Symbol))##N=17319 genes corresponding to #22011 transcripts

match(data.genes,gene2NM$Gene.Symbol)
## match these genes to the dataset's genes 
NMid=c()
NMgeneid=c()
for(i in 1:22088){
nn=grep(paste0("^",data.genes[i],"$"),gene2NM$Gene.Symbol)
NMid=c(NMid,gene2NM$RefSeq[nn])
NMgeneid=c(NMgeneid,gene2NM$Gene.Symbol[nn])
}
# There are N=16392 matches
NMdataids=unique(cbind(NMid,NMgeneid))#N=16200 are unique gene to NM ids

## Load mRNAs for genes from NCBI genome wide FASTA file
library("Biostrings")
FASTA=readDNAStringSet("GCF_000001405.32_GRCh38.p6_rna.fna")
transcript.lengths.key=FASTA@ranges@width

annot.nm=FASTA@ranges@NAMES
# Clean up NM_IDs, the other list doesn't have the decimal points. 
annot.nm=strsplit(annot.nm," ")
annot.nm=sapply(annot.nm, "[", 1)
###Get rid of decimal points since the IDs for the samples do not have them. 
annot.nm=substr(annot.nm,0,nchar(annot.nm)-2)
data.transcript.lengths=transcript.lengths.key[match(NMdataids[,1],annot.nm)]

#Now to find the median transcript length for each gene's transcript set
data.transcript.medians=tapply(as.numeric(data.transcript.lengths),NMdataids[,2],FUN=median,na.rm=TRUE)

##Calculate the RPKMs using the median length of the transcripts for any gene's set of transcript
##Reads per kilobase per million reads.
darmanis.match=darmanis[!is.na(match(data.genes,NMdataids[,2])),]
data.transcript.medians=data.transcript.medians[match(dimnames(darmanis.match)[[1]],NMdataids[,2])]
data.transcript.medians=data.transcript.medians/1000

# calculating the reads per kilobase of transcript
human.rpkm=matrix(ncol=466,nrow = 16097)
for(i in 1:16097){
  human.rpkm[i,]=human.counts.match[i,]/data.transcript.medians[i]
}
dimnames(human.rpkm)=dimnames(darmanis.match)

#calculate total number of mapped reads per sample 
total.reads=apply(darmanis[c(1:22085,22088),],2,sum)### Not including "no feature" but the rest is counted in the sum
total.reads=total.reads/1e6

for(i in 1:466){
  human.rpkm[,i]=human.rpkm[,i]/total.reads[i]
}

human.rpkm=human.rpkm[!is.na(human.rpkm[,1]),]# Left with 15500 genes that aren't NA

##Let's see if these values make any sense: 

plot(apply(human.rpkm,1,mean), ylab="Average Calculated RPKMs for Each Gene")
gene.avg=apply(human.rpkm,1,mean)
hist(gene.avg[gene.avg>1])
min(gene.avg) #0
median(gene.avg)# 2.11
max(gene.avg)# 2432
mean(gene.avg) #11.83

#5) Get the averages by cell type 
human.rpkm.avg=matrix(nrow=15500,ncol=9)
for (i in 1:15500){
  nn=tapply(human.rpkm[i,],Type2ID[1,],FUN=mean)
  human.rpkm.avg[i,]=nn
}

##carry over the cell names and the gene names 
colnames(human.rpkm.avg)=unique(sort(Type2ID[1,]))
rownames(human.rpkm.avg)=rownames(human.rpkm)

human.rpkm.abun=human.rpkm.avg[which(human.rpkm.avg[,6]>1),]

#N=7796 genes with any counts at all 
##N=3218 genes in which microglia have at least 1 rpkm
#N=2775 genes in which microglia have at least 5 rpkm 
rownames(human.rpkm.abun)=rownames(human.rpkm.avg)[which(human.rpkm.avg[,6]>1)] 

#6) Sort by microglia fold change 
##calculate fold change by the average of the averages across the cell types
microglia.fc=human.rpkm.abun[,5]/apply(human.rpkm.abun[,c(1,7,8)],1,mean)

### Fold change only using the oligodensdrocytess,astrocytes, and neurons as comparisons
##Sort genes from largest to smallest fc
microglia.signal=human.rpkm.abun[,5][order(microglia.fc,decreasing=TRUE)]
microglia.genes=rownames(human.rpkm.abun)[order(microglia.fc,decreasing=TRUE)]
nonmicroglia.signal=apply(human.rpkm.abun[,c(1,7,8)],1,mean)[order(microglia.fc,decreasing=TRUE)]
## adjust for cell proportion 
microglia.fc=microglia.fc[order(microglia.fc,decreasing=TRUE)]

### Save this info to text file

write.table(cbind(microglia.fc,microglia.signal,nonmicroglia.signal),file="Human microglia distinctive list 4-28-16CS", quote=FALSE)

##7) Check the differences/similarirites between MiSeq and HiSeq 
summary(as.factor(Type2ID[1,Type2ID[3,]=="MiSeq"]))
### MiSeq sample summary: 
#astrocytes   fetal_quiescent fetal_replicating           neurons 
#3               109                25                 1 

summary(as.factor(Type2ID[1,Type2ID[3,]=="Hiseq"]))

#astrocytes   endothelial     fetal_quiescent     hybrid        microglia           neurons 
#59               20                1               46               16              130 
#oligodendrocytes OPC 
#38               18 

###Since without the fetal samples, this is only affecting 4 samples, that's pretty insignificant, so I'm not going to bother looking into this much further
###But if we use the fetal cells later on, we will have to look into this more. 

##Split different methods into different matrices
human.rpkm.avgMi=human.rpkm.avg[,Type2ID[3,]=="MiSeq"]
human.rpkm.avgHi=human.rpkm.avg[,Type2ID[3,]=="Hiseq"]

### Checking out the consistency among samples 
human.rpkm.microglia=human.rpkm[,which(dimnames(human.rpkm)[[2]]=="microglia")]
microglia.corr=cor(human.rpkm.microglia)
for(i in 1:16){
  microglia.corr[i,i]=NA
}
hist(apply(microglia.corr,2,mean,na.rm=TRUE))
human.rpkm.microglia.concord=human.rpkm.microglia[,which(apply(microglia.corr,2,mean,na.rm=TRUE)>median(apply(microglia.corr,2,mean,na.rm=TRUE)))]

### Now using this data for comparison with Zhang
microglia.fc=human.rpkm.microglia.concord[,5]/apply(human.rpkm.microglia.concord[,c(1,7,8)],1,mean)
### Fold change only using the oligodensdrocytess,astrocytes, and neurons as comparisons
##Sort genes from largest to smallest fc
microglia.signal=human.rpkm.microglia.concord[,5][order(microglia.fc,decreasing=TRUE)]
microglia.genes=rownames(human.rpkm.microglia.concord)[order(microglia.fc,decreasing=TRUE)]
nonmicroglia.signal=apply(human.rpkm.microglia.concord[,c(1,7,8)],1,mean)[order(microglia.fc,decreasing=TRUE)]
## adjust for cell proportion 
microglia.fc=microglia.fc[order(microglia.fc,decreasing=TRUE)]


### Doing the same with neurons since they should be the most easily consistent0
human.rpkm.neuron=human.rpkm[,which(dimnames(human.rpkm)[[2]]=="neurons")]

neuron.corr=cor(human.rpkm.neuron)
for(i in 1:16){
  neuron.corr[i,i]=NA
}
hist(apply(neuron.corr,2,mean,na.rm=TRUE))
median(apply(neuron.corr,2,mean,na.rm=TRUE))

human.rpkm.neuron.concord=human.rpkm.neuron[,which(apply(neuron.corr,2,mean,na.rm=TRUE)>median(apply(neuron.corr,2,mean,na.rm=TRUE)))]



