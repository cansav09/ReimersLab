## Examine GC content of probes and signal values 
gse="GSE5099"

xx=readLines("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL96")

library(installr)
library(seqinr)
library(affxparser)
library(affy)
library(affyPLM)
library(oligo)
library(pd.hugene.2.0.st)
library(oligoData)
library(pd.hg.u133a)
library(pd.hg.u133b)

################ Retrieving Probe Level Data ############

#### Load the .CEL data #############
raw.data=getGEOSuppFiles(gse)
setwd(raw.data.file)
untar(grep("_RAW.tar",dir(),value=TRUE))

chp.files=grep(".chp.gz",dir(),value=TRUE)
cel.files=grep(".CEL.gz",dir(),value=TRUE)

#### Unzip all the CEL files
for(ii in 1:length(cel.files)){
gunzip(cel.files[ii])
}
A.cel.files=grep("A.CEL",dir(),value=TRUE)
B.cel.files=grep("B.CEL",dir(),value=TRUE)

##### Get the signal intensities for the probes
rawDataA <- read.celfiles(A.cel.files)
rawDataB <- read.celfiles(B.cel.files)

signalsA=oligo::intensity(rawDataA)
signalsB=oligo::intensity(rawDataB)

jpeg("Signal Distribution Boxplots non-normalized hgu133a.jpeg")
boxplot(signalsA)
dev.off()

jpeg("Signal Distribution Boxplots non-normalized hgu133b.jpeg")
boxplot(signalsB)
dev.off()

#### Get the probe sequences
pmSeqA <- pmSequence(rawDataA)
phA=rawDataA@phenoData
pmSeqAlog2=log2(pm(rawDataA))

pmSeqB <- pmSequence(rawDataB)
phB=rawDataB@phenoData
pmSeqBlog2=log2(pm(rawDataB))


data(affyExpressionFS)
pmSequence(affyExpressionFS)

hgu133A.seq=c()
hgu133A.gc=c()
for(ii in 1:length(pmSeqA)){
hgu133A.seq=c(hgu133A.seq,as.character(pmSeqA[[ii]]))
hgu133A.gc=c(hgu133A.gc,GC(unlist(strsplit(hgu133A.seq[[ii]],""))))
}

hgu133B.seq=c()
hgu133B.gc=c()
for(ii in 1:length(pmSeqB)){
  hgu133B.seq=c(hgu133B.seq,as.character(pmSeqB[[ii]]))
  hgu133B.gc=c(hgu133B.gc,GC(unlist(strsplit(hgu133B.seq[[ii]],""))))
}

jpeg("Hgu133A probe GC% distribution.jpeg")
hist(hgu133A.gc,xlab="GC ratio")
dev.off()

jpeg("Hgu133B probe GC% distribution.jpeg")
hist(hgu133B.gc,xlab="GC ratio")
dev.off()

### Create the jpegs for the cel images 
for (ii in 1:length(A.cel.files))
{
  jpeg(paste0("image",A.cel.files[ii],".jpg"))
  image(rawDataA[,ii],main=rownames(phA@data)[ii])
  dev.off()
}
for (ii in 1:length(B.cel.files))
{
  jpeg(paste0("image",B.cel.files[ii],".jpg"))
  image(rawDataB[,ii],main=rownames(phB@data)[ii])
  dev.off()
}

##### Download the probe sequence data
probe.dataA=fitProbeLevelModel(rawDataA)
probe.dataB=fitProbeLevelModel(rawDataB)

feat = rawDataA@featureData
featureNames(rawDataA)

##### Normalize the Data 
rma.dataA=rma(rawDataA)
rma.dataB=rma(rawDataB)

##### Plot normalized data 

setwd("/Users/cansav091/Desktop/Neurogenomics Data/Martinez M0M1M2 Lists/AffyGCContentU133")

jpeg("Signal Distribution Boxplots RMA hgu133a.jpeg")
boxplot(rma.dataA,las=2)
dev.off()

jpeg("Signal Distribution Boxplots RMA hgu133b.jpeg")
boxplot(rma.dataB,las=2)
dev.off()

avg.sigA=apply(pmSeqAlog2,1,mean)
avg.sigB=apply(pmSeqBlog2,1,mean)
for(ii in 1:15){
jpeg(paste0("log2 PM probes vs GC content hgu133a Sample",ii,".jpeg"))
plot(hgu133A.gc,pmSeqAlog2[,ii]/avg.sigA,xlab="GC ratio",ylab="Log2 Signal/Avg Log2 Signal")
dev.off()

jpeg(paste0("log2 PM probes vs GC content hgu133b Sample",ii,".jpeg"))
plot(hgu133B.gc,pmSeqBlog2[,ii]/avg.sigB,xlab="GC ratio",ylab="Log2 Signal/Avg Log2 Signal")
dev.off()
}

## Number of probes and how many probesets have that number of probes
#8    10    11    13    14    15    16 
#1     1 21721     4     4     2   482 

probesets=unlist(strsplit(names(probe.dataA@probe.coefs),"at"))
probesets=paste0(probesets[seq(from=1,to=length(probesets),by=2)],"at")











