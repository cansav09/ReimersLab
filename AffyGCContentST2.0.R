## Examine GC content of probes and signal values 
gse="GSE76737"

library(installr)
library(seqinr)
library(affxparser)
library(affy)
library(affyPLM)
library(oligo)
library(pd.hugene.2.0.st)
library(oligoData)
library(BSgenome.Hsapiens.UCSC.hg19)

input="/Users/cansav091/Desktop/Neurogenomics Data/Healy M0 M1 M2Lists/Original Data"
output="/Users/cansav091/Desktop/Neurogenomics Data/Healy M0 M1 M2Lists/AffyGCContentST2.0"
raw.data.file=paste0("/Users/cansav091/Desktop/Neurogenomics Data/Healy M0 M1 M2Lists/Original Data/",gse,"_RAW")

setwd(input)

################ Retrieving Probe Level Data ############

#### Load the .CEL data #############
raw.data=getGEOSuppFiles(gse)
setwd(raw.data.file)

chp.files=grep(".chp.gz",dir(),value=TRUE)
cel.files=grep(".CEL.gz",dir(),value=TRUE)

#### Unzip all the CEL files
for(ii in 1:length(cel.files)){
gunzip(cel.files[ii])
}
chp.files=grep(".chp",dir(),value=TRUE)
cel.files=grep(".CEL",dir(),value=TRUE)

##### Get the signal intensities for the probes
rawData <- read.celfiles(cel.files)
signals=oligo::intensity(rawData)

setwd(output)
jpeg("Signal Distribution Boxplots non-normalized hg ST 2.0.jpeg")
boxplot(signals)
dev.off()

pmSeq <- pmSequence(rawData)
ph=rawData@phenoData
pmSeqlog2=log2(pm(rawData))

st20.seq=c()
st20.gc=c()
for(ii in 704472:length(pmSeq)){
  st20.seq=c(st20.seq,as.character(pmSeq[[ii]]))
  st20.gc=c(st20.gc,GC(unlist(strsplit(st20.seq[[ii]],""))))
  }
setwd(output)
jpeg("GC Distribution hg ST 2.0.jpeg")
hist(st20.gc,las=2)
dev.off()


### Create the jpegs for the cel images 
for (ii in 1:length(cel.files))
{
  jpeg(paste0("image",cel.files[ii],".jpg"))
  image(rawData[,ii],main=rownames(ph@data)[ii])
  dev.off()
}

##### Download the probe sequence data
feat = rawData@featureData

##### Normalize the Data 
rma.data=rma(rawData)

setwd(output)
jpeg("Signal Distribution Boxplots RMA hg ST 2.0.jpeg")
boxplot(rma.data,las=2)
dev.off()

avg.sig=apply(pmSeqlog2,1,mean)
for(ii in 1:30){
  jpeg(paste0("log2 PM probes vs GC content hgu133a Sample",ii,".jpeg"))
  plot(st20.gc,pmSeqlog2[,ii]/avg.sig,xlab="GC ratio",ylab="Log2 Signal/Avg Log2 Signal")
  dev.off()
}



