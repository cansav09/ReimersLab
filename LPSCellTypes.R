## GOAL: Load GEO data files from RNA-Seq study for LPS treated microglia, neurons, and astrocytes
name="LPSCellTypes"
input="/Users/cansav091/Desktop/Neurogenomics Data/CellTypesAndLPSSrnivasan/Original Data"
outDir=paste0(mainDir,"/",name)
outDirDay=paste0(outDir,"/",as.character(Sys.Date()))
if(dir.exists(outDir)==FALSE){
  dir.create(outDir)
}

################ Load original GEO data file with sample info#############################
library(GEOquery)
library(Biobase)
library(psych)
library(CondReg)
library(Biostrings)
library(calibrate)
library(MASS)
#########Upload the soft file that I downloaded from GEO with the sample info##############
gse="GSE75246" ### Put GSE accession number here
data.name="Srnivasan LPS Data"

############### Download SOFT file from GEO #########################
setwd(input)
if(file.exists(paste0(gse,".soft"))==FALSE){
  getGEOfile(gse,destdir=input)
  gunzip(grep(paste0(gse,".soft.gz"),dir(),value=TRUE))
}
geo.soft <- getGEO(filename=grep("GSE",grep(".soft$",dir(),value=TRUE),value=TRUE))
## This part may take some time depending on how many samples and how large the file is. 

## Sample id's from the 
sample.id=gse.file@header$sample_id

######################Fetch the cell type of each sample#####################
tmp=gse.file@gsms$GSM1947162@header$title
tmp=strsplit(tmp," - ")
cell.type=c()
for (ii in 1:28){
 cell.type[ii]=tmp[[ii]][2]
}
##Fetch the treatment type of each sample
treat.type=c()
for (ii in 1:28){
  treat.type[ii]=tmp[[ii]][3]
}
geo.sample.info=cbind(geo.soft@header$sample_id,cell.type,treat.type)

########################### Look for annotation for the data #################
gpl=as.vector(geo.soft@header$platform_id)
data.type=c()
for(ii in 1:length(gpl)){
  xx=grep("<tr valign=\"top\"><td nowrap>Title</td>" ,readLines(paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=",gpl[ii]),n=300))+1
  data.type[ii]=readLines(paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=",gpl[ii]))[xx]
}

############# Download GEO files and load into a matrix ############
geo.data=matrix(ncol=nrow(geo.sample.info)/length(gpl),nrow=length(geo.soft@gsms[[1]]@dataTable@table$VALUE))
gene.name=c()

setwd(input)
if(length(grep(paste0(gse,"series_matrix.txt.gz"),dir(input)))<1){
  getGEO(gse,destdir=input)
  gunzip(grep("txt.gz",dir(),value=TRUE))
}
fnames=dir()

lps.data = matrix(nr=length(read.table(fnames[1],skip=5)[,1]),nc=28)
tmp= read.table(fnames[1],skip=5)
dimnames(lps.data)[[1]]<-tmp[,1]

######Load in the 2nd column of data from each file (third column is counts if that is needed)
for (ii in 1:28) {
  tmp= read.table(fnames[ii],skip=5)
  if (length(tmp[,1]) != N) {
    warning('sample ',ii,' does not have 28059 lines')
    next
    }
  lps.data[,ii] <- tmp[,2] ## reading in the RPKM's only
}

############### Reset the working directory to the normal file#########################
setwd('..')
dimnames(lps.data)[[2]]<-sub('_.*C','.C',sub('.txt','',fnames))## carry over the sample IDs
saveRDS(lps.data, file='lps.treated.cell.types.RDS')

###############Get the Gene names for these things#########################
library(org.Mm.eg.db)

EGIDs <- keys(org.Mm.eg.db)
genesymbol.key=unlist(mget(EGIDs, org.Mm.egSYMBOL))

tmp=genesymbol.key[match(dimnames(lps.data)[[1]],EGIDs)]
gene.names=unlist(tmp)

dimnames(lps.data)[[1]]=gene.names
write.csv(genesymbol.key,file="GeneName-ProbeGPL17021.csv",quote=FALSE)
##### There were 28059 gene entries but only 26,815 have gene names in EGID ########
lps.data=lps.data[!is.na(dimnames(lps.data)[[1]]),]
###### Removed those without matching gene names

###############Let's see if these values make any sense: ######################################################

plot(apply(lps.data,1,mean), ylab="Average Calculated RPKMs for Each Gene")
gene.avg=apply(lps.data,1,mean)
plot(gene.avg,ylim=c(0,250))
summary(gene.avg)

## Getting the average reads for each sample
sample.avg=apply(lps.data,2,mean)
par(mar=c(8,3,3,3))
barplot(sample.avg,las=2,ylim=c(0,4.0))

## Looking at the average reads for each sample
barplot(sample.avg[order(group)],names.arg=sort(group),las=2)

barplot(tapply(sample.avg,group,mean),las=2)

##Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
##0.0000    0.0149    0.4123    1.7160    1.7320 2510.0000 

#Step 4) Sort out only the abundant genes and non-outliers
lps.abund=lps.data[which(gene.avg>1),]
gene.names.abund=gene.names[which(gene.avg>1)]

### Looking at the highest expressed genes
gene.avg[order(gene.avg,decreasing=TRUE)][1:100]


## Make this into a plot for all the samples' skewness
library(ggplot2)
jpeg("Skewness.jpeg")
qplot(sk.data,as.factor(group),ylab="Sample Name",main="Sample Skewness",xlab="Skewness")
dev.off()

### Checking the skewness for each samples' data 
library(e1071)
skewness(lps.abund[,1])


sk.data=c()
for (ii in 1:28){
  jpeg(paste("LPS Data RNA-Seq Distribution ",group[ii],ii,".jpeg"))
  hist(log2(lps.data[,ii]))
  sk.data[ii]=skewness(lps.abund[,ii])
  dev.off()
}

######################################################################
## Check out how much these sample do or do not look like each other. 
cor.abund.lps=cor(lps.abund[1:28,])
install.packages("GGally")
library(GGally)

for (ii in 1:6){
  xx=cor(lps.abund[which(lps.avg[,ii]<200),which(group==unique(group)[1])])
  jpeg(paste("Correlation matrix low abundance genes",unique(group)[ii],".jpeg"))
  pairs(lps.abund[which(lps.avg[,ii]<200),which(group==unique(group)[ii])],main=paste("Correlations of",unique(group)[ii],"R=",substr(mean(xx),1,5)),pch=20,xlim=c(0,10),ylim=c(0,10))
  dev.off()
}


##################################################################
### Making a hierarchical cluster of the sample types

jpeg("Cluster Analysis of Samples w IDs.jpeg",width=600,height=500)
plot(hclust(dist(t(lps.abund)[1:28,1:10087])),xlab="Samples")
dev.off()
library(dplyr)
library(ggplot2)
library(plyr)
library(dendextend)

color.group=mapvalues(group,from = c(unique(group)), to = c("blue","red","brown","green","black","purple"))

jpeg("Cluster Analysis of Samples w IDs.jpeg",width=700,height=500)
par(mar=c(10,3,3,3))
dend <- t(lps.abund)[1:28,1:10087] %>%  scale %>% dist %>% hclust %>% as.dendrogram
dend %>% plot
labels(dend)
dend %>% set("labels_col", color.group[order(unlist(dend))]) %>% set("labels_cex", 1)  %>% plot(main = "Cluster of Cell Types/Treatments")

legend(175,c(unique(group)),lty=c(1,1), lwd=c(2.5,2.5),col=c(color.group),cex=.75) 
dev.off()


##################### Get the averages by cell type ################################
group=paste(treat.type, cell.type)
lps.avg=matrix(nrow=length(dimnames(lps.abund)[[1]]),ncol=6)
for (ii in 1:length(dimnames(lps.abund)[[1]])){
  tmp=tapply(lps.abund[ii,],group,FUN=mean)
  lps.avg[ii,]=tmp
}
rownames(lps.avg)=gene.names.abund
colnames(lps.avg)=c(unique(group))

### Summaries by cell/treatment type##########################################
summaries.cell=matrix(ncol=16,nrow=6)
for(ii in 1:6) {
  summaries.cell[ii,]=c(summary(lps.avg[,ii]),paste(rownames(lps.avg)[order(lps.avg[,ii],decreasing=TRUE)[1:10]]))
}
colnames(summaries.cell)=c("Min.", "1st Qu.","Median","Mean", "3rd Qu.","Max.",rep("Top Genes",10)) 
rownames(summaries.cell)=sort(unique(group))

cell.type.reads= apply(lps.avg,2,mean)
barplot(cell.type.reads,las=2,ylim=c(0,4.0))
par(mar=c(8,3,3,3))
  
##carry over the cell names and the gene names 
colnames(lps.avg)=unique(sort(group))
rownames(lps.avg)=rownames(lps.abund)

##Astrocyte gene avg correlations with lps and control
cor(lps.avg[,1],lps.avg[,4])## R=.978
plot(lps.avg[,1],lps.avg[,4],xlim=c(0,300),ylim=c(0,300),xlab="Astrocyte lps avg RPKM",ylab="Astrocyte control avg RPKM")

##Microglia gene avg correlations with lps and control
cor(lps.avg[,2],lps.avg[,5])## R=.989
plot(lps.avg[,2],lps.avg[,5],xlim=c(0,300),ylim=c(0,300),xlab="Microglia lps avg RPKM",ylab="Microglia control avg RPKM")

##Neuron gene avg correlations with lps and control
cor(lps.avg[,3],lps.avg[,6])## R=.997
plot(lps.avg[,3],lps.avg[,6],xlim=c(0,300),ylim=c(0,300),xlab="Neuron lps avg RPKM",ylab="Neuron control avg RPKM")


############## Using ANOVA to explore possible interaction effects ############################
lps.abund=t(lps.abund)
lps.abund=as.data.frame(lps.abund)

lps.data.anova=matrix(nrow=ncol(lps.abund),ncol=6)
for (ii in 1:ncol(lps.abund)){
 lps.data.anova[ii,1:4]= summary(aov(lps.abund[,ii]~cell.type*treat.type,data=lps.abund))[[1]][["Pr(>F)"]]

 }

#######  Benjamini Hochberg for multiple correction ################################################
lps.data.anova[,4]=p.adjust(lps.data.anova[,1],method = "hochberg")
lps.data.anova[,5]=p.adjust(lps.data.anova[,2],method = "hochberg")
lps.data.anova[,6]=p.adjust(lps.data.anova[,3],method = "hochberg")

signif.cell=lps.abund[,which(lps.data.anova[,4]<.05)]
# 4739 genes are significant after p value correction 

### post hoc on genes that are significant 
cell.tukey.pval=matrix(ncol=3,nrow=ncol(signif))
treat.tukey.pval=matrix(ncol=1,nrow=ncol(signif))
interaction.tukey.pval=matrix(ncol=15,nrow=ncol(signif))
cell.tukey.diff=matrix(ncol=3,nrow=ncol(signif))
treat.tukey.diff=matrix(ncol=1,nrow=ncol(signif))
interaction.tukey.diff=matrix(ncol=15,nrow=ncol(signif))

for (ii in 1:ncol(signif)){
  xx=aov(signif[,ii]~cell.type*treat.type,data=lps.abund)
  tmp=TukeyHSD(x=xx, c('cell.type','treat.type','cell.type:treat.type'), conf.level=0.95)
  cell.tukey.pval[ii,]=t(tmp$cell.type[,4])
  treat.tukey.pval[ii,]=t(tmp$treat.type[,4])
  interaction.tukey.pval[ii,]=t(tmp$'cell.type:treat.type'[,4])
  cell.tukey.diff[ii,]=t(tmp$cell.type[,1])
  treat.tukey.diff[ii,]=t(tmp$treat.type[,1])
  interaction.tukey.diff[ii,]=t(tmp$'cell.type:treat.type'[,1])
}
colnames(cell.tukey.pval)=dimnames(tmp$cell.type)[[1]]
colnames(treat.tukey.pval)=dimnames(tmp$treat.type)[[1]]
colnames(interaction.tukey.pval)=dimnames(tmp$'cell.type:treat.type')[[1]]
rownames(cell.tukey.pval)=dimnames(signif)[[2]]
rownames(treat.tukey.pval)=dimnames(signif)[[2]]
rownames(interaction.tukey.pval)=dimnames(signif)[[2]]

colnames(cell.tukey.diff)=dimnames(tmp$cell.type)[[1]]
colnames(treat.tukey.diff)=dimnames(tmp$treat.type)[[1]]
colnames(interaction.tukey.diff)=dimnames(tmp$'cell.type:treat.type')[[1]]
rownames(cell.tukey.diff)=dimnames(signif)[[2]]
rownames(treat.tukey.diff)=dimnames(signif)[[2]]
rownames(interaction.tukey.diff)=dimnames(signif)[[2]]



##Gathering data for genes that are changed by lps in microglia############
## Want to use only genes for which microglia have a big enough expression of but according to each cell type###
tmp=lps.avg[which(apply(cbind(lps.avg[,2],lps.avg[,5]),1,mean)>1),]

## Comparing lps microglia to all other cells. 
microglia.fc=tmp[,2]/(apply(tmp[,c(1,3,4,5,6)],1,sum)/5)
microglia.fc=microglia.fc[which(cell.tukey.pval[,1]<.05|cell.tukey.pval[,3]<.05)]### 4739 genes are greater than 1.5 FC
microglia.genes=cbind(cell.tukey.pval[,2:3],microglia.fc,tmp[match(attr(microglia.fc,"names"), rownames(tmp)),])
dimnames(microglia.genes)[[2]]=c("microglia-astrocyte","neuron-microglia","fold.change",sort(unique(group)))
write.csv(microglia.genes,file="microglia.lps.vs.othercells.csv",quote=FALSE)

## Comparing lps microglia to control microglia

microglia.lps.fc=tmp[,2]/tmp[,5]
microglia.lps.fc=microglia.lps.fc[which(interaction.tukey.pval[,8]<.05)]### 6740 genes are greater than 1.5 FC
microglia.lps.genes=cbind(interaction.tukey.pval[which(interaction.tukey.pval[,8]<.05),8],microglia.lps.fc,tmp[match(attr(microglia.lps.fc,"names"), rownames(tmp)),])
dimnames(microglia.lps.genes)[[2]]=c("lps vs vehicle","fold.change",sort(unique(group)))
write.csv(microglia.lps.genes,file="microglia.lps.vs.microglia.control.csv",quote=FALSE)

hist(microglia.lps.fc)

##############################################################################################################
### Do the same thing for Neurons!

tmp=lps.avg[which(apply(cbind(lps.avg[,3],lps.avg[,6]),1,mean)>1),]## 9996 are expressed in neuron

## Comparing lps neuron to all other cells. 
neuron.fc=tmp[,3]/(apply(tmp[,c(1,2,4,5,6)],1,sum)/5)
neuron.fc=neuron.fc[which(cell.tukey.pval[,2]<.05|cell.tukey.pval[,3]<.05)]### 4428 genes are greater than 1.5 FC
neuron.genes=cbind(cell.tukey.pval[which(cell.tukey.pval[,2]<.05|cell.tukey.pval[,3]<.05),2:3],neuron.fc,tmp[match(attr(neuron.fc,"names"), rownames(tmp)),])
dimnames(neuron.genes)[[2]]=c("neuron-astrocyte","neuron-microglia","fold.change",sort(unique(group)))
write.csv(neuron.genes,file="neuron.lps.vs.othercells.csv",quote=FALSE)

## Comparing lps neuron to control neuron
neuron.lps.fc=tmp[,3]/tmp[,6]
neuron.lps.fc=neuron.lps.fc[which(interaction.tukey.pval[,12]<.05)]### 7025 genes are greater than 1.5 FC
neuron.lps.genes=cbind(interaction.tukey.pval[which(interaction.tukey.pval[,12]<.05),12],neuron.lps.fc,tmp[match(attr(neuron.lps.fc,"names"), rownames(tmp)),])
dimnames(neuron.lps.genes)[[2]]=c("lps vs vehicle","fold.change",sort(unique(group)))
write.csv(neuron.lps.genes,file="neuron.lps.vs.neuron.control.csv",quote=FALSE)

########################################################################################
### Do the same thing for astrocytes!
tmp=lps.avg[which(apply(cbind(lps.avg[,1],lps.avg[,6]),1,mean)>1),]## 9810 are expressed in astrocyte

## Comparing lps astrocyte to all other cells. 
astrocyte.fc=tmp[,1]/(apply(tmp[,c(1,2,4,5,6)],1,sum)/5)
astrocyte.fc=astrocyte.fc[which(cell.tukey.pval[,1]<.05|cell.tukey.pval[,1]<.05)]### 4702 genes are greater than 1.5 FC
astrocyte.genes=cbind(cell.tukey.pval[which(cell.tukey.pval[,1]<.05|cell.tukey.pval[,1]<.05),1:2],astrocyte.fc,tmp[match(attr(astrocyte.fc,"names"), rownames(tmp)),])
dimnames(astrocyte.genes)[[2]]=c("microglia-astrocyte","astrocyte-neuron","fold.change",sort(unique(group)))
write.csv(astrocyte.genes,file="astrocyte.lps.vs.othercells.csv",quote=FALSE)

## Comparing lps astrocyte to control astrocyte
astrocyte.lps.fc=tmp[,1]/tmp[,6]
astrocyte.lps.fc=astrocyte.lps.fc[which(interaction.tukey.pval[,3]<.05)]### 7064 genes are greater than 1.5 FC
astrocyte.lps.genes=cbind(interaction.tukey.pval[which(interaction.tukey.pval[,3]<.05),3],astrocyte.lps.fc,tmp[match(attr(astrocyte.lps.fc,"names"), rownames(tmp)),])
dimnames(astrocyte.lps.genes)[[2]]=c("lps vs vehicle","fold.change",sort(unique(group)))
write.csv(astrocyte.lps.genes,file="astrocyte.lps.vs.astrocyte.control.csv",quote=FALSE)


##Compare the lps differential list to the one derived from the martinez et al dataset
mart.M1=read.csv("M1 Specific Genes Data.csv")

tmp=match(mart.M1$Gene.Symbol,toupper(rownames(microglia.lps.genes)))

sum(!is.na(tmp))

matched.fc.lps=microglia.lps.genes[match(mart.M1$Gene.Symbol,toupper(rownames(microglia.lps.genes))),1]
matched.fc.lps=matched.fc.lps[!is.na(matched.fc.lps)]

matched.genes.mart=mart.M1[!is.na(tmp),2]
matched.fc.mart=mart.M1[!is.na(tmp),6]

plot(log2(matched.fc.mart),log2(matched.fc.lps), xlab="log2 Fold Change Martinez Microarray",ylab="log2 Fold Change Srinvasan RNA-Seq")
abline(log2(matched.fc.mart),log2(matched.fc.lps))
cor(log2(matched.fc.mart),log2(matched.fc.lps)) ## R=.25

## Note that Martinez is IFN gamma + lps and also is macrophages
## 66 genes of the 154 genes in the Martinez list show up as differentially changed in Martinez et al
