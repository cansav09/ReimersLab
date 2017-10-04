### Download Healy et al Macrophage and Microglia data and sort it. Then find DE genes


name="HealyDataGEOfile"
mainDir="/Users/cansav091/Desktop/Neurogenomics Data/Healy M0 M1 M2Lists/"
input="/Users/cansav091/Desktop/Neurogenomics Data/Healy M0 M1 M2Lists/Original Data"
output="/Users/cansav091/Desktop/Neurogenomics Data/Healy M0 M1 M2Lists/HealyDataGEOfile"

############ Directory Commands ##################
# Type this before creating any output files
setwd(output)

### Need these packages: 
library(Biobase)
library(GEOquery)
library(psych)
library(CondReg)
library(Biostrings)
library(calibrate)
library(MASS)

gse="GSE76737" ### Put GSE accession number here
data.name="Healy Data"

############### Download SOFT file from GEO #########################
setwd(input)
if(file.exists(paste0(gse,".soft"))==FALSE){
  getGEOfile(gse,destdir=input)
  gunzip(paste0(gse,".soft.gz"))
}
geo.soft <- getGEO(filename=paste0(gse,".soft"))
## This part may take some time depending on how many samples and how large the file is. 
############### Get sample info from GEO soft files##################
geo.sample.info=matrix(nrow=length(geo.soft@header$sample_id),ncol=4)
for (ii in 1:length(geo.soft@header$sample_id)){
  geo.sample.info[ii,1:2]=c(geo.soft@header$sample_id[ii],geo.soft@gsms[[ii]]@header$title)
}
geo.groups=unlist(strsplit(geo.sample.info[,2],"_"))

geo.sample.info[,3]=geo.groups[seq(from=1,to=length(geo.groups),by=3)]
geo.sample.info[,4]=geo.groups[seq(from=2,to=length(geo.groups),by=3)]


###################### Retrieving and Matching Annotation ##########
### Check that there is annotation for this at Bioconductor#########
gpl=as.vector(geo.soft@header$platform_id)
data.type=c()
ii=1
xx=grep("<tr valign=\"top\"><td nowrap>Title</td>" ,readLines(paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=",gpl[ii]),n=300))+1
data.type[ii]=readLines(paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=",gpl[ii]))[xx]


geo.data=matrix(nrow=0,ncol=15)
gene.name=c()
################ This loop will repeat for each platform type here ##########
  ################ Retrieve GEO data matrix file########################
  setwd(input)

  if(length(grep(gpl[jj],dir(),value=TRUE))<1){
    getGEO(gse,destdir=getwd())
    gunzip(grep("txt.gz",dir(),value=TRUE))
  }
  data.file=paste0(gse,"_series_matrix.txt")
  begins=which(readLines(data.file)=="!series_matrix_table_begin")
  geo.data=as.matrix(read.table(data.file,skip=begins+1,sep="\t",row.names=1,fill=T))
  geo.probe.ids=rownames(geo.data)[-nrow(geo.data)]
  sample.IDs=as.vector(t(read.table(data.file,skip=begins,sep="\t",nrows=1)[-1]))
  colnames(geo.data)=sample.IDs
  
  ######## Load the annotation for this microarray using Bioconductor #########################
  array.type=unlist(strsplit(data.type,"\\[|\\]"))[2]## If this is microarray data it will download microarray annotation
  array.type=tolower(gsub("\\_","",array.type))
  array.type=gsub("-","",array.type)

  library(hugene20sttranscriptcluster.db)
  xx=hugene20sttranscriptclusterSYMBOL
  hugene20sttranscriptcluster
  probeids <- mappedkeys(xx)
  yy=as.data.frame(xx[probeids])
  gene.name=yy[match(probeids,yy$probe_id),] #This tells you in order of the queried entries,what the coordinates are of the second group of entries only  the first match is given 
  xx=match(geo.probe.ids,gene.name$probe_id)
  
  gene.name=gene.name[xx[!is.na(xx)],]
  geo.data=geo.data[!is.na(xx),]



####### Get rid of genes/probes that aren't very abundant ##########
####### Get rid of bottom quartile genes ###########################
# (renamed gene.avg to be called probe.avg)
probe.avg=apply(geo.data,1,mean)
geo.abund=geo.data[which(probe.avg> summary(probe.avg)[2]),]
gene.abund=gene.name[which(probe.avg> summary(probe.avg)[2]),2]
gene.abund=gene.abund[!is.na(gene.abund)]
geo.abund=geo.abund[!is.na(gene.abund),]

groups=as.factor(paste(geo.sample.info[,3],geo.sample.info[,4]))

#xx=c(which(groups=="microglia M0"),which(groups=="microglia M1"))
#geo.abund=geo.abund[!is.na(gene.abund),xx]
#groups=groups[xx]
#### If there are multiple probes per gene, use Hotelling's to find DE genes #########
  ######### Hotelling's MANOVA #####################################
if(ncol(geo.abund)<nrow(geo.abund)){
  geo.abund=t(geo.abund)
}  
  fstat=matrix(ncol=1,nrow=ncol(geo.abund))#### A bunch of empty vectors to store stats in. 
  pval=matrix(ncol=1,nrow=ncol(geo.abund))
  for(ii in 1:length(gene.abund)){ 
      xx=aov(geo.abund[,ii]~groups)## Summarizes aov object into the report stats.
      xx=summary(xx)
      fstat[ii,]=xx[[1]][["F value"]][1] ### Stores F statistic in the same vector as genes with multiple probes
      pval[ii,]=xx[[1]][["Pr(>F)"]][1] ## Stores p value for 1 probe sets with the other pvalues
  }
  length(which(pval[,1]<.05)) #Before FDR correction
  #length(which(pval[,2]<.05)) 
  
  pval.fdr=pval
  pval.fdr[,1]=p.adjust(pval[,1],method="hochberg") 
  #pval.fdr[,2]=p.adjust(pval[,2],method="hochberg") 
  

#################################################################
######Put basic info in a summary file ##########################
#basic.info=as.vector(c(geo.soft@header$geo_accession,geo.soft@header$summary[1],"hgu133A and hgu133b",table(groups),length(unique(groups)),summary(probe.avg)[1:6],nrow(geo.data),sum(!is.na(gene.name)),length(unique(gene.name)),length(which(pval[,1]<.05)),length(which(pval[,2]<.05)),length(which(pval.fdr[,1]<.05)),length(which(pval.fdr[,2]<.05))))
#names(basic.info)=c("GSE","Experiment Description:","Data/Microarray Type:",sort(unique(geo.groups)),"Number of Total Samples:","Min:","1st Qu:","Median:","Mean:", "3rd Qu:","Max:","Number of probes:","Number of Probes with Gene Names:","Number of Genes:","Cell Type:Number of p<.05 Genes before FDR","Treatment:Number of p<.05 Genes before FDR","Cell Type:Number of p<.05 Genes FDR","Treatment:Number of p<.05 Genes FDR")
  
basic.info=as.vector(c(geo.soft@header$geo_accession,geo.soft@header$summary[1],"hgu133A and hgu133b",table(groups),length(unique(groups)),summary(probe.avg)[1:6],nrow(geo.data),sum(!is.na(gene.name)),length(unique(gene.name)),length(which(pval[,1]<.05)),length(which(pval.fdr[,1]<.05))))
names(basic.info)=c("GSE","Experiment Description:","Data/Microarray Type:",sort(unique(groups)),"Number of Total Samples:","Min:","1st Qu:","Median:","Mean:", "3rd Qu:","Max:","Number of probes:","Number of Probes with Gene Names:","Number of Genes:","Number of p<.05 Genes before FDR","Number of p<.05 Genes FDR")
setwd(output)
write.table(basic.info,file=paste(data.name,"Basic Info File M0 vs M1 Microglia.txt"),quote=FALSE)

groups=as.factor(as.character(groups))
############################################################################
############## Find average for each probe/gene by groups###################
geo.avg=matrix(nrow=length(gene.abund),ncol=length(unique(groups)))
for (ii in 1:ncol(geo.abund)){
  tmp=tapply(geo.abund[,ii],groups,FUN=mean)
  geo.avg[ii,]=t(tmp)
}
############################################################################
########### Calculate fold changes across groups ###########################
dimnames(geo.avg)[[2]]=unique(groups)
dimnames(geo.avg)[[1]]=gene.abund

nn=length(unique(groups))
columns.names=c()
comp=matrix(nrow=2,ncol=0)
for(ii in 1:nn){
  tmp=matrix(ncol=nn,nrow=length(gene.abund))
  for(jj in 1:nn){
    tmp[,jj]=geo.avg[,ii]/geo.avg[,jj]
    columns.names=c(columns.names,paste0(unique(groups)[ii],"\\",unique(groups)[jj]))
    comp=cbind(comp,rbind(ii,jj))
  }
  if(ii==1){
    geo.fc=tmp
  }else{
   geo.fc=cbind(geo.fc,tmp)
  }
}
geo.fc=geo.fc[,-which(comp[1,]==comp[2,])]
colnames(geo.fc)=columns.names[-which(comp[1,]==comp[2,])]
comp=comp[,-which(comp[1,]==comp[2,])]
rownames(geo.fc)=gene.abund

########### Put fold changes into a .csv file ##############################
setwd(output)
write.csv(geo.fc,paste(data.name,"Fold Changes M0 vs M1 Microglia.csv"),quote=FALSE)

#### Descriptive statistics for the probesets for each gene ###########
  signif=cbind(fstat,pval,pval.fdr,geo.avg)
  rownames(signif)=gene.abund
  #colnames(signif)=c("F-stat Cell Type","F-stat Treatment","p value Cell Type","p value Treatment","FDR adj p val Cell Type","FDR adj p val Treatment",unique(as.factor(paste(geo.sample.info[,3],geo.sample.info[,4]))))
  colnames(signif)=c("F-stat","p value","FDR adj",unique(groups))
  

########Put the significance info into a csv file #######################
setwd(output)
write.csv(signif[order(signif[,4]),],file=paste(data.name,"Genes M0 vs M1 Microglia.csv"))


