#### Automatic GEO data microarray sorting and annotating
#### Search and replace 'geo' with what you want the data set to be called and replace the following objects with the necessary information
#### The rest should run automatically
FC=.2
name="MartinezDataGEOfile"
output=paste0("/Users/cansav091/Desktop/Neurogenomics Data/Martinez M0M1M2 Lists/",name)
input="/Users/cansav091/Desktop/Neurogenomics Data/Martinez M0M1M2 Lists/Original Data"

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

gse="GSE5099" ### Put GSE accession number here
data.name="Macrophage Martinez data"

############### Download SOFT file from GEO #########################
setwd(input)
if(file.exists(paste0(gse,".soft"))==FALSE){
  getGEOfile(gse,destdir=input)
  gunzip(paste0(gse,".soft.gz"))
}
geo.soft <- getGEO(filename=paste0(gse,".soft"))
## This part may take some time depending on how many samples and how large the file is. 

############### Get sample info from GEO soft files##################
geo.sample.info=matrix(nrow=length(geo.soft@header$sample_id),ncol=3)
for (ii in 1:length(geo.soft@header$sample_id)){
  geo.sample.info[ii,1:2]=c(geo.soft@header$sample_id[ii],geo.soft@gsms[[ii]]@header$title)
}
warning("The variable 'groups' needs to be tailored to the names of your choosing")
print(geo.sample.info[,2])
geo.groups=c(rep("M0 T0",3),rep("M0 T3",3),rep("M0 T7",3),rep("M1",3),rep("M2",3))

geo.sample.info[,3]=geo.groups

###################### Retrieving and Matching Annotation ##########
### Check that there is annotation for this at Bioconductor#########
gpl=as.vector(geo.soft@header$platform_id)
data.type=c()
for(ii in 1:length(gpl)){
xx=grep("<tr valign=\"top\"><td nowrap>Title</td>" ,readLines(paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=",gpl[ii]),n=300))+1
data.type[ii]=readLines(paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=",gpl[ii]))[xx]
}

geo.data=matrix(nrow=0,ncol=15)
gene.name=c()
################ This loop will repeat for each platform type here ##########
for(jj in 1:length(gpl)){
  ################ Retrieve GEO data matrix file########################
  setwd(input)
  if(length(grep(gpl[jj],dir(),value=TRUE))<1){
    getGEO(gse,destdir=getwd())
    gunzip(grep("txt.gz",dir(),value=TRUE))
  }
  data.file=grep(gpl[jj],dir(),value=TRUE)
  data.file=grep("seriesmatrix.txt",data.file,value=TRUE)
  begins=which(readLines(data.file)=="!series_matrix_table_begin")
  geo.data.tmp=as.matrix(read.table(data.file,skip=begins+1,sep="\t",row.names=1,fill=T,))
  geo.probe.ids=rownames(geo.data.tmp)[-nrow(geo.data.tmp)]
  sample.IDs.tmp=as.vector(t(read.table(data.file,skip=begins,sep="\t",nrows=1)[-1]))
  colnames(geo.data.tmp)=sample.IDs.tmp
  
  ######## Load the annotation for this microarray using Bioconductor #########################
  #array.type=unlist(strsplit(data.type[jj],"\\[|\\]"))[2]## If this is microarray data it will download microarray annotation
  #array.type=tolower(gsub("\\_","",array.type))
  #array.type=gsub("-","",array.type)

  #library(paste0(array.type,".db"),character.only=TRUE)
  #xx=eval(parse(text=paste0(array.type,"SYMBOL")))
  #probeids <- mappedkeys(xx)
  #yy=as.data.frame(xx[probeids])
  #gene.name.tmp=yy$symbol[match(probeids,yy$probe_id)] #This tells you in order of the queried entries,what the coordinates are of the second group of entries only  the first match is given 
  #gene.name.tmp=gene.name.tmp[match(geo.probe.ids,probeids)]
  
  ### Download annotation using the GPL platform ################### 
  if(file.exists(paste0(gpl[jj],".annot"))==FALSE){
    getGEOfile(gpl[jj],destdir=input,AnnotGPL=TRUE)
    gunzip(paste0(gpl[jj],".annot.gz"))
  }
  begins=grep("!platform_table_begin",readLines(paste0(gpl[jj],".annot")))
  gpl.data.tmp=read.table(paste0(gpl[jj],".annot"),sep="\t",skip=begins+1,fill=T,quote="",stringsAsFactors = FALSE)
  assign(paste0(tolower(gpl[jj]),".annot"),gpl.data.tmp,envir=.GlobalEnv)
  xx=match(geo.probe.ids, gpl.data.tmp$V1)
  gene.name.tmp=gpl.data.tmp[xx[!is.na(xx)],3]
  
  #Get rid of probes without gene names
  geo.data.tmp=geo.data.tmp[!is.na(gene.name.tmp),]
  gene.name.tmp=gene.name.tmp[!is.na(gene.name.tmp)]
  
  geo.data=rbind(geo.data,geo.data.tmp)
  gene.name=c(gene.name,gene.name.tmp)
}


####### Get rid of genes/probes that aren't very abundant ##########
####### Get rid of bottom quartile genes ###########################
# (renamed gene.avg to be called probe.avg)
probe.avg=apply(geo.data,1,mean)
cutoff=.10
cutoff=sort(probe.avg)[nrow(geo.data)*cutoff]

geo.groups=c(rep("M0 T0",3),rep("M0 T3",3),rep("M0 T7",3),rep("M1",3),rep("M2",3))
geo.abund=geo.data[which(probe.avg> cutoff),]
gene.abund=gene.name[which(probe.avg> cutoff)]

geo.groups=geo.groups[-c(1:6)]
geo.abund=geo.abund[,-c(1:6)]

pval.probe=c()
fstat.probe=c()
for(ii in 1:length(gene.abund)){
  xx=summary(aov(geo.abund[ii,]~as.factor(geo.groups),data=as.data.frame(geo.abund[ii,])))## Summarizes aov object into the report stats.
  fstat.probe[ii]=xx[[1]][["F value"]][1] ### Stores F statistic in the same vector as genes with multiple probes
  pval.probe[ii]=xx[[1]][["Pr(>F)"]][1]
}
#### If there are multiple probes per gene, use Hotelling's to find DE genes #########
if(length(unique(gene.name))!=length(gene.name)){
  ######### Hotelling's MANOVA #####################################
  fstat=c()#### A bunch of empty vectors to store stats in. 
  hotelling.trace=c()
  pval=c()
  num.probes=c() # number of probes per gene
  df1=c()
  df2=c()
  ctl.sd=c()
  for(ii in 1:length(unique(gene.abund))){ 
    xx=which(gene.abund==unique(gene.abund)[ii])  ### Finds row number index for all the probes for a particular gene
    if(length(xx)<2){
      ctl.sd[ii]=sd(geo.abund[xx,which(geo.groups=="M0 T7")])
    }else{
      ctl.sd[ii]=sd(apply(geo.abund[xx,],2,mean)[which(geo.groups=="M0 T7")])
    }
    num.probes[ii]=p=length(xx) ### Takes note of how many probes are in each set. Sets "p" as the dimensions
    xx=as.matrix(geo.abund[xx,]) ## Takes the data from the row indices indicating previously indicated 
    n=ncol(xx) ## Number of samples which is always 12 for this dataset.
    
    if(p<2){ ## When there is only 1 probe for a gene, This set of code does normal ANOVA and by passes the rest of the code to the next loop
      xx=summary(aov(xx~as.factor(geo.groups),data=as.data.frame(xx)))## Summarizes aov object into the report stats.
      fstat[ii]=xx[[1]][["F value"]][1] ### Stores F statistic in the same vector as genes with multiple probes
      pval[ii]=xx[[1]][["Pr(>F)"]][1] ## Stores p value for 1 probe sets with the other pvalues
      df1[ii]=NA 
      df2[ii]=NA 
      hotelling.trace[ii]=NA
      next
    }else{ ## For genes with more than one probe: 
      xx=manova(t(xx)~as.factor(geo.groups))
      ss.resid =crossprod(xx$residuals) ## Finds the crossproduct of residuals matrix
      ss.effect=crossprod(xx$effects[-1,]) ## Finds the crossproduct of effect matrix, leaves out the intercept, keeps only the effects for each group
      
      if(1/kappa(ss.resid)<1e-10) {## If estimated condition number for the inverse matrix is too small (which is what happens to all the p>>n scenarios; When there are so many probes per gene)
        A1=ss.effect%*%condreg(ss.resid,8)$invS # Those gene sets regularized by the function condreg from Won et al, 2013
        
      }else{ # if the inverse of the matrix can be found, the variance regularization step is skipped
        A1=ss.effect%*%solve(ss.resid) ### 
      }
      
      m=length(xx$assign) ### This comes out to number of groups +1 for intercept
      u=(n-m-p-1)/2 ### Here is the traditional u is calculated. 
      u=ifelse(u<0,0,u) ## However, when p>>n u becomes negative, which throws off the df calculation,
      #so here if it turns out that u is negative, we use 0 instead. (I don't know if this is a good way to do it.)
      t=(abs(p-m+1)-1)/2 ## t gets abnormally large when p>>n, don't know if this is something that needs to be accounted for. 
      s=min(c(p,m-1))### will usually be 3, however, for probesets where there are only two probes, then s=2
      df1[ii]=s*(2*t+s+1)### will be between 6 and 24 for "normal" p<n data. Our data however has a mean df1 of 28 and a max of 321. 
      ###Don't know how these super large df1's should be accounted for. 
      df2[ii]=2*(s*u+1) #### will be between 2 and 14 for our data. 
      ### with a really large p, this "should" become negative, however, that doesn't make sense as a df.
      ## So here I've limited it to be no lower than 2. But, again, don't know how this should really be dealt with. 
      hotelling.trace[ii]=tr(A1) ### Hotelling stat is based on the trace of the matrix found from ss.effects*ss.residuals^-1
      fstat[ii]=(tr(A1)/s)*(df2[ii]/df1[ii])## Adjust F statistic to be proportional to the degrees of freedom for effects and residuals
      pval[ii]=1-pf(fstat[ii], df1[ii], df2[ii]) ## calculates the p value based on the f distribution and df parameters
    }
  }
  length(which(pval<.05)) #Before FDR correction
  pval.fdr=p.adjust(pval,method="hochberg") 
}

#################################################################
######Put basic info in a summary file ##########################
basic.info=as.vector(c(geo.soft@header$geo_accession,geo.soft@header$summary[1],"hgu133A and hgu133b",table(geo.groups),length(unique(geo.groups)),summary(probe.avg)[1:6],nrow(geo.data),sum(!is.na(gene.name)),length(unique(gene.name)),length(which(pval<.05)),length(which(pval.fdr<.05))))
names(basic.info)=c("GSE","Experiment Description:","Data/Microarray Type:",sort(unique(geo.groups)),"Number of Total Samples:","Min:","1st Qu:","Median:","Mean:", "3rd Qu:","Max:","Number of probes:","Number of Probes with Gene Names:","Number of Genes:","Number of p<.05 Genes before FDR","Number of p<.05 Genes FDR")
setwd(outDir)
write.table(basic.info,file=paste(data.name,"Basic Info File M0 M1 M2.txt"),quote=FALSE)

############################################################################
############## Find average for each probe/gene by groups###################
geo.avg=matrix(nrow=nrow(geo.abund),ncol=length(unique(geo.groups)))
for (ii in 1:length(geo.abund[,1])){
  tmp=tapply(geo.abund[ii,],geo.groups,FUN=mean)
  geo.avg[ii,]=tmp
}

############################################################################
########### Calculate fold changes across groups ###########################
dimnames(geo.avg)[[2]]=unique(geo.groups)
dimnames(geo.avg)[[1]]=gene.abund

nn=length(unique(geo.groups))
columns.names=c()
comp=matrix(nrow=2,ncol=0)
for(ii in 1:nn){
  tmp=matrix(ncol=nn,nrow=length(gene.abund))
  for(jj in 1:nn){
    tmp[,jj]=geo.avg[,ii]/geo.avg[,jj]
    columns.names=c(columns.names,paste0(unique(geo.groups)[ii],"\\",unique(geo.groups)[jj]))
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
write.csv(geo.fc,paste(data.name,"Fold Changes M0 M1 M2.csv"),quote=FALSE)

#### Descriptive statistics for the probesets for each gene ###########
if(exists("hotelling.trace")==TRUE){
  sd.probes=c(1:length(unique(gene.abund)))
  avg.group.signal=matrix(nrow=length(unique(gene.abund)),ncol=length(unique(geo.groups)))
  for(ii in 1:length(unique(gene.abund))){
    xx=which(gene.abund==unique(gene.abund)[ii])
    if(length(xx)<2){
      avg.group.signal[ii,]=geo.avg[xx,]
      sd.probes[ii]=NA
      next
    }else{
      avg.group.signal[ii,]=apply(geo.avg[xx,],2,mean)
      sd.probes[ii]=sd(apply(geo.abund[xx,],1,mean))
    }
  }
}
if(exists("hotelling.trace")){
  signif=cbind(ctl.sd,hotelling.trace,fstat,pval,pval.fdr,num.probes,sd.probes,avg.group.signal)
  rownames(signif)=unique(gene.abund)
  colnames(signif)=c("Ctl SD","Hotelling Trace","F-stat","p value","FDR adj p val","Num of Probesets","SD of probsets",unique(geo.groups))
}else{
  signif=cbind(fstat,pval,pval.fdr,avg.group.signal)
  rownames(signif)=unique(gene.abund)
  colnames(signif)=c("F-stat","p value","FDR adj p val",unique(groups))
}

########Put the significance info into a csv file #######################
setwd(output)
write.csv(signif[order(signif[,4]),],file=paste(data.name,"Genes M0 M1 M2.csv"))

#################### Plots #############################################
genes.abund.edit=c()
for(ii in 1:length(strsplit(gene.abund,"///"))){
genes.abund.edit[ii]=unlist(strsplit(gene.abund,"///")[[ii]][1])
}

setwd(output)
jpeg("FC vs ANOVA -log10pval Martinez Data.jpeg",width=1500,height=1500)
xx=cor(abs(geo.fc[,1]),-log10(pval.probe))
plot(abs(geo.fc[,1]),-log10(pval.probe),pch="",xlab="Absolute Value for M1/M0 FC",cex.lab=2 )
textxy(abs(geo.fc[,1]),-log10(pval.probe),genes.abund.edit,offset=0,cex=1.3)
abline(h=1.30,col="red")
abline(v=1.30,col="blue")
dev.off()


ctl.sd.probe=apply(geo.abund[,which(geo.groups=="M0 T7")],1,sd)
jpeg("FC vs Control SD zoomed.jpeg",width=1500,height=1500)
plot(abs(geo.fc[,1]),ctl.sd.probe,pch="",xlab="FC",ylab="Control SD",cex.lab=2,xlim=c(0,5))
textxy(abs(geo.fc[,1]),ctl.sd.probe,genes.abund.edit,offset=0,cex=1.3)
abline(h=1.30,col="red")
abline(v=1.30,col="blue")
dev.off()






