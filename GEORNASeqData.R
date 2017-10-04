#### Automatic GEO data RNA-Seq sorting and annotating
#### Search and replace 'geo' with what you want the data set to be called and replace the following objects with the necessary information
#### The rest 'should' run automatically

### Need these packages: 
library(Biobase)
library(GEOquery)
library(psych)
library(CondReg)
library(Biostrings)

gse="GSE"#####" ### Put GSE accession number here
data.name="#####" ### Put what you want the data to be called here.

############### Download SOFT file from GEO #########################
getGEOfile(gse,destdir=getwd())
gunzip(paste0(gse,".soft.gz"))
geo.soft <- getGEO(filename=paste0(gse,".soft"))
## This part may take some time depending on how many samples and how large the file is. 

###################### Print the summary so you can alter the code if needed ###########
geo.soft@header$summary[1]
geo.soft@header$overall_design[1]

groups=c() ## Name your groups with somehting that will match the sample names

###################### Retrieving and Matching Annotation ##########
### Check that there is annotation for this at Bioconductor#########
gpl=geo.soft@header$platform_id
data.type=grep("<tr valign=\"top\"><td nowrap>Title</td>" ,readLines(paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=",gpl),n=300))+1
data.type=readLines(paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=",gpl))[data.type]

############# Get sample info from SOFT files #########################
geo.sample.info=matrix(nrow=#########,ncol=3)
for (ii in 1:length(geo.data[1,])){
  geo.sample.info[ii,1:3]=c(geo.soft@header$sample_id[ii],geo.soft@gsms[[ii]]@header$title[ii],geo.soft@gsms[[ii]]@header$organism_ch1[[ii]])
}
geo.groups=as.vector(nrow(geo.sample.info))
for (ii in 1:length(groups)){
  geo.groups[grep(groups[ii],geo.sample.info[,2])]=groups[ii]
  
  gene.name=rownames(geo.data)
}


if(length(grep("seq",data.type))<0){
#############For RNA-Seq Data####################################
  data.type=unlist(strsplit(data.type,"justify|td"))[3]## If this is microarray data it will download microarray annotation
  data.type=substr(data.type,3,nchar(data.type)-2)
  
  suppfiles.list=getGEOSuppFiles(gse)
  setwd(paste0(getwd(),"/",gse))
  fnames=dir()
  if(length(fnames)==1){
    geo.data=read.csv(fnames[1],row.names = 1)
  }else{
    tmp= read.csv(fnames[1])
    geo.data = matrix(nr=nrow(tmp),nc=length((geo.soft@gsms)))
    rownames(geo.data)<-tmp[,1]
    for ( ii in 1:length((geo.soft@gsms))) {
      tmp= read.csv(fnames[ii])
      geo.data[,ii] <- tmp[,2] 
    }
  }
  setwd('..')
  #colnames(geo.data)<-sub('_.*C','.C',sub('.csv','',fnames))
  saveRDS(geo.data, file=paste0(data.name,'RDS'))
  
  ##### Checks if data is counts and converts to fpkms ###################
  if(integer(geo.data[,1])){
    setwd("/Users/cansav091/Desktop/Neurogenomics Data/Transcriptomes")
    
    ############### Download mouse annotation ###########################
    if(!is.na(grep("Mus musculus",geo.sample.info[,3]))){
      RNA=readDNAStringSet("Mouse.mrna.fa")
      library(org.Mm.eg.db)
      geneids=as.data.frame(org.Mm.egSYMBOL[mappedkeys(org.Mm.egSYMBOL)])
      rnaids=as.list(org.Mm.egACCNUM[mappedkeys(org.Mm.egACCNUM)])
      xx=match(geneids$gene_id,mappedkeys(org.Mm.egACCNUM))
      geneids=geneids[!is.na(xx),2]
      rnaids=rnaids[xx[!is.na(xx)]]
    }
    ############### Download human annotation ###########################
    if(!is.na(grep("Homo sapiens",geo.sample.info[,3]))){
      RNA=readDNAStringSet("Human.mrna.fa")
      library(org.Hs.eg.db)
      geneids=as.data.frame(org.Hs.egSYMBOL[mappedkeys(org.Hs.egSYMBOL)])
      rnaids=as.list(org.Hs.egACCNUM[mappedkeys(org.Hs.egACCNUM)])
      xx=match(geneids$gene_id,mappedkeys(org.Hs.egACCNUM))
      geneids=geneids[!is.na(xx),2]
      rnaids=rnaids[xx[!is.na(xx)]]
    }
    ##################### Match RNAs to Gene IDs ######################
    rna2gene=c()
    for(ii in 1:length(geneids)){
      rna2gene=c(rna2gene,rep(geneids[ii],length(rnaids[ii][[1]])))
    }
    ################# Access the lengths of each transcript############
    rna2gene=cbind(rna2gene,unlist(rnaids))
    xx=match(substr(RNA@ranges@NAMES,0,nchar(RNA@ranges@NAMES)-2),rna2gene[,2])
    rna.widths=cbind(rna2gene[xx[!is.na(xx)],],RNA@ranges@width[!is.na(xx)])
    
    ############## Identify the median transcript of each gene in kb ############
    gene.transcript.length=tapply(as.numeric(rna.widths[,3]),rna.widths[,1],median)/1000
    ###### Keep only the genes that we have lengths for #########################  
    xx=match(gene.name,names(gene.transcript.length))
    gene.transcript.length=gene.transcript.length[xx[!is.na(xx)]]
    geo.data=geo.data[!is.na(xx),]
    
    ####### Divide by per million total reads and kb of transcript lengths #######
    total.reads=apply(geo.data,2,sum)/1e6
    geo.rpkms=as.matrix(geo.data)%*%diag(1/total.reads)
    transcript.length=diag(1/as.vector(gene.transcript.length))
    geo.rpkms=t(geo.rpkms)%*%transcript.length
    probe.avg=gene.avg
  }
}
}else{
  warning("Not recognized as RNA-Seq data")
}

####### Get rid of genes/probes that aren't very abundant ##########
####### Get rid of bottom quartile genes ###########################
# (renamed gene.avg to be called probe.avg)
probe.avg=apply(geo.data,1,mean)
geo.abund=geo.data[which(probe.avg> summary(probe.avg)[2]),]
geo.gene.names.abund=gene.name[which(probe.avg> summary(probe.avg)[2])]

#### If there are multiple probes per gene, use Hotelling's to find DE genes #########
if(nrow(geo.data)!=length(gene.name)){
  ######### Hotelling's MANOVA #####################################
  fstat=c()#### A bunch of empty vectors to store stats in. 
  hotelling.trace=c()
  pval=c()
  num.probes=c() # number of probes per gene
  df1=c()
  df2=c()
  for(ii in 1:length(unique(geo.gene.names.abund))){ 
    xx=which(geo.gene.names.abund==unique(geo.gene.names.abund)[ii])  ### Finds row number index for all the probes for a particular gene
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
      
      if(kappa(ss.resid)<1e-20) {## If estimated condition number for the inverse matrix is too small (which is what happens to all the p>>n scenarios; When there are so many probes per gene)
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
}else{
  fstat=c()
  pval=c()
  geo.abund=t(geo.abund)
  for(ii in 1:length(geo.gene.names.abund)){
    xx=summary(aov(geo.abund[,ii]~as.factor(groups),data=as.data.frame(geo.abund)))
    fstat[ii]=xx[[1]][["F value"]][1] ### Stores F statistic in the same vector as genes with multiple probes
    pval[ii]=xx[[1]][["Pr(>F)"]][1] ## Stores p value for 1 probe sets with the other pvalues
  }
  length(which(pval<.05)) #Before FDR correction
  pval.fdr=p.adjust(pval,method="hochberg") ##After FDR, 15 significant genes p<.05
}


#################################################################
######Put basic info in a summary file ##########################
basic.info=as.vector(c(geo.soft@header$geo_accession,geo.soft@header$summary[1],data.type,geo.soft@header$web_link,table(groups),nrow(geo.sample.info),summary(probe.avg),nrow(geo.data),sum(!is.na(gene.name)),length(unique(gene.name)),length(which(pval<.05)),length(which(pval.fdr<.05))))
names(basic.info)=c("GSE","Experiment Description:","Data/Microarray Type:","Web Link:",sort(unique(groups)),"Number of Total Samples:","Min:","1st Qu:","Median:","Mean:", "3rd Qu:","Max:","Number of probes:","Number of Probes with Gene Names:","Number of Genes:","Number of p<.05 Genes before FDR","Number of p<.05 Genes FDR")
write.table(basic.info,file=paste(data.name,"Basic Info File.txt"),quote=FALSE)

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
dimnames(geo.avg)[[1]]=geo.gene.names.abund

nn=length(unique(groups))
columns.names=c()
for(ii in 1:nn){
  tmp=matrix(ncol=nn,nrow=length(geo.gene.names.abund))
  for(jj in 1:nn){
    tmp[,jj]=geo.avg[,ii]/geo.avg[,jj]
    columns.names=c(columns.names,paste0(unique(geo.groups)[ii],"\\",unique(geo.groups)[jj]))
  }
  if(ii==1){
    geo.fc=tmp
  }else{
    geo.fc=cbind(geo.fc,tmp)
  }
}
colnames(geo.fc)=columns.names
rownames(geo.fc)=geo.gene.names.abund

########### Put fold changes into a .csv file ##############################
write.csv(geo.fc,paste(data.name,"Fold Changes.csv"),quote=FALSE)

#### Descriptive statistics for the probesets for each gene ###########
if(exists(hotelling.trace)){
  sd.probes=c()
  avg.group.signal=matrix(nrow=length(unique(geo.gene.names.abund)),ncol=length(unique(groups)))
  for(ii in 1:length(unique(geo.gene.names.abund))){
    xx=which(geo.gene.names.abund==unique(geo.gene.names.abund)[ii])
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
if(exists(hotelling.trace)){
  signif=cbind(hotelling.trace,fstat,pval,pval.fdr,num.probes,sd.probes,avg.group.signal)
  rownames(signif)=unique(geo.gene.names.abund)
  colnames(signif)=c("Hotelling Trace","F-stat","p value","FDR adj p val","Num of Probesets","SD of probsets",unique(groups))
}else{
  signif=cbind(fstat,pval,pval.fdr,avg.group.signal)
  rownames(signif)=unique(geo.gene.names.abund)
  colnames(signif)=c("F-stat","p value","FDR adj p val",unique(groups))
}

########Put the significance info into a csv file #######################
write.csv(signif[order(signif[,4]),],file=paste(data.name,"Genes.csv"))

#################### Plots #############################################
jpeg(paste(data.name,"Probes/Genes vs. -log10pvals.jpeg"))
plot(num.probes,-log10(pval.fdr))
dev.off()

