### Objecive: Find out which genes are differentially expressed in social stress situations

input="/Users/cansav091/Desktop/Neurogenomics Data/TungStressData/Original Data"
output="/Users/cansav091/Desktop/Neurogenomics Data/TungStressData/GEOMicroarrayScriptTung"

### Need these packages: 
library(Biobase)
library(GEOquery)
library(psych)
library(CondReg)
library(Biostrings)
library(beepr)
library(calibrate)
gse="GSE34129" ### Put GSE accession number here
data.name="Tung Stress Data Macaque"

############### Download SOFT file from GEO #########################
setwd(input)
if(file.exists(paste0(gse,".soft"))==FALSE){
getGEOfile(gse,destdir=getwd())
gunzip(paste0(gse,".soft.gz"))
geo.soft <- getGEO(filename=paste0(gse,".soft"))
}## This part may take some time depending on how many samples and how large the file is. 

if(length(grep("SuperSeries",geo.soft@header$relation))>1){
  xx=grep("SuperSeries",geo.soft@header$relation,value=TRUE)
  gses=substr(xx,17,nchar(xx))
  for(ii in 2:length(gses)){
    if(file.exists(paste0(gses[ii],".soft"))=="FALSE"){
      getGEOfile(gses[ii],destdir=input)
      gunzip(paste0(gses[ii],".soft.gz"))
    }
    assign(paste0("geo.soft",ii),getGEO(gses[ii],filename=paste0(gses[ii],".soft")),envir=.GlobalEnv)
  }
}

gpl=c()
groups=c()
######### Manually run this code for each geo.soft file ##########
  ff=1
  gpl=c(gpl,geo.soft1@header$platform_id,geo.soft2@header$platform_id)
  geo.sample.info=matrix(nrow=120,ncol=3)
  for (ii in 1:length(geo.soft1@header$sample_id)){
    geo.sample.info[ii,1:3]=c(geo.soft1@header$sample_id[ii],grep("dominance",geo.soft1@gsms[[ii]]@header$characteristics_ch1,value=TRUE),gpl[ff])
  }
    for (ii in 1:length(geo.soft2@header$sample_id)){
      geo.sample.info[(ii+100),1:3]=c(geo.soft2@header$sample_id[ii],grep("dominance",geo.soft2@gsms[[ii]]@header$characteristics_ch1,value=TRUE),gpl[ff])
    }
   }

  groups=as.factor(gsub(":","",gsub(" ","",geo.sample.info[,2])))

  ################ Retrieve GEO data matrix file########################
  data.file=paste0(gse,"-",gpl[ff],"_series_matrix.txt")
      if(file.exists(data.file)==FALSE){
      getGEO(gses[ff],destdir=getwd())
      gunzip(paste0(gse,"-",gpl[ff],"_series_matrix.txt.gz"))
      }
  
  begins=which(readLines(data.file)=="!series_matrix_table_begin")
  geo.data=as.matrix(read.table(data.file,skip=begins,sep="\t",nrows=length(readLines(data.file))-1,row.names=1,fill=T))
  geo.probe.ids=rownames(geo.data)[-1]
  sample.IDs=geo.data[1,]
  geo.data=geo.data[2:nrow(geo.data),]
  geo.data=apply(geo.data,2,as.numeric)
  rownames(geo.data)=geo.probe.ids
  
###################### Retrieving and Matching Annotation ##########
### Check that there is annotation for this at Bioconductor#########
data.type=grep("<tr valign=\"top\"><td nowrap>Title</td>" ,readLines(paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=",gpl[ff]),n=300))+1
data.type=readLines(paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=",gpl[ff]))[data.type]
  
    ################ For Microarray Data #####################################
    array.type=unlist(strsplit(data.type,">"))[2]## If this is microarray data it will download microarray annotation
    array.type=gsub("</td","",array.type)
    if(ff<3){
      array.type="illuminaHumanv4"
      library(illuminaHumanv4.db)
    }
    ######## Load the annotation for this microarray #########################
    xx=eval(parse(text=paste0(array.type,"SYMBOL")))
    probeids <- mappedkeys(xx)
    yy=as.data.frame(xx[probeids])
    gene.name=yy$symbol[match(geo.probe.ids,yy$probe_id)] #This tells you in order of the queried entries,what the coordinates are of the second group of entries only  the first match is given 
    ### Get rid of the data that don't have gene names that match to them
    if(sum(!is.na(gene.name))/length(geo.probe.ids)<.8){
     warning("Only ", print(paste0(round(sum(!is.na(gene.name))/length(geo.probe.ids),3)," of the probes have genes matching with this annotation")))
    }        
    geo.data=geo.data[!is.na(gene.name),]
    gene.name=gene.name[!is.na(gene.name)]

    #######################################################################
    ############# Get sample info from SOFT files #########################
    assign(paste0(gpl[ff],".data"),geo.data,envir=.GlobalEnv)
    assign(paste0(gpl[ff],".gene.name"),gene.name,envir=.GlobalEnv)
    assign(paste0(gpl[ff],".groups"),groups,envir=.GlobalEnv)
    assign(paste0(gpl[ff],".sample.info"),geo.sample.info,envir=.GlobalEnv)
    
    saveRDS(geo.data,file=paste0(name," ",gpl[ff]," Data.RDS"))

ff=2
  geo.data=eval(parse(text=paste0(gpl[ff],".data")))
  gene.name=eval(parse(text=paste0(gpl[ff],".gene.name")))
  groups=eval(parse(text=paste0(gpl[ff],".groups")))
  geo.sample.info=eval(parse(text=paste0(gpl[ff],".sample.info")))
  
  ####### Get rid of genes/probes that aren't very abundant ##########
  ####### Get rid of bottom quartile genes ###########################
  # (renamed gene.avg to be called probe.avg)
  probe.avg=apply(geo.data,1,mean)
  cutoff=.10
  cutoff=sort(probe.avg)[nrow(geo.data)*.10]
  
  geo.abund=geo.data[which(probe.avg> cutoff),]
  gene.abund=gene.name[which(probe.avg> cutoff)]

  assign(paste0(tolower(gpl[ff]),".abund"),geo.abund,envir=.GlobalEnv)
  assign(paste0(tolower(gpl[ff]),".gene.abund"),gene.abund,envir=.GlobalEnv)
  
  #### If there are multiple probes per gene, use Hotelling's to find DE genes #########
  if(nrow(geo.abund)>ncol(geo.abund)){
  geo.abund=t(geo.abund)
  }
  fstat=c()
  pval=c()
  groups=as.factor(groups[1:120])
  for(ii in 1:length(gene.abund)){
    xx=aov(geo.abund[,ii]~groups,data=as.data.frame(geo.abund))
    xx=summary(xx)
    fstat[ii]=xx[[1]][["F value"]][1] ### Stores F statistic in the same vector as genes with multiple probes
    pval[ii]=xx[[1]][["Pr(>F)"]][1] ## Stores p value for 1 probe sets with the other pvalues
  }
  beep()
  length(which(pval<.05)) #Before FDR correction
  pval.fdr=p.adjust(pval,method="hochberg") ##After FDR, 15 significant genes p<.05
  
  
  #################################################################
  ######Put basic info in a summary file ##########################
  basic.info=as.vector(c(geo.soft@header$geo_accession,geo.soft@header$summary[1],table(groups),nrow(geo.sample.info),summary(probe.avg)[1:6],nrow(geo.data),sum(!is.na(gene.name)),length(unique(gene.name)),length(which(pval<.05)),length(which(pval.fdr<.05))))
  names(basic.info)=c("GSE","Experiment Description:",sort(unique(groups)),"Number of Total Samples:","Min:","1st Qu:","Median:","Mean:", "3rd Qu:","Max:","Number of probes:","Number of Probes with Gene Names:","Number of Genes:","Number of p<.05 Genes before FDR","Number of p<.05 Genes FDR")
  setwd(output)
  write.table(basic.info,file=paste(data.name,gpl[ff],"Basic Info File 2.txt"),quote=FALSE)
  
  ############################################################################
  ############## Find average for each probe/gene by groups###################
  geo.avg=matrix(nrow=ncol(geo.abund),ncol=length(unique(groups)))
  for (ii in 1:length(gene.abund)){
    tmp=tapply(geo.abund[,ii],groups,FUN=mean)
    geo.avg[ii,]=tmp
  }
  dimnames(geo.avg)[[2]]=unique(groups)
  dimnames(geo.avg)[[1]]=gene.abund
  
  ############################################################################
  ########### Calculate fold changes across groups ###########################
  
  nn=length(unique(groups))
  columns.names=c()
  comp=matrix(nrow=2,ncol=0)
  for(ii in 1:nn){
    tmp=matrix(ncol=nn,nrow=length(gene.abund))
    for(jj in 1:nn){
      tmp[,jj]=geo.avg[,ii]/geo.avg[,jj]
      columns.names=c(columns.names,paste0(unique(groups[which(geo.sample.info[,3]==gpl[ff])])[ii],"\\",unique(groups)[jj]))
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
  columns.names=columns.names[-which(comp[1,]==comp[2,])]
  
  comp=comp[,-which(comp[1,]==comp[2,])]
  rownames(geo.fc)=gene.abund
  
  geo.fc[which(geo.fc<1)]=-1/geo.fc[which(geo.fc<1)]
  
  overall.fc=apply(geo.avg[,c(1,5)],1,mean)/apply(geo.avg[,c(2,4)],1,mean)
  overall.fc[which(overall.fc<1)]=-1/overall.fc[which(overall.fc<1)]
  
  ########### Put fold changes into a .csv file ##############################
  setwd(output)
  write.csv(geo.fc,paste(data.name,gpl[ff],"Fold Changes.csv"),quote=FALSE)
  
  #### Descriptive statistics for the probesets for each gene ###########
  avg.group.signal=geo.avg
  signif=cbind(fstat,pval,pval.fdr,avg.group.signal,geo.fc)
  rownames(signif)=gene.abund
  colnames(signif)=c("F-stat","p value","FDR adj p val",unique(groups),columns.names)
  
  
  ########Put the significance info into a csv file #######################
  setwd(output)
  write.csv(signif[order(signif[,4]),],file=paste(data.name,gpl[ff],"Genes.csv"))

  #### Use Regression to look at the relationships between domninance ranks #########
  tstat=c()
  pval.rg=c()
  for(ii in 1:length(gene.abund)){
    xx=summary(lm(geo.abund[,ii]~ as.numeric(substr(groups,14,14))))
    tstat[ii]=xx$coefficients[2,3] ### Stores F statistic in the same vector as genes with multiple probes
    pval.rg[ii]=xx$coefficients[2,4] ## Stores p value for 1 probe sets with the other pvalues
  }
  beep()
  length(which(pval.rg<.05)) #Before FDR correction
  pval.fdr.rg=p.adjust(pval.rg,method="hochberg") ##After FDR, 15 significant genes p<.05

  plot(-log10(pval.rg),-log10(pval))
  
  setwd(output)
  jpeg("FC vs Regression -log10pval Tung Data.jpeg",width=1500,height=1500)
  xx=cor(abs(overall.fc),-log10(pval.rg))
  plot(abs(overall.fc),-log10(pval.rg),pch="",xlab="Absolute Value for High Dominance/Low Dominance FC",cex.lab=2 )
  textxy(abs(overall.fc),-log10(pval.rg),gene.abund,offset=0,cex=1.3)
  abline(h=1.30,col="red")
  abline(v=1.30,col="blue")
  dev.off()

  
  
  
  