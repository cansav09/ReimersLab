######### Get p values and fold changes for differential expression analysis of Microarray and RNA-Seq data ##########
######### Determine which genes are DE across treatment groups for each dataset ########################################

setwd("/Users/cansav091/Desktop/Neurogenomics Data/MetaAnalysisSCZ")
authors=read.csv("Author-GSEKey.csv",header=FALSE)
soft.files="/Users/cansav091/Desktop/Neurogenomics Data/MetaAnalysisSCZ/SoftFiles" #Create a directory to store the soft files
gse.data="/Users/cansav091/Desktop/Neurogenomics Data/MetaAnalysisSCZ/Data" #Create a directory to store the matrix files
qc.output="/Users/cansav091/Desktop/Neurogenomics Data/MetaAnalysisSCZ/QCoutput" #Create a directory to store the QC analysis files
de.output="/Users/cansav091/Desktop/Neurogenomics Data/MetaAnalysisSCZ/DEoutput" #Create a directory to store the DE files


gpl=read.csv(file="GPLMetaDataInfo.csv")
gpl=cbind(gpl,matrix(ncol=2,nrow=nrow(gpl)))

gse.de.analysis=matrix(nrow=nrow(gpl),ncol=4)
gse.de.analysis[,1]=as.character(gpl[,1])
gse.de.analysis[,2]=as.character(gpl[,2])


if(file.exists("Groups to DE analyze by.csv")){
  xx=read.csv("Groups to DE analyze by.csv")
  gpl[,9]=xx$X1
}else{
################### Determine what groups are to be tested #######################
for(ii in 1:nrow(gpl)){
  gse.soft=paste0(soft.files,"/",gpl[ii,3])
  setwd(gse.soft)
  gse.soft=read.csv(paste0("Sample Info ",gpl[ii,2]," .csv"),skip=1,fill=TRUE,header=FALSE,colClasses="character")
  groups=rep(NA,nrow(gse.soft))
  
  headers=read.csv(paste0("Sample Info ",gpl[ii,2]," .csv"),nrows=1,fill=TRUE,header=FALSE,colClasses="character")
  headers=grep("Characteristic",headers)

  sample.char=gse.soft[,c(4,headers)]
  remov=c()
  for(kk in 1:length(ncol(sample.char))+1){
    if(length(unique(sample.char[,kk]))==1){
      remov=c(remov,kk)
    }
  }
  if(length(remov)>0){
  sample.char=sample.char[,-remov]
  } 
  if(is.vector(sample.char)==FALSE){
  keywords=as.vector(apply(sample.char,2,as.character))
  }else{
    keywords=as.vector(as.character(sample.char))
  }
  keywords=gsub("\\(","",keywords)
  keywords=gsub("\\)","",keywords)
  keywords=unlist(strsplit(keywords," "))
  keywords=unlist(strsplit(keywords,"_"))
  keywords=unlist(strsplit(keywords,"-"))
  
  numbers=paste(names(table(keywords)),"n=",table(keywords),collapse="  ")
  keywords=unique(keywords)
  print(sample.char)
  keywords=select.list(c(NA,keywords),multiple=TRUE,title=paste0("Pick out the keywords to create groups by",numbers,collapse=" "))

  if(is.na(keywords[1])){
    gse.de.analysis[ii,3]=NA
    next
  }
  for(jj in 1:length(keywords)){
    if(is.vector(sample.char)==FALSE){
    xx=which(matrix(grepl(keywords[jj],as.matrix(sample.char)),ncol=ncol(sample.char)),arr.ind=TRUE)
    groups[unique(xx[,1])]=keywords[jj]
    }else{
      xx=grep(keywords[jj],sample.char)
      groups[unique(xx)]=keywords[jj]
    }
  }
  gpl[ii,8]=paste0(groups,collapse="\t")
  xx=paste0(unique(groups)," n=",table(groups),",")
  gse.de.analysis[ii,4]=paste0(xx,collapse=" ")
}

gpl[,9]=gse.de.analysis[,4]
setwd("/Users/cansav091/Desktop/Neurogenomics Data/MetaAnalysisSCZ")
write.csv(gpl,"Groups to DE analyze by.csv",quote=TRUE)
}

###################### Differential Expression Analysis ######################

for(ii in 2:nrow(gpl)){
  de.output="/Users/cansav091/Desktop/Neurogenomics Data/MetaAnalysisSCZ/DEoutput" #Create a directory to store the DE files

  if(dir.exists(paste0(de.output,"/",gpl[ii,3]))==FALSE){
    dir.create(paste0(de.output,"/",gpl[ii,3]))
  }

  gse=gpl[ii,2]
  
  gse.dir=paste0(gse.data,"/",gpl[ii,3])
  gse.soft=paste0(soft.files,"/",gpl[ii,3])
  gse.de=paste0(de.output,"/",gpl[ii,3])
  
  gpl.name=gpl[ii,1]
  gpl.title=gpl[ii,4]
  
  groups=as.character(gpl[ii,9])
  groups=as.factor(unlist(strsplit(groups,"\t")))

  setwd(gse.dir)
  load(paste0(gse,"-",gpl.name,"matrixAllSamples.RData"))
  
    if(file.exists(paste0("GeneName-Probe",gpl.name,".csv"))){
      gene.key=read.csv(paste0("GeneName-Probe",gpl.name,".csv"))
      gene.name=gene.key$symbol[match(geo.probe.ids,gene.key$probe_id)] #This tells you in order of the queried entries,what the coordinates are of the second group of entries only  the first match is given 
      gene.name=gene.name[match(geo.probe.ids,gene.key$probe_id)]
      
      geo.data=geo.data[!is.na(gene.name),]
      gene.name=gene.name[!is.na(gene.name)]
      }

      ####### Get rid of genes/probes that aren't very abundant ##########
      ####### Get rid of bottom quartile genes ###########################
      # (renamed gene.avg to be called probe.avg)

      probe.avg=apply(geo.data,1,mean)
      cutoff=.10
      cutoff=sort(probe.avg)[nrow(geo.data)*cutoff]
  
      geo.abund=geo.data[which(probe.avg> cutoff),]
      gene.abund=gene.name[which(probe.avg> cutoff)]
      
      if(length(which(is.na(groups)))>1){
      geo.abund=geo.abund[,-which(is.na(groups))]
      groups=groups[-which(is.na(groups))]
      }
      ####### Use ANOVA for DE testing ###########################
      pval.probe=c()
      fstat.probe=c()
      for(jj in 1:length(gene.abund)){
      xx=summary(aov(geo.abund[jj,]~as.factor(groups),data=as.data.frame(geo.abund[jj,])))## Summarizes aov object into the report stats.
      fstat.probe[jj]=xx[[1]][["F value"]][1] ### Stores F statistic in the same vector as genes with multiple probes
      pval.probe[jj]=xx[[1]][["Pr(>F)"]][1]
      }
      beep(sound=2)
      pval.fdr=p.adjust(pval.probe,method="hochberg")
      #### Select probe with largest p value ###########
      
      #### Descriptive statistics for the probe/gene ###########
      geo.avg=matrix(nrow=nrow(geo.abund),ncol=length(unique(groups)))
      geo.sd=matrix(nrow=nrow(geo.abund),ncol=length(unique(groups))+1)
      for (jj in 1:length(geo.abund[,1])){
        geo.avg[jj,]=tapply(geo.abund[jj,],as.character(groups),FUN=mean)
        geo.sd[jj,]=c(tapply(geo.abund[jj,],as.character(groups),FUN=sd),sd(geo.abund[jj,]))
      }
      dimnames(geo.avg)[[1]]=gene.abund
      dimnames(geo.sd)[[1]]=gene.abund
      
      dimnames(geo.avg)[[2]]=unique(groups)
      dimnames(geo.sd)[[2]]=c(paste0(unique(groups),"SD"),"Total SD")
      
      setwd(gse.de)
      write.csv(cbind(geo.avg,geo.sd),file=paste(gpl[ii,4],gpl[ii,1],paste(unique(groups),collapse="-"),"FC and SD's.csv"),quote=FALSE)
      ########### Put fold changes into a .csv file ##############################
      
      nn=length(unique(groups))
      columns.names=c()
      comp=matrix(nrow=2,ncol=0)
      for(kk in 1:nn){
        tmp=matrix(ncol=nn,nrow=length(gene.abund))
        for(jj in 1:nn){
          tmp[,jj]=geo.avg[,kk]/geo.avg[,jj]
          columns.names=c(columns.names,paste0(unique(groups)[kk],"\\",unique(groups)[jj]))
          comp=cbind(comp,rbind(kk,jj))
        }
        if(kk==1){
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
      setwd(gse.de)
      write.csv(geo.fc,paste(gse,gpl[ii,1],paste(unique(groups),collapse="-"),"Fold Changes.csv"),quote=FALSE)
      
      ########### Put fold changes into a .csv file ##############################
      signif=cbind(fstat.probe,pval.probe,pval.fdr,geo.avg)
      rownames(signif)=gene.abund
      colnames(signif)=c("F-stat","p value","FDR adj p val",as.character(unique(groups)))
      setwd(gse.de)
      write.csv(signif,file=paste(gse,gpl[ii,1],paste0(unique(groups),collapse="-"),"ANOVA Stats.csv"))
      
  }
