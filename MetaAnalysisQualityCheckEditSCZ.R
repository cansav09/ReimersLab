#################################### 1) Platform Info Sorting ##########################
#################################### 2) Gene name to probe id key download #############
#################################### 3) Create GC% and MA plots ########################

library(beepr)
library(devtools)
library(installr)
library(seqinr)
library(affxparser)
library(affy)
library(affyPLM)
library(oligo)
library(oligoData)
library(ggplot2)
library(lumi)
library(illuminaio)
library(limma)

soft.files="/Users/cansav091/Desktop/Neurogenomics Data/MetaAnalysisSCZ/SoftFiles" #Create a directory to store the soft files
gse.data="/Users/cansav091/Desktop/Neurogenomics Data/MetaAnalysisSCZ/Data" #Create a directory to store the matrix files
qc.output="/Users/cansav091/Desktop/Neurogenomics Data/MetaAnalysisSCZ/QCoutput" #Create a directory to store the QC analysis files
de.output="/Users/cansav091/Desktop/Neurogenomics Data/MetaAnalysisSCZ/DEoutput" #Create a directory to store the DE files
GSE.list=c("GSE21138","GSE53987","GSE12649")

if(file.exists("GPLMetaDataInfo.csv")){
  setwd("/Users/cansav091/Desktop/Neurogenomics Data/MetaAnalysisSCZ")
  gpl=as.data.frame(read.csv("GPLMetaDataInfo.csv"))[,-1]
  

}else{
  setwd("/Users/cansav091/Desktop/Neurogenomics Data/MetaAnalysisSCZ")
  authors=read.csv("Author-GSEKey.csv",header=FALSE)
  for(ii in 1:length(GSE.list)){
    xx=authors[grep(GSE.list[ii],authors$V1),2]
    gse.dir.names[ii]=paste0(xx," et al- ",GSE.list[ii])
  }
  ############## Create a matrix to hold the platform data ##############
  gpl=strsplit(as.character(meta.gse$Platform)," ")
  n.gpls=c()
  for(ii in 1:length(gpl)){
    n.gpls[ii]=length(gpl[[ii]])/3 # retrieve the number of GPL's each GSE has.
  }
  gpl=unlist(gpl)
  gpl=cbind(gpl[seq(from=1,to=length(gpl),by=3)],gpl[seq(from=3,to=length(gpl),by=3)])
  gpl=cbind(gpl,rep(GSE.list,n.gpls),rep(NA,nrow(gpl)),rep(NA,nrow(gpl)),rep(NA,nrow(gpl)),rep(NA,nrow(gpl))) # Put this all in one dataset
  colnames(gpl)=c("GPL","SampleNumber","GSE","Title","GPL.Annot?","BioconductorPlatform","BioconductorAnnot")
  
  #### Prep bioconductor annotation package list
  web=readLines("https://www.bioconductor.org/packages/release/data/annotation/")
  web=web[-c(1:117)]
  xx=grep("<td><a href=\"html/",web)
  
  bioc.packs=strsplit(as.character(web[xx]),"html")
  bioc.packs=unlist(bioc.packs)[seq(from=2,to=length(unlist(bioc.packs)),by=3)]
  bioc.packs=substr(bioc.packs,2,nchar(bioc.packs)-1)
  names(bioc.packs)=gsub("</td>","",gsub("            <td>","",web[xx+2]))
  
  ########### Run through all of the datasets and gather their GPL titles and whether or not they have GPL.annot files###############
  for(ii in 1:nrow(gpl)){
    gse.dir=paste0(gse.data,"/",gpl[ii,3])
    gse.qc=paste0(qc.output,"/",gse.dir.names[ii])
    if(dir.exists(gse.qc)==FALSE){
      dir.create(gse.qc)
    }
    gpl.name=gpl[ii,1]
    
    gpl.title=c()
    for(jj in 1:nrow(gpl)){
      xx=readLines(paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=",gpl.name))
      gpl.title=xx[grep("Title",xx)+1][1]
      gpl.title=gsub("<td style=\"text-align: justify\">","",gpl.title)
      gpl.title=gsub("</td>","",gpl.title)
      gpl[ii,4]=gpl.title
      ##### Identify where annotation is available
      gpl[ii,5]=ifelse(length(grep(paste0(gpl.name,".annot"),dir(gse.data[ii]),value=TRUE))>0,"Yes","No")
    }
    beep(sound=2)
  } 
  if(ii!=nrow(gpl)){beep(sound=9)}else{beep(sound=5)}
  
  ################################################################################################ 
  ################################ Download Microarray Annotation ################################
  
  terms=tolower(c("Affymetrix","Illumina","Mouse","Human","ST","430","u133","ht12","1.0","2.0"))
  for(ii in 1:nrow(gpl)){
    gse.dir=paste0(gse.data,"/",grep(gpl[ii,3],gse.dir.names,value=TRUE))
    gpl.name=gpl[ii,1]
    gpl.title=gpl[ii,4]
    
    search=terms[!is.na(match(terms,unlist(strsplit(tolower(gpl.title),c(" ",".")))))]
    
    if(length(grep("Seq",gpl.title))<1){
      if(length(grep("Affy",gpl.title))>0){
        ######################## Load Affy Platform Info ##################
        if(is.na(gpl[ii,6])){
          pds=select.list(c("Expand List",bioc.packs[grep(paste0(search,collapse="|"),names(bioc.packs))]), title=paste0("Pick out Platform Design Package: ",gpl.title))
          if(pds==""){
            pds=gpl[ii,6]="None Found"
          }else{
            if(pds=="Expand List"){
              pds=select.list(bioc.packs,multiple=TRUE, title=paste0("Pick out Platform Design Package: ",gpl.title))
              if(pds==""){
                pds=gpl[ii,6]="None Found"
              }
            }
            if(length(pds)>0){
              pds=gpl[ii,6]=bioc.packs[match(pds,bioc.packs)]
            }
          }
        }else{
          pds=gpl[ii,6]
        }
        if(!is.na(match(gpl.title,gpl[-(ii),4]))){
          gpl[which(gpl[,4]==gpl.title),6]=gpl[ii,6]
        }
        if(pds!="None Found"){
          if(length(grep(pds,installed.packages()))<1){
            source("https://bioconductor.org/biocLite.R")
            biocLite(pds)
          }
          library(grep(pds,installed.packages(),value=TRUE)[1],character.only=TRUE)
        }
      }
      ######################## Load gene - probe id key ##################
      if(is.na(gpl[ii,7])){
        gene.ids=select.list(c("Expand List",bioc.packs[grep(paste0(search,collapse="|"),names(bioc.packs))]), title=paste0("Pick out Gene Annotation Package: ",gpl.title))
        if(gene.ids==""){
          gene.ids=gpl[ii,7]="None Found"
        }else{
          if(gene.ids=="Expand List"){
            gene.ids=select.list(bioc.packs,multiple=TRUE, title=paste0("Pick out Gene Annotation Package: ",gpl.title))
            if(gene.ids==""){
              gene.ids=gpl[ii,7]="None Found"
            }
          }
          if(length(gene.ids)>0){
            gene.ids=gpl[ii,7]=bioc.packs[match(gene.ids,bioc.packs)]
          }
          gene.ids=gpl[ii,7]=bioc.packs[match(gene.ids,bioc.packs)]
        }
      }else{
        gene.ids=gpl[ii,7]
      }
      if(!is.na(match(gpl.title,gpl[-(ii),4]))){
        gpl[which(gpl[,4]==gpl.title),7]=gpl[ii,7]
      }
      if(gene.ids!="None Found"){
        if(length(grep(gene.ids,installed.packages()))<1){
          source("https://bioconductor.org/biocLite.R")
          biocLite(gene.ids)
        }
        gene.ids=grep(gene.ids,installed.packages(),value=TRUE)[1]
        library(gene.ids,character.only=TRUE)
        xx=eval(parse(text=paste0(gsub(".db","",gene.ids),"SYMBOL")))
        
        probeids <- mappedkeys(xx)
        gene.name.key=as.data.frame(xx[probeids])
        
        ######################## Save key to csv file #########################
        setwd(gse.dir)
        write.csv(gene.name.key,file=paste0("GeneName-Probe",gpl.name,".csv"),quote=FALSE)
      }
    }
    beep(sound=2)
  } 
  if(ii!=nrow(gpl)){beep(sound=9)}else{beep(sound=5)}
  
  #### Write the package info into the GPL meta info file ########
  gpl=c(gpl[,1:3],rep(NA,nrow(gpl)),gpl[,4:ncol(gpl)])
  for(ii in 1:nrow(gpl)){
    gpl[ii,4]=grep(gpl[ii,3],gse.dir.names,value=TRUE)[1]
  }
  
  setwd("/Users/cansav091/Desktop/Neurogenomics Data/MetaAnalysisSCZ")
  write.csv(gpl,file="GPLMetaDataInfo.csv",quote=FALSE)
}

################################ For Microarray Data:
################################ Normalize data, Check for GC bias, Make MA Plots ################################
for(ii in 1:nrow(gpl)){
  ii=2
  dir.name=gpl[ii,4]
  gse=gpl[ii,3]
  gse.dir=paste0(gse.data,"/",dir.name)
  gse.qc=paste0(qc.output,"/",dir.name)
  if(dir.exists(gse.qc)==FALSE){
    dir.create(gse.qc)
  }
  
  gpl.name=gpl[ii,1]
  gpl.title=gpl[ii,5]
  
  ###################### Remove Bad Samples ######################
  setwd(gse.qc)
  keep=NULL
  if(file.exists(grep("QualityControl",dir()))==TRUE){
  synop=read.csv(grep("QualityControl",dir(),value=TRUE))
  keep=synop$GSM.ID[-grep("Remove",synop$Conclusion)]
  }
    ################################ Affymetrix Microarrays ################################
      if(dir.exists(paste0(gse.dir,"/",gpl[ii,3]))){
        for(mm in 1:length(region)){
         setwd(paste0(gse.dir,"/",gpl[ii,3]))
      }
      if(!is.na(paste0(gpl[ii,7]))){
        library(paste0(gpl[ii,7]),character.only=TRUE)
        library(paste0(gpl[ii,8]),character.only=TRUE)
      }
        cel.files=grep(".CEL",dir(),value=TRUE)
        region=NULL
        
        xx=c()
        for(kk in 1:(length(keep))){
        xx=c(xx,grep(keep[kk],cel.files))
        }
        
        cel.files=cel.files[xx]
        setwd(paste0(gse.dir,"/",gpl[ii,3]))
        if(ii==2){
        region=c("STR","HPC","B46")
        cel.files=grep(region[mm],cel.files,value=TRUE)
        }  
        ##### Load in .CEL files and get the signal intensities for the probes
        rawData=read.celfiles(cel.files)
        cel.data=read.celfiles(cel.files)
        signals=oligo::intensity(rawData)
        samp=gsub(".CEL.gz","",dimnames(signals)[[2]])
        
        #### Get the probe sequences ################ 
        pmSeq=pmSequence(rawData)
        pmSeqlog2=log2(pm(rawData))
        
        #### Calculate GC Content ################  
        zz=round(runif(25000,min=0,max=length(pmSeq)))
        affy.seq=c()
        affy.gc=c()
        for(kk in 1:length(zz)){
          bb=zz[kk]
          affy.seq=c(affy.seq,as.character(pmSeq[[bb]]))
          affy.gc=c(affy.gc,GC(unlist(strsplit(affy.seq[[kk]],""))))
        }
        
        ### Create the jpegs for the cel images ########################
        if(dir.exists(paste0(gse.qc,"/CEL images",region[mm]))==FALSE){
          dir.create(paste0(gse.qc,"/CEL images",region[mm]))
        }
        setwd(paste0(gse.qc,"/CEL images",region[mm]))
        samp=gsub(".gz","",cel.files)
        for (kk in 1:length(cel.files)){
          jpeg(paste0("image",cel.files[kk],".jpg"))
          image(rawData[,kk],main=samp[kk])
          dev.off()
        }
        
        ### Plot GC content vs relative signal
        if(dir.exists(paste0(gse.qc,"/GC ContentOutliersRemoved",region[mm]))==FALSE){
          dir.create(paste0(gse.qc,"/GC ContentOutliersRemoved",region[mm]))
        }
        setwd(paste0(gse.qc,"/GC ContentOutliersRemoved",region[mm]))
        avg.sig=apply(pmSeqlog2,1,mean)
        for(kk in 1:length(cel.files)){
          jpeg(paste0("log2 PM Probes vs GC content ",samp[kk],".jpeg"))
          plot(jitter(affy.gc,factor=2.3),pmSeqlog2[zz,kk]-avg.sig[zz],pch=".",xlab="GC Ratio",ylab="log2(signal)-avg(log2(signal))",ylim=c(-1,1))
          abline(h=0,col="red")
          dev.off()
        }
  
        ########### Calculate trimmed means #####################
        # This step takes a little while
        log.standard = construct.log.standard(intensity(cel.data))
        
        ########### Create Bias display images ######################
        # This step takes a little while as well. 
        setwd(paste0(gse.dir,"/",gpl[ii,3]))

        if(dir.exists(paste0(gse.qc,"/BiasDisplayOutliersRemoved",region[mm]))==FALSE){
          dir.create(paste0(gse.qc,"/BiasDisplayOutliersRemoved",region[mm]))
        }
        signals=oligo::intensity(rawData)
        setwd(paste0(gse.qc,"/BiasDisplayOutliersRemoved",region[mm]))
        for(kk in 1:length(cel.files)){
          jpeg(paste0("Bias Display Plot for sample ",samp[kk],gpl.name,".jpeg"))
          bias.display(signals[,kk],log.standard,mask=1:length(log.standard))
          dev.off()
        }
    ########### Raw Data MA Plots ######################
    if(dir.exists(paste0(gse.qc,"/MA PlotsOutliersRemoved",region[mm]))==FALSE){
      dir.create(paste0(gse.qc,"/MA PlotsOutliersRemoved",region[mm]))
    }
    setwd(paste0(gse.qc,"/MA PlotsOutliersRemoved",region[mm]))
      for(kk in 1:ncol(signals)){
      jpeg(paste0("MA Plot for sample",samp[kk],".jpeg"))
      yy=log2(signals[zz,kk])-log.standard[zz]
      plot(log.standard[zz],yy,main=paste0("MA plot for ",samp[kk]),pch=".",xlab="trim mean of log2(signal)",ylab="log2(signal)-trim mean of log2(signal)",ylim=c(-1.5,1.5))
      abline(h=0, col = "blue")
      dev.off()
        }
    beep(sound=2)
      }
    }





for(ii in 1:nrow(gpl)){
  dir.name=gpl[ii,4]
  gse=gpl[ii,3]
  gse.dir=paste0(gse.data,"/",dir.name)
  gse.qc=paste0(qc.output,"/",dir.name)
  gse.soft=paste0(soft.files,"/",dir.name)
  gpl.name=gpl[ii,1]
  gpl.title=gpl[ii,5]
  
  file.name=unlist(strsplit(as.character(dir.name)," "))[1]
  ###################### Remove Bad Samples ######################
  setwd(gse.qc)
  synop=read.csv(grep("QualityControl",dir(),value=TRUE))
  keep=synop$GSM.ID[-grep("Remove",synop$FinalConclusion)]
  
    ################################ Affymetrix Microarrays ################################
        setwd(paste0(gse.dir,"/",gpl[ii,3]))
        library(paste0(gpl[ii,7]),character.only=TRUE)
        library(paste0(gpl[ii,8]),character.only=TRUE)
  
  for(mm in 1:length(region)){  
        setwd(paste0(gse.dir,"/",gpl[ii,3]))
        cel.files=grep(".CEL",dir(),value=TRUE)
        xx=c()
        for(kk in 1:(length(keep))){
          xx=c(xx,grep(keep[kk],cel.files))
        }
        cel.files=cel.files[xx]
        setwd(paste0(gse.dir,"/",gpl[ii,3]))
          if(ii==2){
            region=c("STR","HPC","B46")
            cel.files=grep(region[mm],cel.files,value=TRUE)
          }
        ##### Load in .CEL files and get the signal intensities for the probes
        rawData=read.celfiles(cel.files)
        cel.data=read.celfiles(cel.files)
        signals=oligo::intensity(rawData)
        samp=gsub(".CEL.gz","",dimnames(signals)[[2]])
        
        ##### Load in .CEL files and get the signal intensities for the probes
        rma.data=rma(rawData)

        setwd(gse.dir)
        gene.names=read.csv(grep(".csv",dir(),value=TRUE))
        setwd(gse.soft)
        sample.info=read.csv(grep("Sample Info",dir(),value=TRUE))
        if(ii==2){
          sample.info=sample.info[grep(tolower(region[mm]),sample.info[,3]),]
        }
        setwd(gse.qc)
        saveRDS(c(gene.names,sample.info,rma.data),paste0(file.name,gse,region[mm],"RMANormalized.RDS"))
  }
#################### REDO QC but with normalized data##############        
        
        ########### Calculate trimmed means #####################
        # This step takes a little while
        rma.data=exprs(rma.data)
        log.standard = construct.log.standard(rma.data,trim=.1)
        
        zz=round(runif(25000,min=0,max=nrow(rma.data)))
        ########### Raw Data MA Plots ######################
        if(dir.exists(paste0(gse.qc,"/NormalizedMAPlotsOutliersRemoved",region[mm]))==FALSE){
          dir.create(paste0(gse.qc,"/NormalizedMAPlotsOutliersRemoved",region[mm]))
        }
        setwd(paste0(gse.qc,"/NormalizedMAPlotsOutliersRemoved",region[mm]))
        for(kk in 1:ncol(rma.data)){
          jpeg(paste0("Normalized MA Plot for sample",samp[kk],".jpeg"))
          yy=log2(rma.data[zz,kk])-log.standard[zz]
          plot(log.standard[zz],yy,main=paste0("MA plot for ",samp[kk]),pch=".",xlab="mean of log2(signal)",ylab="log2(signal)-trim mean of log2(signal)")
          abline(h=0, col = "blue")
          dev.off()
        }
      }
    }


