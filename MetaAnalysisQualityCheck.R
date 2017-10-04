#################################### 1) Platform Info Sorting ##########################
#################################### 2) Gene name to probe id key download #############
#################################### 3) Create GC% and MA plots ########################


############# Be sure to have loaded the "SmudgeMiner" function before running this code #######################
library(beepr)
library(devtools)
library(seqinr)
library(affxparser)
library(affy)
library(affyPLM)
library(oligo)
library(oligoData)
library(ggplot2)
library(lumi)
library(illuminaio)

soft.files="/Users/cansav091/Desktop/Neurogenomics Data/MetaAnalysisM1/SoftFiles" #Create a directory to store the soft files
gse.data="/Users/cansav091/Desktop/Neurogenomics Data/MetaAnalysisM1/Data" #Create a directory to store the matrix files
qc.output="/Users/cansav091/Desktop/Neurogenomics Data/MetaAnalysisM1/QCoutput" #Create a directory to store the QC analysis files
de.output="/Users/cansav091/Desktop/Neurogenomics Data/MetaAnalysisM1/DEoutput" #Create a directory to store the DE files

setwd("/Users/cansav091/Desktop/Neurogenomics Data/MetaAnalysisM1")
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

if(file.exists("GPLMetaDataInfo.csv")){
  setwd("/Users/cansav091/Desktop/Neurogenomics Data/MetaAnalysisM1")
  gpl=as.data.frame(read.csv("GPLMetaDataInfo.csv"))[,-1]
}else{

########### Run through all of the datasets and gather their GPL titles and whether or not they have GPL.annot files###############
for(ii in 1:nrow(gpl)){
gse.dir=paste0(gse.data,"/",gpl[ii,3])
gse.qc=paste0(qc.output,"/",gpl[ii,3])
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
for(ii in 8:nrow(gpl)){
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
            gpl[which(gpl[,4]==gpl.title),6:7]=gpl[ii,6:7]
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
        gpl[which(gpl[,4]==gpl.title),6:7]=gpl[ii,6:7]
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
gpl=cbind(gpl[,1:3],rep(NA,nrow(gpl)),gpl[,4:ncol(gpl)])
for(ii in 1:nrow(gpl)){
  gpl[ii,4]=grep(gpl[ii,3],gse.dir.names,value=TRUE)
}

setwd("/Users/cansav091/Desktop/Neurogenomics Data/MetaAnalysisM1")
write.csv(gpl,file="GPLMetaDataInfo.csv",quote=FALSE)
}

################################ For Microarray Data:
################################ Normalize data, Check for GC bias, Make MA Plots ################################
for(ii in 1:nrow(gpl)){
  dir.name=gpl[ii,4]
  gse=gpl[ii,3]
  gse.dir=paste0(gse.data,"/",dir.name)
  gse.qc=paste0(qc.output,"/",dir.name)
  if(dir.exists(gse.qc)==FALSE){
    dir.create(gse.qc)
  }

  gpl.name=gpl[ii,1]
  gpl.title=gpl[ii,5]
  gpl.title
    if(length(grep("Seq",gpl.title))<1){  
################################ Affymetrix Microarrays ################################
    if(length(grep("Affy",gpl.title))>0){
        if(dir.exists(paste0(gse.dir,"/",gpl[ii,3]))){
            setwd(paste0(gse.dir,"/",gpl[ii,3]))
          }
        if(!is.na(paste0(gpl[ii,7]))){
          library(paste0(gpl[ii,7]),character.only=TRUE)
          library(paste0(gpl[ii,8]),character.only=TRUE)
        }
            if(length(grep(".CEL",dir()))>0){
              cel.files=grep(".CEL",dir(),value=TRUE)
              
              if(ii==3){
              cel.files=grep("chipA.CEL",dir(),value=TRUE)
              }
              if(ii==4){
                cel.files=grep("chipB.CEL",dir(),value=TRUE)
              }
                ##### Load in .CEL files and get the signal intensities for the probes
                cel.files=unique(gsub(".gz","",cel.files))
                rawData=read.celfiles(cel.files)
                if(length(grep("ST",gpl.title))>0){
                  cel.data=read.celfiles(cel.files)
                  signals=oligo::intensity(cel.data)
                }else{
                  cel.data = ReadAffy(filenames=cel.files)
                }
                signals=oligo::intensity(rawData)
                samp=gsub(".CEL","",dimnames(signals)[[2]])
                
                
                ###### Save the Raw Data to a .csv #########
                setwd(gse.qc)
                write.csv(signals,paste0("Raw Data",dir.name,gpl.name,".csv"),quote=FALSE)
                    
                jpeg(paste0("Signal Distribution Boxplots non-normalized affy ",gpl.name,".jpeg"))
                boxplot(signals)
                dev.off()
                    
                #### Get the probe sequences ################ 
                pmSeq=pmSequence(rawData)
                   
                #### Calculate GC Content ################  
                affy.seq=c()
                affy.gc=c()
                for(kk in 801033:length(pmSeq)){
                affy.seq=c(affy.seq,as.character(pmSeq[[kk]]))
                affy.gc=c(affy.gc,GC(unlist(strsplit(affy.seq[[kk]],""))))
                }
                beep()
                
                jpeg(paste0("Probe GC% distribution",GSE.list[ii],gpl.name,".jpeg"))
                  hist(affy.gc,xlab="GC ratio")
                dev.off()
                
                ### Create the jpegs for the cel images ########################
                if(dir.exists(paste0(gse.qc,"/CEL images"))==FALSE){
                  dir.create(paste0(gse.qc,"/CEL images"))
                }
                setwd(paste0(gse.qc,"/CEL images"))
                for (kk in 1:length(cel.files)){
                jpeg(paste0("image",cel.files[kk],".jpg"))
                  image(rawData[,kk],main=rownames(ph@data)[kk])
                dev.off()
                }
                    
                ##### Normalize the Data 
                setwd(gse.qc)
                rma.data=rma(rawData)
                write.csv(rma.data,paste0("RMA Normalized Data",dir.name,paste0(gpl[ii,c(1,3)],collapse="-"),".csv"),quote=FALSE)
                
                ##### Plot normalized data 
                jpeg(paste0("Signal Distribution Boxplots RM affy",GSE.list[ii],paste0(gpl[ii,1],collapse="-"),".jpeg"))
                  boxplot(rma.data,las=2)
                dev.off()
                    
                ### Plot GC content vs relative signal
                if(dir.exists(paste0(gse.qc,"/GC Content"))==FALSE){
                  dir.create(paste0(gse.qc,"/GC Content"))
                }
                setwd(paste0(gse.qc,"/GC Content"))
                avg.sig=apply(pmSeqlog2,1,mean)
                for(kk in 1:length(cel.files)){
                jpeg(paste0("log2 PM probes vs GC content affy Sample",kk,gpl.name,".jpeg"))
                  plot(jitter(affy.gc,factor=2.3),pmSeqlog2[,kk]/avg.sig,pch=".",xlab="GC ratio",ylab="Log2 Signal/vg Log2 Signal")
                  abline(h=1,col="red")
                dev.off()
                }
                beep()
                ########### Calculate trimmed means #####################
                # This step takes a little while
                trimmed.means=c()
                nn=round(ncol(signals)*.05)
                nn2=ncol(signals)-nn
                for(ll in 1:nrow(signals)){
                  trimmed=signals[ll,-c(which(order(signals[ll,])>=nn2),which(order(signals[ll,])<=nn))]
                  trimmed.means[ll]=mean(trimmed)
                }
                beep()
                log.standard = log2(trimmed.means)
                
                
                jpeg(paste0("Histogram Log2 Trimmed Means",dir.name,gpl.name,".jpeg"))
                hist(log2(trimmed.means))
                dev.off()
                
                ########### Create Bias display images ######################
                # This step takes a little while as well. 
                setwd(paste0(gse.dir,"/",gpl[ii,3]))
                
                if(dir.exists(paste0(gse.qc,"/BiasDisplay"))==FALSE){
                  dir.create(paste0(gse.qc,"/BiasDisplay"))
                }
                
                setwd(paste0(gse.qc,"/BiasDisplay"))
                for(kk in 1:length(cel.files)){
                jpeg(paste0("Bias Display Plot for sample ",samp[kk],gpl.name,".jpeg"))
                bias.display(signals[,kk],log.standard,mask=1:length(log.standard))
                dev.off()
                }
            }
          }
################################## Illumina Microarrays ##########################
    if(length(grep("Illumina",gpl.title))>0){
      if(dir.exists(paste0(gse.dir,"/",GSE.list[ii]))){
        setwd(paste0(gse.dir,"/",GSE.list[ii]))
        if(length(grep("_RAW.tar",dir(),value=TRUE))>0){
          untar(grep("_RAW.tar",dir(),value=TRUE))
        }
        if(length(grep(".gz",dir(),value=TRUE))>0){
          for(kk in 1:length(grep(".gz",dir()))){
            gunzip(grep(".gz",dir(),value=TRUE)[kk])
            }
          }
        
        if(length(grep(".bgx",dir()))){
          bgx=readBGX(grep(".bgx",dir(),value=TRUE))
          setwd(gse.dir)
          gene.name.key=cbind(bgx$probes$Probe_Id,bgx$probes$Probe_Id)
          write.csv(gene.name.key,file=paste0("GeneName-Probe",gpl.name[jj],".csv"),quote=FALSE)
          setwd(paste0(gse.dir,"/",GSE.list[ii]))
          }
      }
      setwd(paste0(gse.dir,"/",GSE.list[ii]))
      if(length(grep("non_normalized",dir()))>0){
      headers=as.vector(read.table(grep("non_normalized",dir(),value=TRUE),nrow=1))
      detec.cols=grep("detect",headers,value=TRUE)
      probes="ID_REF"
      signals=grep("detect",headers[-1],invert=TRUE,value=TRUE)
      rawData=lumiR(grep("non_normalized",dir(),value=TRUE),columnNameGrepPattern=list(detection=detec.cols,exprs=sigs),annotationColumn="ID_REF")
      samp=colnames(rawData) 
          }
        }
      ########### Raw Data MA Plots ######################
      if(dir.exists(paste0(gse.qc,"/MA Plots"))==FALSE){
        dir.create(paste0(gse.qc,"/MA Plots"))
        }
        setwd(paste0(gse.qc,"/MA Plots"))

      for(kk in 1:ncol(signals)){
        jpeg(paste0("MA Plot for sample",samp[kk],".jpeg"))
        zz=round(runif(20000,min=1,max=nrow(signals)))
        xx=(log2(signals[zz,kk])+log2(trimmed.means[zz]))/2
        yy=log2(signals[zz,kk]/trimmed.means[zz])
        plot(xx,yy,main=paste0("MA plot for ",samp[kk]),pch=".",xlab="(log2(signal)+log2(trim mean signal))/2",ylab="log2(signal/trim mean signal)")
        abline(h=0, col = "blue")
        dev.off()
      }
        beep(sound=2)
    }
  }
if(ii!=length(GSE.list)){beep(sound=9)}else{beep(sound=5)}










