#### Meta Analysis Dataset Downloading
#### Downloads the SOft Files and Matrix Data (or Raw files) for a list of GSE accession numbers

#### Necessary libraries:
library(GEOquery)
library(beepr)

setwd("/Users/cansav091/Desktop/Neurogenomics Data/MetaAnalysisSCZ")
GSE.list=c("GSE12649","GSE21138","GSE53987")# Load in a list of desired GSE's
authors=read.csv("Author-GSEKey.csv",header=FALSE)

soft.files="/Users/cansav091/Desktop/Neurogenomics Data/MetaAnalysisSCZ/SoftFiles" #Create a directory to store the soft files
gse.data="/Users/cansav091/Desktop/Neurogenomics Data/MetaAnalysisSCZ/Data" #Create a directory to store the matrix files
qc.output="/Users/cansav091/Desktop/Neurogenomics Data/MetaAnalysisSCZ/QCoutput" #Create a directory to store the QC analysis files
de.output="/Users/cansav091/Desktop/Neurogenomics Data/MetaAnalysisSCZ/DEoutput" #Create a directory to store the DE files

if(dir.exists(soft.files)==FALSE){
    dir.create(soft.files)
    dir.create(gse.data)
    dir.create(qc.output)
    dir.create(de.output)
    dir.create(paste0(soft.files,"/","Original Soft Files"))
}

soft.info=matrix(nrow=length(GSE.list),ncol=12)
gse.dir.names=c()
for(ii in 1:length(GSE.list)){
  
  ############### Download SOFT file from GEO #########################
  setwd(paste0(soft.files,"/","Original Soft Files"))
  if(file.exists(paste0(GSE.list[ii],".soft"))==FALSE){
    getGEOfile(GSE.list[ii],destdir=paste0(soft.files,"/","Original Soft Files"))
    gunzip(paste0(GSE.list[ii],".soft.gz"))
  }
  geo.soft <- getGEO(filename=paste0(GSE.list[ii],".soft"))
  
  xx=authors[grep(GSE.list[ii],authors$V1),2]
  gse.dir.names[ii]=paste0(xx," et al- ",GSE.list[ii])
  soft.info[ii,2]=as.character(unlist(xx)[length(xx)])
  
  gse.dir=paste0(soft.files,"/",xx," et al- ",GSE.list[ii])
  gse.dir.names[ii]=paste0(xx," et al- ",GSE.list[ii])
  
  if(dir.exists(gse.dir)==FALSE){
    dir.create(gse.dir)
  }
  
  ## Extract info to put in the metadata matrix (meta as in all the datasets)
  soft.info[ii,1]=geo.soft@header$title
  soft.info[ii,2]=as.character(unlist(xx)[length(xx)])
  soft.info[ii,4]=geo.soft@header$submission_date
  soft.info[ii,9]=ifelse(length(geo.soft@header$pubmed_id)>0,geo.soft@header$pubmed_id,NA)
  soft.info[ii,10]=ifelse(length(geo.soft@header$supplementary_file)>0,"Yes","No")
  soft.info[ii,11]=length(geo.soft@gsms)
  soft.info[ii,8]=ifelse(length(geo.soft@gsms[[1]]@dataTable@table$ID_REF)==0,NA,length(geo.soft@gsms[[1]]@dataTable@table$ID_REF))
  
  setwd(gse.dir)
  if(length(geo.soft@header$supplementary_file)>0){
    write.csv(geo.soft@header$supplementary_file,file=paste("Supp Files urls",GSE.list[ii],".csv"))
  }
  
  ### Take out the info specific to the dataset
  n.char=length(geo.soft@gsms[[1]]@header$characteristics_ch1)
  xx=n.char/length(geo.soft@gsms)
  if(n.char>0){
    xx=n.char/length(geo.soft@gsms) # characteristics per sample
    if(xx==round(xx)){
      geo.sample.info=matrix(nrow=length(geo.soft@header$sample_id),ncol=9+xx)
    }else{
      n.char.samp=c()
      for(jj in 1:length(geo.soft@header$sample_id)){
        n.char.samp=c(n.char.samp,length(geo.soft@gsms[[jj]]@header$characteristics_ch1))
      }
      n.char=max(n.char.samp)
      geo.sample.info=matrix(nrow=length(geo.soft@header$sample_id),ncol=9+n.char)
    }
  }else{
    geo.sample.info=matrix(nrow=length(geo.soft@header$sample_id),ncol=10)
  }
  geo.sample.info[,1]=geo.soft@header$sample_id
  geo.variables=c("organism_ch1","title","label_protocol_ch1","hyb_protocol_ch1","extract_protocol_ch1","treatment_protocol_ch1","data_processing","platform_id")
  
  for(ll in 1:(length(geo.variables))){
    if(length(eval(parse(text=paste0("geo.soft@gsms[[1]]@header$",geo.variables[ll]))))==1){
      for(jj in 1:length(geo.soft@header$sample_id)){
        geo.sample.info[jj,ll+1]=eval(parse(text=paste0("geo.soft@gsms[[jj]]@header$",geo.variables[ll])))[1]
      }
    }
    if(length(eval(parse(text=paste0("geo.soft@gsms[[1]]@header$",geo.variables[ll]))))>1){
      for(jj in 1:length(geo.soft@header$sample_id)){
        geo.sample.info[jj,ll+1]=eval(parse(text=paste0("geo.soft@gsms[[1]]@header$",geo.variables[ll])))[jj]
      }
    }
    if(length(eval(parse(text=paste0("geo.soft@gsms[[1]]@header$",geo.variables[ll]))))<1){
      geo.sample.info[,ll+1]=NA
    }
  }
  if(n.char>0){
    if(xx==round(xx)){
      for(jj in 1:length(geo.soft@header$sample_id)){
        xx=n.char/length(geo.soft@gsms)
        geo.sample.info[jj,10:(9+xx)]=geo.soft@gsms[[1]]@header$characteristics_ch1[((jj*xx)-(xx-1)):(jj*xx)]
      }
    }else{
      for(jj in 1:length(geo.soft@header$sample_id)){
        if(n.char.samp[jj]>0){
          geo.sample.info[jj,10:(9+n.char)]=ifelse(n.char.samp[jj]==n.char,geo.soft@gsms[[jj]]@header$characteristics_ch1,c(geo.soft@gsms[[jj]]@header$characteristics_ch1,rep(NA,n.char-n.char.samp[jj])))
        }
      }
    }
  }else{
    geo.sample.info[jj,10]=NA
  }
  
  if(xx==round(xx)){
    n.char=xx
  }
  n.char=ifelse(n.char<1,1,n.char)
  
  ######### Create a Sample Info Matrix for this dataset
  setwd(gse.dir)
  colnames(geo.sample.info)=c("GSM ID","Organism","Title","Label Protocol","Hyb Protocol","Extraction Protocol","Treatment Protocol","Data Processing","Platform",paste(rep("Characteristic",n.char),1:n.char))
  geo.sample.info=gsub(",","",geo.sample.info)
  write.csv(geo.sample.info,file=paste("Sample Info",GSE.list[ii],".csv"),quote=FALSE)
  
  ###### Get the N's for the different sample variables
  sample.char=table(geo.sample.info[,10:(9+n.char)])
  org.tab=table(geo.sample.info[,2])
  data.tab=table(geo.sample.info[,8])
  gpl.tab=table(geo.sample.info[,9])
  
  ##### Create a crosstabs of the treatment groups
  if(n.char>1){
    write.csv(sample.char,file=paste("Sample Characteristics",GSE.list[ii]),quote=FALSE,row.names=FALSE)
    if(n.char==2){
      cross.tabs=xtabs(~geo.sample.info[,9]+geo.sample.info[,10], data=geo.sample.info[,9:(8+n.char)])
      write.csv(cross.tabs,file=paste("Xtabs of Treatment groups",GSE.list[ii],".csv"),quote=FALSE)
    }
    if(n.char==3){
      cross.tabs=xtabs(~geo.sample.info[,9]+geo.sample.info[,10]+geo.sample.info[,11], data=geo.sample.info[,9:(8+n.char)])
      write.csv(cross.tabs,file=paste("Xtabs of Treatment groups",GSE.list[ii],".csv"),quote=FALSE)
    }
  }
  soft.info[ii,3]=ifelse(length(gpl.tab)>0,paste(names(gpl.tab),"n=",gpl.tab,collapse=" "),NA)
  soft.info[ii,6]=ifelse(length(sample.char)>0,paste(names(sample.char),"n=",sample.char,collapse=" "),NA)
  soft.info[ii,5]=paste(names(org.tab),"n=",org.tab,collapse=" ")
  soft.info[ii,7]=paste(names(data.tab),"n=",data.tab,collapse=" ")
  beep(sound=2)
}
if(ii!=length(GSE.list)){beep(sound=9)}else{beep(sound=5)}


colnames(soft.info)=c("Title","First Author","Platform","Date","Species","Sample Char","Data Type","Total Probe or Gene Count","PubMed Ref","Supp Files?","Total Sample N","FilesList")
rownames(soft.info)=GSE.list
soft.info=gsub(",","",soft.info)

setwd("/Users/cansav091/Desktop/Neurogenomics Data/MetaAnalysisSCZ")
write.csv(soft.info,file="Meta Dataset Info Table2.csv",quote=FALSE)
if(file.exists("Meta Dataset Info Table.csv")){beep(sound=8)}

######################################################################################################
################## Download the data for each corresponding GSE Accession Number #####################
setwd("/Users/cansav091/Desktop/Neurogenomics Data/MetaAnalysisSCZ")
meta.gse=read.csv("Meta Dataset Info Table2.csv")
meta.gse$Total.Probe.or.Gene.Count[which(meta.gse$Total.Probe.or.Gene.Count=="RNA-Seq")]=NA

################## Sort the GPLs by the GSE numbers #####################

gpl=strsplit(as.character(meta.gse$Platform)," ")
n.gpls=c()
for(ii in 1:length(gpl)){
  n.gpls[ii]=length(gpl[[ii]])/3 # retrieve the number of GPL's each GSE has.
}
gpl=unlist(strsplit(as.character(meta.gse$Platform)," "))
gpl=cbind(gpl[seq(from=1,to=length(gpl),by=3)],gpl[seq(from=3,to=length(gpl),by=3)])
gpl=cbind(gpl,rep(GSE.list,n.gpls),rep(NA,nrow(gpl)),rep(NA,nrow(gpl)),rep(NA,nrow(gpl))) # Put this all in one dataset


################## Download the files by GSE number #####################

for(ii in 1:length(gse.dir.names)){#### repeats for each dataset
  gse.dir=paste0(gse.data,"/",gse.dir.names[ii])
  if(dir.exists(gse.dir)==FALSE){ #creates a directory to store the normalized and non normalized data
    dir.create(gse.dir)
  }
  setwd(gse.dir)
  gpl.name=gpl[which(gpl[,3]==GSE.list[ii]),1]
  
  if(length(gpl.name)==1){ # If the dataset only uses one GPL...
    if(file.exists(paste0(gpl.name,".annot"))==FALSE){
      xx=readLines(paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=",gpl.name))
      if(length(xx[grep(paste0(gpl.name,".annot"),xx)])>0){
        suppressWarnings(download.file(url=paste0("ftp://ftp.ncbi.nlm.nih.gov/geo/platforms/",substr(gpl.name,1, nchar(gpl.name)-3),"nnn/",gpl.name,"/annot/",gpl.name,".annot.gz"),destfile=paste0(gpl.name,".annot.gz")))
        if(length(grep(paste0(gpl.name,".annot"),dir()))>0){
          gunzip(paste0(gpl.name,".annot.gz"))
        }
      }
    }
    if(!is.na(meta.gse$Total.Probe.or.Gene.Count[ii])){ ## Essentially if it is microarray
      if(file.exists(paste0(GSE.list[ii],"_series_matrix.txt"))==FALSE){
        matrixfile=paste0(GSE.list[ii],"_series_matrix.txt.gz")
        download.file(url=paste0("ftp://ftp.ncbi.nlm.nih.gov/geo/series/",substr(GSE.list[ii],1, nchar(GSE.list[ii])-3),"nnn/",GSE.list[ii],"/matrix/",matrixfile),destfile=matrixfile)
        gunzip(matrixfile)
      }
      if(meta.gse$Supp.Files.[ii]=="Yes"){
        if(dir.exists(paste0(gse.dir,"/",GSE.list[ii]))==FALSE){
          getGEOSuppFiles(GSE.list[ii]) # Download Raw data
          setwd(paste0(setwd(gse.dir),"/",GSE.list[ii]))
          
          if(length(grep("_RAW.tar",dir()))>0){
            untar(paste0(GSE.list[ii],"_RAW.tar"))
            chp.files=grep(".chp.gz",dir(),value=TRUE)
            cel.files=grep(".CEL.gz",dir(),value=TRUE)
            if(length(cel.files)>1){ # unzip any CEL or CHP files
              for(kk in 1:length(cel.files)){
                gunzip(cel.files[kk])
              }
              if(length(chp.files)>1){
                for(kk in 1:length(cel.files)){
                  gunzip(chp.files[kk])
                }
              }
            }
          }
        }
      }
    }else{ ### If it is RNA-seq data
      if(meta.gse$Supp.Files.[ii]=="Yes"){
        if(dir.exists(paste0(gse.dir,"/",GSE.list[ii]))==FALSE){
          getGEOSuppFiles(GSE.list[ii]) # Download Supp files
          setwd(paste0(setwd(gse.dir),"/",GSE.list[ii]))
          
          if(length(grep(".gz",dir()))>0){  # unzip files
            data.files=grep(".gz",dir(),value=TRUE)
            for(kk in 1:length(data.files)){
              gunzip(data.files[kk])
            }
            data.files=substr(data.files,1,nchar(data.files)-3)
          }
          if(length(grep("_RAW.tar",dir()))>0){ # untar files
            untar(paste0(GSE.list[ii],"_RAW.tar"))
          }
        }
      }
    }
  }
  ############If the GSE has only one GPL for it ##########
  if(length(gpl.name)>1){
    for(kk in 1:length(gpl.name)){
      if(file.exists(paste0(gpl.name[kk],".annot"))==FALSE){############ Download GPL annotation file ##########
        xx=readLines(paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=",gpl.name[jj]))
        if(length(xx[grep(paste0(gpl.name[kk],".annot"),xx)])>0){
          download.file(url=paste0("ftp://ftp.ncbi.nlm.nih.gov/geo/platforms/",substr(gpl.name[kk],1, nchar(gpl.name[kk])-3),"nnn/",gpl.name[kk],"/annot/",gpl.name[kk],".annot.gz"),destfile=paste0(gpl.name[kk],".annot.gz"))
          gunzip(paste0(gpl.name[kk],".annot.gz"))
        }
      }
      
      if(as.numeric(meta.gse$Total.Probe.or.Gene.Count[ii])>0){
        if(file.exists(paste0(GSE.list[ii],"-",gpl.name[kk],"_series_matrix.txt"))==FALSE){
          matrixfile=paste0(GSE.list[ii],"-",gpl.name[kk],"_series_matrix.txt.gz")
          download.file(url=paste0("ftp://ftp.ncbi.nlm.nih.gov/geo/series/",substr(GSE.list[ii],1, nchar(GSE.list[ii])-3),"nnn/",GSE.list[ii],"/matrix/",matrixfile),destfile=matrixfile)
          gunzip(matrixfile)
        }else{
          if(meta.gse$Supp.Files.[ii]=="Yes"){
            if(dir.exists(paste0(gse.dir,"/",GSE.list[ii]))==FALSE){
              getGEOSuppFiles(GSE.list[ii])
              setwd(paste0(setwd(gse.dir),"/",GSE.list[ii]))
              
              if(length(grep(".gz",dir()))>0){
                data.file=grep(".gz",dir(),value=TRUE)
                gunzip(data.file)
                substr(data.file,1,nchar(data.file)-3)
              }
              if(length(grep("_RAW.tar",dir()))>0){
                untar(paste0(GSE.list[ii],"_RAW.tar"))
                if(length(grep(".gz",dir()))>0){
                  data.files=grep(".gz",dir(),value=TRUE)
                  for(ll in 1:length(data.files)){
                    gunzip(data.files[ll])
                  }
                  data.files=substr(data.files,1,nchar(data.files)-3)
                }
              }
            }
          }
        }
      }
    }
  }
  beep(sound=2)
}
if(ii!=length(GSE.list)){beep(sound=9)}else{beep(sound=5)}

files.list=matrix(ncol=length(GSE.list),nrow=2)
files.list[1,]=GSE.list
for(ii in 1:length(gse.dir.names)){ #### repeats for each dataset
  gse.dir=paste0(gse.data,"/",gse.dir.names[ii])
  setwd(gse.dir)
  files=paste0(dir(),collapse="\n")
  if(meta.gse$Supp.Files.[ii]=="Yes"){
    if(dir.exists(paste0(gse.dir,"/",GSE.list[ii]))==FALSE){
      getGEOSuppFiles(GSE.list[ii])
    }
    setwd(paste0(gse.dir,"/",GSE.list[ii]))
    files=paste0(files,paste0(GSE.list[ii],"/",dir(),collapse="\n"),collpase="\n")
    files.list[2,ii]=files
  }
}

setwd("/Users/cansav091/Desktop/Neurogenomics Data/MetaAnalysisM1")
write.table(files.list,file="DatasetFilesList.txt",quote=FALSE)

save.image(file = "MetaAnalysisSCZ.RData")
