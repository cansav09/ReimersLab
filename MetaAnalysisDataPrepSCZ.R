#################### Annotation ########################

library(beepr)
library(devtools)

soft.files="/Users/cansav091/Desktop/Neurogenomics Data/MetaAnalysisSCZ/SoftFiles" #Create a directory to store the soft files
gse.data="/Users/cansav091/Desktop/Neurogenomics Data/MetaAnalysisSCZ/Data" #Create a directory to store the matrix files
qc.output="/Users/cansav091/Desktop/Neurogenomics Data/MetaAnalysisSCZ/QCoutput" #Create a directory to store the QC analysis files
de.output="/Users/cansav091/Desktop/Neurogenomics Data/MetaAnalysisSCZ/DEoutput" #Create a directory to store the DE files

####################  Load in GSE numbers and info ###################
setwd("/Users/cansav091/Desktop/Neurogenomics Data/MetaAnalysisSCZ")
GSE.list=c("GSE21138","GSE53987","GSE12649")

###################### Load in the normalized data #####################
gse.dataset.info=matrix(nrow=nrow(gpl),ncol=7)

for(ii in 1:nrow(gpl)){
  ii=1
    matrix.files=c()
    gse=gpl[ii,3]
    gse.dir=paste0(gse.data,"/",grep(gpl[ii,3],gse.dir.names,value=TRUE)[1])
    gse.soft=paste0(soft.files,"/",grep(gpl[ii,3],gse.dir.names,value=TRUE))
    gse.de=paste0(de.output,"/",gpl[ii,3])
    gpl.name=gpl[ii,1]
    gpl.title=gpl[ii,6]
    
    gse.files=unlist(strsplit(files.list[ii],"\n"))
    
    if(dir.exists(gse.qc)==FALSE){
        dir.create(gse.qc)
    }
    sample.files=grep("GSM",dir(),value=TRUE)
    setwd(gse.dir)
        if(length(grep("Seq",gpl.title))>0){
            setwd(paste0(gse.dir,"/",gse))
            sample.files=grep("GSM",dir(),value=TRUE)
            if(length(grep("RAW.tar",dir()))>0){
                for(ll in 1:length(grep("RAW.tar",dir()))){
                  untar(grep("RAW.tar",dir(),value=TRUE)[ll])
                }
              }
            if(length(grep(".gz",dir()))>0){
                    xx=length(grep(".gz",dir()))
                    for(ll in 1:xx){
                        gunzip(grep(".gz",dir(),value=TRUE)[ll])
                    }
                sample.files=gsub(".gz","",sample.files)
                xx=read.table(sample.files[1],fill=TRUE,colClasses="character")
                
                ll=1
                num.nas=0 # Identify at what row the data start
                while(num.nas<1){
                    ll=ll+1
                    pp=suppressWarnings(as.numeric(xx[ll,]))
                    num.nas=length(which(!is.na(pp)))
                }
                start.line=ll
                
                geo.data=read.table(sample.files[1],skip=ll-1,fill=TRUE)
                cols=as.vector(t(read.table(sample.files[1],nrow=1,skip=ll-2,colClasses="character")))
                vars=cols
                remov=select.list(vars,multiple=TRUE, title="Which variables do you NOT want in the final dataset?")
                remov=match(remov,cols)
                
                xx=read.table(sample.files[2],skip=ll,fill=TRUE)
                if(length(remov)>0){
                xx=xx[,-remov]
                geo.data=geo.data[,-remov]
                }
                if(nrow(xx)==nrow(geo.data)){
                    for(pp in 1:ncol(xx)){
                        mtches=match(as.character(geo.data[,pp]),as.character(xx[,pp]))
                        if(length(which((mtches==1:nrow(xx))=="FALSE"))<1){
                            remov=c(remov,pp)
                        }
                    }
                    remov=ifelse(length(remov)>0,remov,NULL)
                    n.cols=length(cols)-length(remov)
                    geo.data=cbind(geo.data,matrix(ncol=((length(sample.files)*n.cols)-(ncol(geo.data)-1)),nrow=nrow(geo.data)))
                    for(nn in 2:length(sample.files)){
                        yy=read.table(sample.files[nn],skip=ll,fill=TRUE)[,-remov]
                        geo.data[,((nn*n.cols-(n.cols-1))+length(remov)):((nn*n.cols)+length(remov))]=yy
                    }
                    samps=rep(paste0("Sample",2:length(sample.files)),each=n.cols)
                    vars=rep(cols[-remov])
                    cols=c(cols,paste0(samps,vars))
                    dimnames(geo.data)[[2]]=cols
                    
                    sample.IDs=sample.files
                    gene.names=select.list(cols[1:length(xx)],multiple=TRUE, title="What IDs or genesymbols do you want represented?")
                
                setwd(gse.dir)
                write.csv(geo.data,file=paste0(gse,"-",gpl.name,"Data.csv"))
                save(list=c("gene.names","sample.IDs","geo.data"),file=paste0(gse,"-",gpl.name,"matrix.RData"),envir=.GlobalEnv)
                }else{ 
                warning(gse.dir.names[ii]," sample file ",sample.files[nn]," does not have ",nrow(xx)," rows")}
                }
        }
      
        ############################## Upload data matrix ####################
      if(length(grep(".gz",sample.files,value=TRUE))>0){
        for(ll in 1:length(grep(".gz",sample.files,value=TRUE))){
          gunzip(grep(".gz",sample.files,value=TRUE)[ll])
        }
      }
            
      if(length(grep(gpl[ii,3],gpl[,3]))>1){
              if(length(grep(paste0(gpl.name,"_series_matrix.txt"),dir(gse.dir)))>0){
                setwd(gse.dir) 
                matrix.files=grep(paste0(gpl.name,"_series_matrix.txt"),dir(gse.dir),value=TRUE)
              }
              if(dir.exists(paste0(gse.dir,"/",gse))){
                if(length(grep(paste0(gpl.name,"_series_matrix.txt"),dir(paste0(gse.dir,"/",gse))))>0){
                  setwd(paste0(gse.dir,"/",gse))
                  matrix.files=grep(paste0(gpl.name,"_series_matrix.txt"),dir(paste0(gse.dir,"/",gse)),value=TRUE)
                }
              }
            }else{
              ##### When there's one gpl
              if(length(grep("_series_matrix.txt",dir(gse.dir)))==1){
                setwd(gse.dir) 
                matrix.files=grep("series_matrix.txt",dir(gse.dir),value=TRUE)
                }
              }
            if(length(matrix.files)>0){
              for(kk in 1:length(matrix.files)){
              begins=which(readLines(matrix.files)=="!series_matrix_table_begin")
              geo.data=as.matrix(read.table(matrix.files,skip=begins+1,sep="\t",row.names=1,fill=T))
              if(dim(geo.data)[1]<2){
                next()
              }
              geo.probe.ids=rownames(geo.data)[-nrow(geo.data)]
              sample.IDs=as.vector(t(read.table(matrix.files,skip=begins,sep="\t",nrows=1)[-1]))
              colnames(geo.data)=sample.IDs
              
              gse.dataset.info[ii,7]=ifelse(length(grep("_at",geo.probe.ids))>10,"Affy Probes",gse.dataset.info[ii,7])
              gse.dataset.info[ii,7]=ifelse(length(grep("ILMN",geo.probe.ids))>10,"Illumina Probes",gse.dataset.info[ii,7])
              gse.dataset.info[ii,7]=ifelse(length(grep("A_51_P100021",geo.probe.ids))>10,"Agilent Probes",gse.dataset.info[ii,7])
              gse.dataset.info[ii,7]=ifelse(length(grep("NM_",geo.probe.ids))>10,"RefSeq Transcripts",gse.dataset.info[ii,7])
              gse.dataset.info[ii,7]=ifelse(length(grep("ENST",geo.probe.ids))>10,"Ensemble Transcripts",gse.dataset.info[ii,7])
              gse.dataset.info[ii,7]=ifelse(length(grep("1070000",geo.probe.ids))>6,"Affy ST Probes",gse.dataset.info[ii,7])
              gse.dataset.info[ii,7]=ifelse(length(which(!is.na(match(c("actb","gapdh","cd8","bdnf"),tolower(geo.probe.ids)))))>0,"Gene Names",gse.dataset.info[ii,7])
              
              if(is.na(gse.dataset.info[ii,7])){
              gse.dataset.info[ii,7]=select.list(c("Probe IDs", "Transcript IDs", "Gene Names","NA"),multiple=TRUE, title=paste("Are these Probe IDs, Transcript IDs, or Gene Names?",paste(head(geo.probe.ids),collapse=" ")))
              }
              ###### Identify IDs, load in the GeneName-Probe key if it exists 
              if(length(grep(paste0("GeneName-Probe",gpl.name),dir(gse.dir)))>0){
                setwd(gse.dir)
                gene.names=read.table(grep(paste0("GeneName-Probe",gpl.name),dir(gse.dir),value=TRUE))
              ###### Save Data to RDS file
                save(list=c("gene.names","geo.probe.ids","sample.IDs","geo.data"),file=paste0(gse,"-",gpl.name,"matrix.RData"),envir=.GlobalEnv)
              }else{
              ###### Save Data to RDS file
              save(list=c("geo.probe.ids","sample.IDs","geo.data"),file=paste0(gse,"-",gpl.name,"matrix.RData"),envir=.GlobalEnv)
              }
            }
        if(length(grep(paste0(gpl.name,".annot"),dir()))>0){
            gpl.soft=getGEO(filename=paste0(gpl.name,".annot")[1])
            annot=gpl.soft@dataTable@table
            colnames(annot)=gpl.soft@dataTable@columns
            save(list=c("annot","geo.probe.ids","sample.IDs","geo.data"),file=paste0(gse,"-",gpl.name,"matrix.RData"),envir=.GlobalEnv)
          }
        }
  beep(sound=2)
}
if(ii!=length(GSE.list)){beep(sound=9)}else{beep(sound=5)}

        
        