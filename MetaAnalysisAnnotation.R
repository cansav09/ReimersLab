#################### Annotation ########################

library(beepr)
library(devtools)

soft.files="/Users/cansav091/Desktop/Neurogenomics Data/MetaAnalysisM1/SoftFiles" #Create a directory to store the soft files
gse.data="/Users/cansav091/Desktop/Neurogenomics Data/MetaAnalysisM1/Data" #Create a directory to store the matrix files
qc.output="/Users/cansav091/Desktop/Neurogenomics Data/MetaAnalysisM1/QCoutput" #Create a directory to store the QC analysis files
de.output="/Users/cansav091/Desktop/Neurogenomics Data/MetaAnalysisM1/DEoutput" #Create a directory to store the DE files

####################  Load in GSE numbers and info ###################
setwd("/Users/cansav091/Desktop/Neurogenomics Data/MetaAnalysisM1")


################## Create matrix that holds gpl info ###################
gpl=strsplit(as.character(meta.gse$Platform)," ")
n.gpls=c()
for(ii in 1:length(gpl)){
    n.gpls[ii]=length(gpl[[ii]])/3 # retrieve the number of GPL's each GSE has.
}
gpl=unlist(strsplit(as.character(meta.gse$Platform)," "))
gpl=cbind(gpl[seq(from=1,to=length(gpl),by=3)],gpl[seq(from=3,to=length(gpl),by=3)])
gpl=cbind(gpl,rep(GSE.list,n.gpls),rep(NA,nrow(gpl)),rep(NA,nrow(gpl)),rep(NA,nrow(gpl))) # Put this all in one dataset
colnames(gpl)=c("GPL","SampleNumber","GSE","Title","Files","Bioconductor?")


###################### Load in the normalized data #####################
gse.dataset.info=matrix(nrow=nrow(gpl),ncol=7)

for(ii in 1:nrow(gpl)){
  dir.name=grep(gpl[ii,3],gse.dir.names,value=TRUE)
  
  gse.dir=paste0(gse.data,"/",dir.name)
  gse.qc=paste0(qc.output,"/",dir.name)
  if(dir.exists(gse.qc)==FALSE){
    dir.create(gse.qc)
  }
  
  gpl.name=gpl[which(gpl[,3]==gpl[ii,3]),1]
  gpl.title=gpl[ii,4]
  
  setwd(gse.dir)
        if(length(grep("Seq",gpl.title))>0){
            setwd(paste0(gse.dir,"/",gpl[ii,3]))
            sample.files=grep("GSM",dir(),value=TRUE)
            if(length(sample.files)>0){
                if(length(grep(".gz",sample.files,value=TRUE))>0){
                    for(ll in 1:length(grep(".gz",sample.files,value=TRUE))){
                        gunzip(grep(".gz",sample.files,value=TRUE)[ll])
                    }
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
                    gene.names=select.list(cols[1:length(xx)],multiple=TRUE, title=paste0("What IDs or genesymbols do you want represented?"))
                } 
                setwd(gse.dir)
                write.csv(geo.data,file=paste0(gpl[ii,3],"-",gpl[ii,1],"Data.csv"))
                save(list=c("gene.names","sample.IDs","geo.data"),file=paste0(gpl[ii,3],"-",gpl[ii,1],".RData"),envir=.GlobalEnv)
                }else{
                    warning(dir.name," sample file ",sample.files[nn]," does not have ",nrow(xx)," rows")}
        }
        ############################## Upload data matrix ####################
        ##### When there's more than one gpl
            if(length(gpl.name)>0){
              if(length(grep(paste0(gpl[ii,1],"_series_matrix.txt"),gse.dir))>0){
                setwd(gse.dir) 
              }
              if(dir.exists(paste0(gse.dir,"/",gpl[ii,3]))){
                if(length(grep(paste0(gpl[ii,1],"_series_matrix.txt"),gse.dir))>0){
                  setwd(paste0(gse.dir,"/",gpl[ii,3]))
                }
              }
            }else{
              ##### When there's one gpl
              if(length(grep(paste0(gpl[ii,1],"_series_matrix.txt"),gse.dir))>0){
                setwd(gse.dir) 
              }
              if(dir.exists(paste0(gse.dir,"/",gpl[ii,3]))){
                if(length(grep(paste0("series_matrix.txt"),gse.dir))>0){
                  setwd(paste0(gse.dir,"/",gpl[ii,3]))
                }
              }
            }  
            matrix.files=grep("series_matrix.txt",dir(),value=TRUE)
            if(length(matrix.files)>0){
              for(kk in 1:length(matrix.files)){
              begins=which(readLines(matrix.files)=="!series_matrix_table_begin")
              geo.data=as.matrix(read.table(matrix.files,skip=begins+1,sep="\t",row.names=1,fill=T))
              geo.probe.ids=rownames(geo.data)[-nrow(geo.data)]
              sample.IDs=as.vector(t(read.table(matrix.files,skip=begins,sep="\t",nrows=1)[-1]))
              colnames(geo.data)=sample.IDs
              
              ###### Identify IDs
              gse.dataset.info[ii,7]=select.list(c("Probe IDs", "Transcript IDs", "Gene Names?"),multiple=TRUE, title=paste0("Are the above Probe IDs, Transcript IDs, or Gene Names?:",paste(head(geo.probe.ids),collapse=" ")))
              
              ###### Save Data to RDS file
              save(list=c("geo.probe.ids","sample.IDs","geo.data"),file=paste0(gpl[ii,3],"-",gpl[ii,1],"matrix.RData"),envir=.GlobalEnv)
              }
            }
        if(length(grep(paste0(gpl[ii,1],".annot"),dir()))>0){
            gpl.soft=getGEO(filename=paste0(gpl[ii,1],".annot"))
            annot=gpl.soft@dataTable@table
            colnames(annot)=gpl.soft@dataTable@columns
            
            }
            beep(sound=2)
        }beep(sound=8)
        
        
        