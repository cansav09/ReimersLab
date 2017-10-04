## Demographic downloader


setwd("/Users/cansav091/Desktop/Neurogenomics Data/MetaAnalysisSCZ")
GSE.list=c("GSE12649","GSE21138","GSE53987")# Load in a list of desired GSE's
authors=read.csv("Author-GSEKey.csv",header=FALSE)
gpl=as.data.frame(read.csv("GPLMetaDataInfo.csv"))[,-1]

for(ii in 1:length(GSE.list)){
  ii=1
  
setwd(paste0("/Users/cansav091/Desktop/Neurogenomics Data/MetaAnalysisSCZ/SoftFiles/",gpl$GSE.Dir[ii]))
samples=read.csv(grep("Sample Info",dir(),value=TRUE))
sample.IDs=samples$GSM.ID

for(jj in 1:length(sample.IDs)){
web.info=readLines(paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=",sample.IDs[jj]))

start=grep("Characteristics",web.info)
web.info=web.info[start:(start+100)]
stop=grep("tr valign",web.info)

web.info=web.info[1:stop[1]+1]
web.info=unlist(strsplit(web.info,"<br>"))
web.info[1]=substr(web.info[1],33,nchar(web.info[1]))

if(jj==1){
  first.row=t(web.info[-length(web.info)])
  demographic=matrix(ncol=length(first.row),nrow=length(sample.IDs))
  rownames(demographic)=sample.IDs
  
  xx=unlist(strsplit(first.row,":"))
  if(length(xx)>ncol(first.row)){
  colnames(demographic)=xx[seq(from=1,to=length(xx),by=2)]
  demographic[1,]=xx[seq(from=2,to=length(xx),by=2)]
  }else{
    demographic[1,]=first.row
  }

}else{
  row.info=t(web.info[-length(web.info)]) 
  if(length(xx)>length(first.row)){
  demographic[jj,]=unlist(unlist(strsplit(row.info,":")))[seq(from=2,to=length(xx),by=2)]
  }else{
    demographic[jj,]=row.info
  }
}

}
saveRDS(demographic,file=paste0(GSE.list[ii],"DemographicData.RDS"))
setwd("/Users/cansav091/Desktop/Neurogenomics Data/MetaAnalysisSCZ")
saveRDS(demographic,file=paste0(GSE.list[ii],"DemographicData.RDS"))

}



demo1=readRDS("GSE53987DemographicData.RDS")
demo2=readRDS("GSE21138DemographicData.RDS")
demo3=readRDS("GSE12649DemographicData.RDS")


