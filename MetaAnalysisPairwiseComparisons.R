###### Meta Analysis Pairwise comparisons

gpl=read.csv("GPLDataInfoGroupsEdit.csv")[,-1]
gpl=cbind(gpl,matrix(ncol=2,nrow=nrow(gpl)))


datasets=c()
datasetspvals=c()
for(ii in 1:nrow(gpl)){
  gse.de=paste0(de.output,"/",gpl[ii,4])
  if(dir.exists(gse.de)){
  setwd(gse.de)
  filename=grep(gpl[ii,1],dir(),value=TRUE)
  if(length(filename)==1){
    obj=unlist(strsplit(grep(gpl[ii,1],filename,value=TRUE)," "))[1]
    assign(tolower(obj),read.csv(grep("FC and SD's.csv",filename,value=TRUE)),envir=.GlobalEnv)
    assign(paste0(tolower(obj),"pval"),read.csv(grep("ANOVA",filename,value=TRUE)),envir=.GlobalEnv)
    datasets=c(datasets,tolower(obj))
    datasetspvals=c(datasetspvals,paste0(tolower(obj),"pval"))
  }
  if(length(filename)>1){
  obj=unlist(strsplit(grep(gpl[ii,1],filename,value=TRUE)," "))[1]
  assign(paste0(tolower(obj),ii),read.csv(grep("FC and SD's.csv",filename,value=TRUE)),envir=.GlobalEnv)
  assign(paste0(tolower(obj),ii,"pval"),read.csv(grep("ANOVA",filename,value=TRUE)),envir=.GlobalEnv)
  datasets=c(datasets,paste0(tolower(obj),ii))
  datasetspvals=c(datasetspvals,paste0(tolower(obj),ii,"pval"))
    }
  }
}



for(ii in 1:length(datasets)){
  gse1pval=eval(parse(text=paste0(datasets[ii],"pval")))
  gse1=eval(parse(text=datasets[ii]))
  
  assign(paste0(datasetspval[ii]),tapply(gse1pval$p.value,gse1pval$X,min),envir=.GlobalEnv)
  
  tapply(gse1pval$p.value,gse1pval$X,min)
  assign(paste0(),tapply(gse1pval$p.value,gse1pval$X,min),envir=.GlobalEnv)
  
}






for(ii in 1:length(datasets)){
  gse1pval=eval(parse(text=paste0(datasets[ii],"pval")))
  gse1=eval(parse(text=datasets[ii]))
    for(jj in 1:length(datasets)){
      if(ii==jj){next}
      gse2pval=eval(parse(text=paste0(datasets[ii],"pval")))
      gse2=eval(parse(text=datasets[jj]))
      
      
      xx=match(gse1$X,gse2$X)
      gse1=gse1[is.na(xx),]
      gse1=gse1[xx[is.na(xx)],]
      
    }
  }









