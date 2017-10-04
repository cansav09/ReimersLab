





#################################### Preliminary Quality Checks ##########################

gse.dir=paste0(gse.data,"/",gse.dir.names[ii])
gse.de=paste0(de.output,"/",gse.dir.names[ii])

geo.data=matrix(nrow=0,ncol=)
gene.name=c()
################ This loop will repeat for each platform type here ##########
################ Retrieve GEO data matrix file########################
setwd(input)
if(length(grep(gpl[jj],dir(),value=TRUE))<1){
  getGEO(gse,destdir=getwd())
  gunzip(grep("txt.gz",dir(),value=TRUE))
}
data.file=grep(gpl[jj],dir(),value=TRUE)
data.file=grep("seriesmatrix.txt",data.file,value=TRUE)
begins=which(readLines(data.file)=="!series_matrix_table_begin")
geo.data=as.matrix(read.table(data.file,skip=begins+1,sep="\t",row.names=1,fill=T,))
geo.probe.ids=rownames(geo.data)[-nrow(geo.data)]
sample.IDs=as.vector(t(read.table(data.file,skip=begins,sep="\t",nrows=1)[-1]))
colnames(geo.data)=sample.IDs


fnames=dir()
if(length(fnames)>1){
  geo.data = matrix(nr=as.numeric(meta.gse$Total.Probe.or.Gene.Count[ii]),nc=length(fnames))
  xx= read.table(fnames[1])
  rownames(geo.data)<-xx[,1]
  for ( kk in 1:length(fnames)) {
    tmp= read.table(fnames[kk])
    geo.data[,kk] <- xx[,2] 




