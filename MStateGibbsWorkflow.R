##### Goal: See how the gibbs data looks using which ever gene list you choose to put in.
## Run martinezDataGEOfile Code before doing this. 
## This workflow extracts the data for all the probes that match the gene list you give it and then uses the
## one with the highest SD to represent the gene and then does PCA

name="MStateGibbsWorkflow"
mainDir=paste0("/Users/cansav091/Desktop/Neurogenomics Data/martinez m0m1m2 Lists/",name)
input="/Users/cansav091/Desktop/Neurogenomics Data/martinez m0m1m2 Lists/Original Data"
outDir=paste0(mainDir,"/",name)
outDirDay=paste0(outDir,"/",as.character(Sys.Date()))
if(dir.exists(outDir)==FALSE){
    dir.create(outDir)
}
if(dir.exists(outDirDay)==FALSE){
    dir.create(outDirDay)
}

library(calibrate)
library(mASS)
library(Biobase)
library(GEOquery)
# Import Gibbs et al data if it hasn't been imported yet
if(exists("qnormed")==FALSE){
load("Gibbs.cortex.qnormed.Rdata")
}
#Import annotation for Gibbs data if it hasn't been imported yet
setwd(input)
gse="GSE15745"
if(file.exists("GSE15745.soft")==FALSE){
getGEOfile(gse,destdir=getwd())
gunzip(grep("GSE15745.soft.gz",dir(),value=TRUE))
geo.soft.gibbs <- getGEO(filename="GSE15745.soft")
}

setwd(input)
mstategenes=read.csv("m1 Genes Hotel p<.05 FC>1.3.csv")[,1]
mstategenes=intersect(mstategenes,gene.name)
warning("Name the m State what you'd like it to be called in output")
mstate="m1"


##### Retrieving all columns of data for all probes that correspond to the gene lists#####
mstategib=matrix(nrow=291,ncol=0)#Need an empty vector to put the data I want in
mstategib.gene.name=c()#Need an empty vector for the gene names that correspond to the matrix
for (i in 1:length(mstategenes)){
    nn=which(gene.name==mstategenes[i])#gives indices for gene.name that match the names from m1gene names
    mstategib=cbind(mstategib,qnormed[,nn])#Add the correspondin columns to the growing list of the dataset
    mstategib.gene.name=c(mstategib.gene.name,gene.name[nn])#Create a vector with the gene names for the new dataset
}
#Saving these data to a text file
setwd(outDirDay)
write.table(mstategib,file=paste(mstate,"Probes Gibbs Data",date()) ,quote=FALSE,sep="\t")

mstatepca=prcomp(mstategibb,scale=TRUE)
cor.mstate=cor(mstatepca$x[,1],mstategib)

##Plot the loadings to look at the distribution of the loadings throughout the geneset
jpeg(paste(mstate,"Geneset PCA Loading Distribution Gibbs",date()),width=700,height=600)
plot(mstatepca$rotation,cex.lab=1.5)
textxy(mstatepca$rotation[,1],mstatepca$rotation[,2],mstategib.gene.name)## Delete this if you don't want the gene names on the plot
dev.off()

setwd(outDir)
########### Download Annotation for Gene(s) of Interest ######### 
gene.family=c("GABA","Glutamate")
for(jj in 1:length(gene.family)){
  HUGO=readLines("http://www.genenames.org/cgi-bin/genefamilies/")
  

##################### Edit list ##############################
hugo.lists=HUGO[grep(gene.family[jj],HUGO)]
hugo.lists=unlist(strsplit(hugo.lists,c("\"")))
sub.cat=hugo.lists[seq(from=3, to=length(hugo.lists),by=3)]
sub.cat=substr(sub.cat,2,nchar(sub.cat)-4)
hugo.lists=hugo.lists[seq(from=2,length(hugo.lists),by=3)]

##### Download gene family information from HUGO website ###################################
gene.family.info=data.frame()
for(ii in 1:length(hugo.lists)){
  download.file(paste0("http://www.genenames.org/cgi-bin/genefamilies/",hugo.lists[ii],"/download/branch") ,destfile=paste0("HUGO Gene Family Info for ",sub.cat[ii],".txt"))
  xx=read.table(paste0("HUGO Gene Family Info for ",sub.cat[ii],".txt"),fill=TRUE,sep="\t")[-1,]
  gene.family.info=rbind(gene.family.info,xx)
}
colnames(gene.family.info)=c("HGNC ID","Approved Symbol","Approved Name","Status","Previous Symbols","Synonyms","Chromosome","Accession Numbers","RefSeq IDs","Gene Family Tag","Gene family","ID")

################## PCA for gene family of interest ####################################
gene.family.data=matrix(nrow=291,ncol=0)
gene.names.oi=as.character(gene.family.info$`Approved Symbol`)
genes.oi.loadings=c()
gene.oi.match=c()
gene.fam.probes=c()

for (ii in 1:length(gene.names.oi)){
  nn=grep(paste0("^",gene.names.oi[ii],"$"),gene.name)
  gene.family.data=cbind(gene.family.data,qnormed[,nn])
  gene.oi.match=c(gene.oi.match,gene.name[nn])
  gene.fam.probes=c(gene.fam.probes,rep(gene.names.oi[ii],length(nn)))
} 
means=apply(gene.family.data,2,mean)
# Found the mean difs for each sample 
mean.difs=apply(gene.family.data,1,function(x) x-means)
# Found the covariance of each 
covmat=cov(gene.family.data)
# Found the Eigenvalues for the covariance matrix
eigs=eigen(covmat)

assign(paste0(gene.family[jj],".pca"),prcomp(gene.family.data),envir = .GlobalEnv)
assign(paste0(gene.family[jj],".Gibbs"),gene.family.data,envir = .GlobalEnv)
assign(paste0(gene.family[jj],".gene.names"),gene.oi.match,envir = .GlobalEnv)
assign(paste0(gene.family[jj],".eigs"),eigs,envir = .GlobalEnv)

}


eigs.sort=eigs$vectors[,order(eigs$values,decreasing=TRUE)]


write.csv(cbind(Glutamate.gene.names,Glutamate.pca$rotation[,1],Glutamate.eigs$values,Glutamate.eigs$vectors[,order(Glutamate.eigs$values,decreasing=TRUE)][,1],apply(Glutamate.Gibbs,2,sd)),file="Glutamate Gibbs Loadings.csv")
write.csv(cbind(GABA.gene.names,GABA.pca$rotation[,1],GABA.eigs$values,GABA.eigs$vectors[,order(GABA.eigs$values,decreasing=TRUE)][,1],apply(GABA.Gibbs,2,sd)),file="GABA Gibbs Loadings.csv")








