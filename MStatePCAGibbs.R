#See how the Gibbs et al data look using M1 Genes identified from Martinez et al data
name="MStatePCAGibbs"
outDir=paste0(mainDir,"/",name)
outDirDay=paste0(mainDir,"/",name,"/",as.character(Sys.Date()))
outDirGraphsDay=paste0(outDir,"/Graphs","/",as.character(Sys.Date()))
if(dir.exists(outDir)==FALSE){
  dir.create(outDir)
}
if(dir.exists(outDirDay)==FALSE){
  dir.create(outDirDay)
  dir.create(outDirGraphsDay)
}

#Need these libraries
library(calibrate)
library(MASS)

load("Gibbs.cortex.qnormed.Rdata")
Mstategenes=read.csv("M1 Genes Hotel p<.05 FC>1.3.csv")[,1] # Import M1 genes
Mstategenes=intersect(Mstategenes,gene.name) # find which M1 genes are represented in Gibbs data
Mstate="M1" # Just for titles of jpegs and etc. 

##### Retrieving all columns of data for all probes that correspond to the gene lists for m1, M20, and M12
Mstategib=matrix(nrow=291,ncol=0)#Need an empty vector to put the data I want in
Mstategib.gene.name=c()#Need an empty vector for the gene names that correspond to the matrix
for (i in 1:length(Mstategenes)){
    nn=which(gene.name==Mstategenes[i])#gives indices for gene.name that match the names from m1gene names
    Mstategib=cbind(Mstategib,qnormed[,nn])#Add the correspondin columns to the growing list of the dataset
    Mstategib.gene.name=c(Mstategib.gene.name,gene.name[nn])#Create a vector with the gene names for the new dataset
}
#Saving these data to a text file
write.table(Mstategib,file=paste(Mstate,"Probes Gibbs Data",date()) ,quote=FALSE,sep="\t")

Mstatepca=prcomp(Mstategibb,scale=TRUE) 
cor.Mstate=cor(Mstatepca$x[,1],Mstategib) 

##Plot the loadings to look at the distribution of the loadings throughout the geneset
jpeg(paste(Mstate,"Geneset PCA Loading Distribution Gibbs",date()),width=700,height=600)
plot(Mstatepca$rotation,cex.lab=1.5)
textxy(Mstatepca$rotation[,1],Mstatepca$rotation[,2],Mstategib.gene.name)## Delete this if you don't want the gene names on the plot
dev.off()


########### Download Annotation for Gene of Interest ######### 
HUGO=readLines("http://www.genenames.org/cgi-bin/genefamilies/")
gene.family="GABA"

##################### Edit list ##############################
hugo.lists=HUGO[grep(gene.family,HUGO)]
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
write.table(gene.family.info,file=paste0(gene.family," Overall Gene Family Info HUGO.txt"))


################## PCA for gene family of interest ####################################
gene.family.data=matrix(nrow=291,ncol=0)
gene.names.oi=as.character(gene.family.info$`Approved Symbol`)
genes.oi.loadings=c()
gene.oi.match=c()
gene.fam.probes=c()

for (ii in 1:length(gene.names.oi)){
  nn=grep(paste0("^",gene.names.oi[ii],"$"),gene.name)
  gene.family.data=cbind(gene.family.data,qnormed[,nn])
  gene.oi.match=c(gene.oi.match,gene.names.oi[nn])
  gene.fam.probes=c(gene.fam.probes,rep(gene.names.oi[ii],length(nn)))
} 
gene.fam.pca=prcomp(gene.family.data,scale=TRUE)
cor.gene=cor(gene.fam.pca$x[,1],gene.family.data)

jpeg(paste(Mstate,gene.family,"PC scatterplot.jpeg"))
hist(cor.gene,main="Distribution of correlations")
dev.off()

####### Testing correlations with randomly selected housekeping genes with M1 PCA #############
house.genes=read.table("HousekeepingGenes.txt",sep="\t")
house.genes=as.character(house.genes$V1)
house.genes=substr(house.genes,1,nchar(house.genes)-1)
house.data=matrix(nrow=291,ncol=0)
house.gene.match=c() 
house.probes=c()

#for(jj in 1:100){
  #house.gene.sel=house.genes[floor(runif(15,min=1,max=length(house.genes)))]
  for (ii in 1:length(house.genes)){
    nn=grep(paste0("^",house.gene.sel[ii],"$"),gene.name)
    house.data=cbind(house.data,qnormed[,nn])
    house.gene.match=c(house.gene.match,house.gene.sel[nn])
    house.probes=c(house.probes,rep(house.gene.sel[ii],length(nn)))
  } 
  #house.pca=prcomp(house.data,scale=TRUE)
  #cor.house[jj]=cor(house.pca$x[,1],Mstatepca$x[,1])
#}
house.pca=prcomp(house.data,scale=TRUE)
jpeg(paste(Mstate,"and Housekeeping Genes","PC scatterplot.jpeg"))
cor.house=cor(house.pca$x[,1],house.data)
plot(house.pca$x[,1],Mstatepca$x[,1],ylab=paste0(Mstate,"PC"),xlab= paste("Housekeeping PC"),sub=paste("R=",round(cor(house.pca$x[,1],Mstatepca$x[,1]),3)),cex.lab=1.5)
dev.off()

################## Testing correlations with randomly selected genes and M1 PCA#############
cor.rand=c()
xx=rep((1:10)^2,10)

for(jj in 1:100){
rand.gene.sel=gene.name[floor(runif(xx[jj],min=1,max=length(gene.name)))]
rand.data=matrix(nrow=291,ncol=0)
for (ii in 1:length(rand.gene.sel)){
  nn=grep(paste0("^",rand.gene.sel[ii],"$"),gene.name)
  rand.data=cbind(rand.data,qnormed[,nn])
} 
rand.pca=prcomp(rand.data,scale=TRUE)
cor.rand[jj]=cor(rand.pca$x[,1],Mstatepca$x[,1])
}
plot(abs(cor.rand),xx,xlab="Absolute R squared value",ylab="Number of probes selected for PCA")



################## Testing correlations with randomly selected genes corr with randomly selected genes #############
cor.rand=c()
xx=rep((1:10)^2,10)
for(jj in 1:100){
  rand.gene.sel=gene.name[floor(runif(15,min=1,max=length(gene.name)))]
  rand.gene.sel2=gene.name[floor(runif(15,min=1,max=length(gene.name)))]
  rand.data=matrix(nrow=291,ncol=0)
  rand.data2=matrix(nrow=291,ncol=0)
  for (ii in 1:length(rand.gene.sel)){
    nn=grep(paste0("^",rand.gene.sel[ii],"$"),gene.name)
    rand.data=cbind(rand.data,qnormed[,nn])
    nn2=grep(paste0("^",rand.gene.sel2[ii],"$"),gene.name)
    rand.data2=cbind(rand.data2,qnormed[,nn2])
  } 
  rand.pca=prcomp(rand.data,scale=TRUE)
  rand.pca2=prcomp(rand.data2,scale=TRUE)
  cor.rand[jj]=cor(rand.pca$x[,1],rand.pca2$x[,1])
}
plot(abs(cor.rand),xx,xlab="Absolute R squared value",ylab="Number of probes selected for PCA")


#### Created normally distributed data and tested PCA of truly random data ######
tru.cor.rand=c()
for(jj in 1:100){
tru.rand.data=matrix(nrow=291,ncol=0)
tru.rand.data2=matrix(nrow=291,ncol=0)
for(ii in 1:xx[jj]){
mm=rnorm(291,mean=8.175,sd=1.62)
tru.rand.data=cbind(tru.rand.data,mm)
tt=rnorm(291,mean=8.175,sd=1.62)
tru.rand.data2=cbind(tru.rand.data2,tt)
}
tru.rand.pca=prcomp(tru.rand.data,scale=TRUE)
tru.rand.pca2=prcomp(tru.rand.data2,scale=TRUE)
tru.cor.rand[jj]=cor(tru.rand.pca$x[,1],tru.rand.pca2$x[,1])
}
plot(abs(tru.cor.rand), xx,xlab="Absolute R squared value",ylab="Number of random variables used for PCA")





############## Trying to get rid of meaningless variation##########
library(ppcor) 

pcor(cbind(rand.pca$x[,1],rand.pca2$x[,1],house.pca$x[,1]))

pcor(cbind(gene.fam.pca$x[,1],Mstatepca$x[,1],house.pca$x[,1]))


GibbsTable=cbind(gene.names.oi,probes.repre,genes.oi.loadings,individ.cor.m0,individ.cor.m1,individ.cor.m2)
write.table(GibbsTable,file="Genes of Interest Gibbs Data Correlations 5-2-16CS",quote=FALSE,sep="\t",col.names = c("GibbsGeneNames","GeneFamily","GenePCALoadings","CorrM0PCA","CorrM1PCA","CorrM2PCA"))

read.table("Genes of Interest Gibbs Data Correlations 1-19-17CS")




plot(gene.fam.pca$x[,1],house.pca$x[,1])


summary(apply(qnormed,2,mean))

meandif=matrix(ncol=ncol(rand.data),nrow=nrow(rand.data))
xbar=apply(rand.data,2,mean)
for(ii in 1:ncol(rand.data)){
  meandif[,ii]= rand.data[,ii] - xbar[ii]
}
cov(meandif)

rr=summary(lm(rand.data~ house.pca$x[,1]))

pvals=c()
for(ii in 1:length(rr[[1]]))
  rr[[1]]$coefficients[2,4]

