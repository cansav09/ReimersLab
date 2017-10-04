#Investigating the PCA Correlation problem. Using Raw data for Gibbs data 

source("https://bioconductor.org/biocLite.R")
biocLite("illuminaio")
####### Testing correlations with randomly selected housekeping genes with m1 PCA #############
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
#cor.house[jj]=cor(house.pca$x[,1],mstatepca$x[,1])
#}
house.pca=prcomp(house.data,scale=TRUE)
setwd(outDirDay)
jpeg(paste(mstate,"and Housekeeping Genes","PC scatterplot.jpeg"))
cor.house=cor(house.pca$x[,1],house.data)
plot(house.pca$x[,1],mstatepca$x[,1],ylab=paste0(mstate,"PC"),xlab= paste("Housekeeping PC"),sub=paste("R=",round(cor(house.pca$x[,1],mstatepca$x[,1]),3)),cex.lab=1.5)
dev.off()

################## Testing correlations with randomly selected genes and m1 PCA#############
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
    cor.rand[jj]=cor(rand.pca$x[,1],mstatepca$x[,1])
}
plot(abs(cor.rand),xx,xlab="Absolute R squared value",ylab="Number of probes selected for PCA")


############## Trying to get rid of meaningless variation##########
library(ppcor)

pcor.test(gene.fam.pca$x[,1],house.pca2$x[,1],house.pca$x[,1])
pcor(cbind(gene.fam.pca$x[,1],mstatepca$x[,1],house.pca$x[,1]))

GibbsTable=cbind(gene.names.oi,probes.repre,genes.oi.loadings,individ.cor.m0,individ.cor.m1,individ.cor.m2)
write.table(GibbsTable,file="Genes of Interest Gibbs Data Correlations 5-2-16CS",quote=FALSE,sep="\t",col.names = c("GibbsGeneNames","GeneFamily","GenePCALoadings","Corrm0PCA","Corrm1PCA","Corrm2PCA"))

read.table("Genes of Interest Gibbs Data Correlations 1-19-17CS")

plot(gene.fam.pca$x[,1],house.pca$x[,1])

summary(apply(qnormed,2,mean))

meandif=matrix(ncol=ncol(rand.data),nrow=nrow(rand.data))
xbar=apply(rand.data,2,mean)
for(ii in 1:ncol(rand.data)){
    meandif[,ii]= rand.data[,ii] - xbar[ii]
}


################## Testing correlations with randomly selected genes corr with randomly selected genes #############
cor.rand=c()
xx=rep((1:10)^2,10)
for(jj in 1:100){
  rand.gene.sel=gene.name[floor(runif(xx[jj],min=1,max=length(gene.name)))]
  rand.gene.sel2=gene.name[floor(runif(xx[jj],min=1,max=length(gene.name)))]
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


#### Turn normal data into rank data
rank.data=qnormed
for(ii in 1:ncol(qnormed)){
  rank.data[,ii]=order(qnormed[,ii])
}

#### Use the rank data to do the same exp with using random sets of probes 
cor.rand=c()
xx=rep((1:10)^2,10)
for(jj in 1:100){
  rand.gene.sel=gene.name[floor(runif(xx[jj],min=1,max=length(gene.name)))]
  rand.gene.sel2=gene.name[floor(runif(xx[jj],min=1,max=length(gene.name)))]
  rand.data=matrix(nrow=291,ncol=0)
  rand.data2=matrix(nrow=291,ncol=0)
  for (ii in 1:length(rand.gene.sel)){
    nn=grep(paste0("^",rand.gene.sel[ii],"$"),gene.name)
    rand.data=cbind(rand.data,rank.data[,nn])
    nn2=grep(paste0("^",rand.gene.sel2[ii],"$"),gene.name)
    rand.data2=cbind(rand.data2,rank.data[,nn2])
  } 
  rand.pca=prcomp(rand.data,scale=TRUE)
  rand.pca2=prcomp(rand.data2,scale=TRUE)
  cor.rand[jj]=cor(rand.pca$x[,1],rand.pca2$x[,1])
}

jpeg("Nonparametric PCA rand cor.jpeg")
plot(xx,cor.rand)
dev.off()
boxplot(rand.pca$x[,1],as.factor(region),names=c("FCTX","TCTX"))



xx=rep((1:10)^2,10)
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

pcor.test(gene.fam.pca$x[,1],house.pca2$x[,1],house.pca$x[,1])
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

## Second part of this question: is this magical PC consistent for a given sample across trials? 
xx=rep((1:10)^2,100)
rand.pca.x=matrix(nrow=291,ncol=1000)
rand.pca.load=list()
for(jj in 1:1000){
  rand.gene.sel=gene.name[floor(runif(xx[jj],min=1,max=length(gene.name)))]
  rand.data=matrix(nrow=291,ncol=0)
    for (ii in 1:length(rand.gene.sel)){
    nn=grep(paste0("^",rand.gene.sel[ii],"$"),gene.name)
    rand.data=cbind(rand.data,qnormed[,nn])
    pca.tmp=prcomp(rand.data,scale=TRUE)
    } 
  rand.pca.x[,jj]=pca.tmp$x[,1]
  rand.pca.load[[jj]]=pca.tmp$rotation
  }


boxplot(rand.pca$x[,1],as.factor(region),names=c("FCTX","TCTX"))
plot(rand.pca$x[,1][which(region=="TCTX")],rand.pca2$x[,1][which(region=="TCTX")])
textxy(rand.pca2$x[,1][which(region=="FCTX")],dimnames(qnormed)[[1]][which(region=="FCTX")])## Delete this if you don't want the gene names on the plot
dev.off()

### Checking if the housekeeping genes as a partial correlation control
pcor.test(gene.fam.pca$x[,1],Mstatepca$x[,1], house.pca$x[,1])
cor(cbind(gene.fam.pca$x[,1],rand.pca.x))[1:5]

rpca=matrix(nrow=291,ncol=1000)
for(ii in 1:1000){
  rpca[,ii]=prcomp(rnorm(291))$x[,1]
}

setwd(outDirGraphs)
jpeg("Mstate Correlations with GABR Gibbs Data.jpeg")
barplot(sort(t(cor(Mstatepca.edit$x[,1],gene.family.data))),names.arg=gene.oi.match[order(t(cor(Mstatepca.edit$x[,1],gene.family.data)))],las=2,cex.names=.8)
dev.off()


# Test PCA Example 
# Randomly chose 100 genes 
rand.gene.sel=gene.name[floor(runif(100,min=1,max=length(gene.name)))]
rand.data=matrix(nrow=291,ncol=0)
for (ii in 1:length(rand.gene.sel)){
  nn=grep(paste0("^",rand.gene.sel[ii],"$"),gene.name)
  rand.data=cbind(rand.data,qnormed[,nn])
  pca.tmp=prcomp(rand.data,scale=TRUE)
} 




setwd(outDirGraphs)
### Looking at each step of PCA 
rand.data
jpeg("Hist rand probe data correlations.jpeg")
hist(cor(rand.data)[which(cor(rand.data)!=1)])
dev.off()
# Found the variable means 
means=apply(rand.data,2,mean)
# Found the mean difs for each sample 
mean.difs=apply(rand.data,1,function(x) x-means)
jpeg("Hist mean.difs correlations.jpeg")
hist(cor(mean.difs))
dev.off()
hist(apply(cor(mean.difs),1,mean))
# Found the covariance of each 
covmat=cov(rand.data)
# Found the Eigenvalues for the covariance matrix
eigs=eigen(covmat)
hist(eigs$vectors)
# Found the Z scores for each sample for each variable
Zscor=apply(rand.data,1,function(x) (x-apply(rand.data,2,mean))/sqrt(apply(rand.data,2,sd)))
jpeg("Hist Zscor correlations.jpeg")
hist(cor(Zscor))
dev.off()
# Sort Eigenvalues
eigs.sort=eigs$vectors[,order(eigs$values,decreasing=TRUE)]
eigs$values/sum(eigs$values) #percentage variance due to each. 
# Calculate Loadings for each variables
loads=eigs$vectors*sqrt(eigs$values)
# Calculate PCA scores for the first PC
pc1scor=c(1:291)
for(ii in 1:291){
pc1scor[ii]=sum(mean.difs[,ii]*eigs.sort[,1])
}
pc2scor=c(1:291)
for(ii in 1:291){
  pc2scor[ii]=sum(mean.difs[,ii]*eigs.sort[,2])
}

cor(prcomp(rand.data)$x[,1],prcomp(rand.data)$x[,2])
head(prcomp(rand.data)$x[,1])




