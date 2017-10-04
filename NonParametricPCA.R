# Non parametric PCA Example

name="NonParametricPCA"
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


# Test PCA Example 
# Randomly chose 100 genes 
rand.gene.sel=gene.name[floor(runif(100,min=1,max=length(gene.name)))]
rand.data=matrix(nrow=291,ncol=0)
for (ii in 1:length(rand.gene.sel)){
  nn=grep(paste0("^",rand.gene.sel[ii],"$"),gene.name)
  rand.data=cbind(rand.data,qnormed[,nn])
  pca.tmp=prcomp(rand.data,scale=TRUE)
} 
### Looking at each step of PCA 
data=synapse.gene.data
rank.data=data
for(ii in 1:ncol(data)){
  rank.data[,ii]=order(data[,ii])
}

# Found the variable means 
means=apply(rank.data,2,mean)
# Found the mean difs for each sample 
mean.difs=apply(rank.data,1,function(x) x-means)
hist(cor(mean.difs))
hist(apply(cor(mean.difs),1,mean))
# Found the covariance of each 
covmat=cov(rank.data)
# Found the Eigenvalues for the covariance matrix
eigs=eigen(covmat)
hist(eigs$vectors)
# Found the Z scores for each sample for each variable
Zscor=apply(rank.data,1,function(x) (x-apply(rank.data,2,mean))/sqrt(apply(rank.data,2,sd)))
hist(cor(Zscor))
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

plot(pc1scor,Mstatepc1scor)
cor(pc1scor,Mstatepc1scor)

Mstatepc1scor=pc1scor






