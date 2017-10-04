
### Try out the synaptic markers found in the Counts et al, 2014 Alzheimers paper. 

counts.syn.markers=read.table("CountsEtAlSynapseMarkers.txt")
counts.syn.markers=toupper(as.character(counts.syn.markers$V1))

syn.marker.data=matrix(nrow=291,ncol=0)
syn.marker.match=c()
syn.marker.probes=c()
for (ii in 1:length(counts.syn.markers)){
  nn=grep(paste0("^",counts.syn.markers[ii],"$"),gene.name)
  syn.marker.data=cbind(syn.marker.data,qnormed[,nn])
  syn.marker.match=c(syn.marker.match,gene.name[nn])
  syn.marker.probes=c(syn.marker.probes,rep(counts.syn.markers[ii],length(nn)))
} 



pca.syn=prcomp( syn.marker.data)
plot(pca.syn$rotation)

###### Check how this PCA correlates with M1 PCA 
plot(pca.syn$x[,1],Mstatepca.edit$x[,1])
plot(pca.syn$x[,1],prcomp(rand.data)$x[,1])


###### Check if this also has inflated correlations 
xx=rep((1:10)^2,100)
rand.pca.x=matrix(nrow=291,ncol=1000)
rand.pca.load=list()
for(jj in 460:1000){
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


plot(xx, cor())

