###goal: Look at genes of interest and their relationship to M state gene sets 
### The output of this script is 
## 1) Three graphs that show the correlation of the combined composite score for all the probes that match the gene of interest
## 2) A table that shows 
##    a) The probes that were matched to the gene of interest that was inputted
##    b) Their respective loadings within their gene family
##    c) The correlation of each individual probe with each m state PCA composite variable

## Pre step:  Run the m0 m1 m2 gibbs code first to have the m state PCAs
name="GeneOfInterestGibbsQueryMstate"
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
#### Put a list of gene symbols that you want to lookup in the parantheses below using regexpr

genes.of.interest=c('^GRI[ANK]*','GABR[A-Z]','GLRA*')

genes.of.interest.gibbs=matrix(nrow=291,ncol=0)## This is where the gibbs data for a probe that matches the gene of interest is placed
gene.names.oi=c()## This variable shows the gene symbol as it is in the gibbs dataset
genes.oi.loadings=c()## This variable will show the loading for that particular probe for a given gene in the gibbs dataset
probes.repre=c()## This variable shows the the gene symbol entry (from "genes.of.interest) that matched the individual probe
individ.cor.m0=c()### These variables will show the correlation of each individual probe to the PCA scores for the respective "M states"
individ.cor.m1=c()
individ.cor.m2=c()
for (i in 1:length(genes.of.interest)){
  nn=grep(genes.of.interest[],gene.name)
  genes.of.interest.gibbs=cbind(genes.of.interest.gibbs,qnormed[,nn])
  gene.names.oi=c(gene.names.oi,gene.name[nn])
  nnp=prcomp(qnormed[,nn])
  probes.repre=c(probes.repre,rep(gene.of.interest[i],length(nn)))
  genes.oi.loadings=c(genes.oi.loadings,nnp$rotation[,1])
  individ.cor.m0=c(individ.cor.m0,cor(genes.of.interest.gibbs,m0p$x[,1]))
  individ.cor.m1=c(individ.cor.m1,cor(genes.of.interest.gibbs,m1p$x[,1]))
  individ.cor.m2=c(individ.cor.m2,cor(genes.of.interest.gibbs,m2p$x[,1]))
###Make graph and correlation for M0
jpeg(paste0(genes.of.interest[i],"M0corr.jpeg",date()))
corM=cor(nnp$x[,1],m0p$x[,1])
plot(nnp$x[,1],m0p$x[,1] ,ylab="M0 composite score",xlab= paste(gene.of.interest[i],"gene composite score"),sub=paste("R=",corM),cex.lab=1.5)
dev.off()
###Make graph and correlation for M1
jpeg(paste0(genes.of.interest[i],"M1corr.jpeg",date()))
corM=cor(nnp$x[,1],m1p$x[,1])
plot(nnp$x[,1],m1p$x[,1] ,ylab="M1 composite score",xlab= paste(gene.of.interest[i],"gene composite score"),sub=paste("R=",corM),cex.lab=1.5)
dev.off()
###Make graph and correlation for M2
jpeg(paste0(genes.of.interest[i],"M2corr.jpeg",date()))
corM=cor(nnp$x[,1],m2p$x[,1])
plot(nnp$x[,1],m2p$x[,1] ,ylab="M2 composite score",xlab= paste(gene.of.interest[i],"gene composite score"),sub=paste("R=",corM),cex.lab=1.5)
dev.off()
}


### Save all the data into a table
write.table(cbind(gene.names.oi,probes.repre,genes.oi.loadings,individ.cor.m0,individ.cor.m1,individ.cor.m2),file=paste("Genes of Interest Gibbs Data Correlations",date()),quote=FALSE,sep="\t",col.names = c("GibbsGeneNames","GeneFamily","GenePCALoadings","CorrM0PCA","CorrM1PCA","CorrM2PCA"))
