# Investigating current M state gene list with Cell specific expression data from Zhang

name="CompareZhangwMState"
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

Mstate.data=read.csv("M1 Genes Hotel p<.05 FC>1.3.csv")
zhang=read.csv("barreslab_rnaseq.csv")


xx=match(Mstate.data$X,toupper(zhang$Gene.symbol))
Mstate.data=Mstate.data[!is.na(xx),] # 127 of 144 are in the Zhang list 
zhang.mstate=zhang[xx[!is.na(xx)],]


##### Don't want to have genes in my Mstate list that are highly expressed in other cell types 
plot(apply(zhang.mstate[,c(3:6,8)],1,max),Mstate.data$M1)
textxy(apply(zhang.mstate[,c(3:6,8)],1,max),Mstate.data$M1,Mstate.data$X)
which(apply(zhang.mstate[,c(3:6,8)],1,max)>2)

### How does the PC's look if I remove these genes? 

Mstate.data=Mstate.data[which(apply(zhang.mstate[,c(3:6,8)],1,max)<2),]

Mstategib.edit=matrix(nrow=291,ncol=0)#Need an empty vector to put the data I want in
Mstategib.gene.name.edit=c()#Need an empty vector for the gene names that correspond to the matrix
for (ii in 1:length(Mstate.data$X)){
  nn=which(gene.name==Mstate.data$X[ii])#gives indices for gene.name that match the names from m1gene names
  Mstategib.edit=cbind(Mstategib.edit,qnormed[,nn])#Add the correspondin columns to the growing list of the dataset
  Mstategib.gene.name.edit=c(Mstategib.gene.name.edit,gene.name[nn])#Create a vector with the gene names for the new dataset
}

dimnames(Mstategib.edit)[[2]]=Mstategib.gene.name.edit

#Previous plot with no removals of genes based on Zhang data
jpeg("Plot of PCA with NO removal of other cell expressing genes.jpeg")
plot(Mstatepca$rotation)
textxy(Mstatepca$rotation[,1],Mstatepca$rotation[,2],Mstategib.gene.name)
dev.off()

setwd(outDirGraphs)
Mstatepca.edit=prcomp(Mstategib)
jpeg("Plot of PCA after removal of other cell expressing genes.jpeg")
plot(Mstatepca.edit$rotation)
textxy(Mstatepca.edit$rotation[,1],Mstatepca.edit$rotation[,2],Mstategib.gene.name)
dev.off()

