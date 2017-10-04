setwd("/Users/cansav091/Desktop/Neurogenomics Data/Martinez M0M1M2 Lists/")
load("Gibbs.cortex.qnormed.Rdata")

output="/Users/cansav091/Desktop/Neurogenomics Data/Martinez M0M1M2 Lists/GABAGluGibbs"
input="/Users/cansav091/Desktop/Neurogenomics Data/Martinez M0M1M2 Lists/Original Data"

################# Download Annotation for Gene(s) of Interest ############### 
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
  assign(paste0(gene.family[jj],".gene.family.info"),gene.family.info,envir=.GlobalEnv)
}

#######################################################################################
################## PCA for gene family of interest ####################################
for(jj in 1:length(gene.family)){
  gene.family.info=eval(parse(text=paste0(gene.family[jj],".gene.family.info")))
  gene.names.oi=as.character(gene.family.info$`Approved Symbol`)
  
  gene.family.data=matrix(nrow=nrow(qnormed),ncol=0)
  gene.family.names=c()
  gene.fam.probes=c()
  for (ii in 1:length(gene.names.oi)){
    nn=grep(paste0("^",gene.names.oi[ii],"$"),toupper(gene.name))
    gene.family.data=cbind(gene.family.data,qnormed[,nn])
    gene.family.names=c(gene.family.names,gene.name[nn])
    gene.fam.probes=c(gene.fam.probes,rep(gene.names.oi[ii],length(nn)))
  } 
  assign(paste0(gene.family[jj],".pca"),prcomp(gene.family.data,scale=TRUE),envir = .GlobalEnv)
  assign(paste0(gene.family[jj],".gibbs"),gene.family.data,envir = .GlobalEnv)
  assign(paste0(gene.family[jj],".gene.names"),gene.family.names,envir = .GlobalEnv)
  assign(paste0(gene.family[jj],".subgroups"),unique(gene.family.info$`Gene Family Tag`))
}

setwd(input)
microglia.activ=read.csv("M0,M1,M2 MANOVA Gene List FDR= 0.05 FC cutoff= 1.2 Microglia M1 Markers.csv")
microglia.activ=microglia.activ[which(microglia.activ$FDR.adj.p.val<.05),]

microglia.activ.gibbs=matrix(nrow=nrow(qnormed),ncol=0)
microglia.activ.genes=c()
for (ii in 1:length(microglia.activ$X)){
  nn=grep(paste0("^",microglia.activ$X[ii],"$"),toupper(gene.name))
  microglia.activ.gibbs=cbind(microglia.activ.gibbs,qnormed[,nn])
  microglia.activ.genes=c(microglia.activ.genes,gene.name[nn])
} 

microglia.activ.pca=prcomp(microglia.activ.gibbs,scale=TRUE)
#######################################################################################

plot(Glutamate.pca$x[,1]~GABA.pca$x[,1])
xx=lm(Glutamate.pca$x[,1]~GABA.pca$x[,1])
abline(xx)
abline(a=1.362048e-15 ,b=  1.492178 )

##################### What genes are these PCs correlated with? ##########################
setwd(output)
if(ncol(qnormed)<nrow(qnormed)){
qnormed=t(qnormed)
}
GABA.gene.cor=matrix(ncol=ncol(qnormed),nrow=2)
dimnames(GABA.gene.cor)[[2]]=gene.name
for(ii in 1:ncol(qnormed)){
  xx=cor.test(qnormed[,ii],GABA.pca$x[,1])
  GABA.gene.cor[1,ii]=xx$estimate
  GABA.gene.cor[2,ii]=xx$p.value
}
GABA.gene.cor[2,]=p.adjust(GABA.gene.cor[2,],method="hochberg")

#GABA.gene.cor=GABA.gene.cor[,which(GABA.gene.cor[2,]<.05)]
write.csv(t(GABA.gene.cor[,order(abs(GABA.gene.cor[2,]))][,1:100]),file="GABA PC with Other Gene Correlations Gibbs Data.csv")


Glutamate.gene.cor=matrix(ncol=ncol(qnormed),nrow=2)
dimnames(Glutamate.gene.cor)[[2]]=gene.name
for(ii in 1:ncol(qnormed)){
  xx=cor.test(qnormed[,ii],Glutamate.pca$x[,1])
  Glutamate.gene.cor[1,ii]=xx$estimate
  Glutamate.gene.cor[2,ii]=xx$p.value
}
Glutamate.gene.cor[2,]=p.adjust(Glutamate.gene.cor[2,],method="hochberg")

#Glutamate.gene.cor=Glutamate.gene.cor[,which(Glutamate.gene.cor[2,]<.05)]
write.csv(t(Glutamate.gene.cor[,order(abs(Glutamate.gene.cor[2,]))][,1:100]),file="Glutamate PC with Other Gene Correlations Gibbs Data.csv")

microglia.activ.gene.cor=matrix(ncol=ncol(qnormed),nrow=2)
dimnames(microglia.activ.gene.cor)[[2]]=gene.name
for(ii in 1:ncol(qnormed)){
  xx=cor.test(qnormed[,ii],microglia.activ.pca$x[,1])
  microglia.activ.gene.cor[1,ii]=xx$estimate
  microglia.activ.gene.cor[2,ii]=xx$p.value
}
microglia.activ.gene.cor[2,]=p.adjust(microglia.activ.gene.cor[2,],method="hochberg")

setwd(input)
read.csv("HousekeepingGenes.txt",sep="\t")
house.genes=read.table("HousekeepingGenes.txt",sep="\t")
house.genes=as.character(house.genes[,1])
house.genes=substr(house.genes,1,nchar(house.genes)-1)

house.data=matrix(nrow=nrow(qnormed),ncol=0)
house.gene.gibbs=c()
for (ii in 1:length(house.genes)){
  nn=grep(paste0("^",house.genes[ii],"$"),gene.name)
  house.data=cbind(house.data,qnormed[,nn])
  house.gene.gibbs=c(house.gene.gibbs,house.genes[nn])
}

house.pca=prcomp(house.data[,],scale=TRUE)
house.pca=prcomp(house.data[,runif(20,min=0,max=ncol(house.data))],scale=TRUE)


house.gene.cor=matrix(ncol=ncol(qnormed),nrow=2)
dimnames(house.gene.cor)[[2]]=gene.name
for(ii in 1:ncol(qnormed)){
  xx=cor.test(qnormed[,ii],house.pca$x[,1])
  house.gene.cor[1,ii]=xx$estimate
  house.gene.cor[2,ii]=xx$p.value
}
house.gene.cor[2,]=p.adjust(house.gene.cor[2,],method="hochberg")


###################### Basic stats on correlations of each PC with the transcriptome ##########
hist(abs(Glutamate.gene.cor[1,]))
mean(abs(Glutamate.gene.cor[1,])) ### R=.37

hist(abs(GABA.gene.cor[1,]))
mean(abs(GABA.gene.cor[1,])) ### R=.38

hist(abs(microglia.activ.gene.cor[1,]))
mean(abs(microglia.activ.gene.cor[1,]))### R=.36

hist(abs(house.gene.cor[1,])) 
mean(abs(house.gene.cor[1,])) ### R=.34


plot(GABA.gene.cor[1,],Glutamate.gene.cor[1,])
cor(GABA.gene.cor[1,],Glutamate.gene.cor[1,])
plot(-log10(GABA.gene.cor[2,]),-log10(Glutamate.gene.cor[2,]))

plot(GABA.gene.cor[1,],microglia.activ.gene.cor[1,])
cor(GABA.gene.cor[1,],microglia.activ.gene.cor[1,])
plot(-log10(GABA.gene.cor[2,]),-log10(microglia.activ.gene.cor[2,]))

plot(Glutamate.gene.cor[1,],microglia.activ.gene.cor[1,])
cor(Glutamate.gene.cor[1,],microglia.activ.gene.cor[1,])
plot(-log10(Glutamate.gene.cor[2,]),-log10(microglia.activ.gene.cor[2,]))

plot(-log10(Glutamate.gene.cor[2,]),-log10(house.gene.cor[2,]))
plot(-log10(GABA.gene.cor[2,]),-log10(house.gene.cor[2,]))
plot(-log10(microglia.activ.gene.cor[2,]),-log10(house.gene.cor[2,]))

plot(Glutamate.gene.cor[1,],house.gene.cor[1,])
plot(GABA.gene.cor[1,],house.gene.cor[1,])
plot(microglia.activ.gene.cor[1,],house.gene.cor[1,])









