# GO Ontology downloader for M1 gene List

name="GOontologyM1"
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

############ Directory Commands ##################
# Type this before creating any output files
setwd(outDirDay)
# Type this before creating any graphs
setwd(outDirGraphsDay)
# Type this before referencing something in the main directory
setwd(mainDir)
# Type this before referencing original data
setwd(input)


library("GO.db")
ls("package:GO.db")
GOterms.info <- as.list(GOTERM)

library("org.Hs.eg.db")
library("TxDb.Hsapiens.BioMart.igis")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
xx=org.Hs.egSYMBOL
EZgenes <- mappedkeys(xx)
EZgenes=as.data.frame(xx[EZgenes])

gene.names=EZgenes[match(M1$Gene.Symbol,EZgenes$symbol),] #This tells you in order of the queried entries,what the coordinates are of the second group of entries only  the first match is given 
gene.names=gene.names[!is.na(gene.names[,2]),]

xx=org.Hs.egGO
GOterms <- mappedkeys(xx)
GOterms=as.data.frame(xx[GOterms])

M1GOterms=data.frame()
for(ii in 1:nrow(gene.names)){
 xx=grep(paste0('^',gene.names$gene_id[ii],'$'),GOterms$gene_id)
 xx=cbind(rep(gene.names$symbol[ii],length(xx)),GOterms[xx,])
 M1GOterms=rbind(M1GOterms,xx)
}

M1GOterms=M1GOterms[which(M1GOterms$Ontology=="BP"),]
M1GOterms.info=GOterms.info[match(M1GOterms$go_id,names(GOterms.info))]
write.csv(M1GOterms,file="GO Ontology for M1 Genes.csv")
write.csv(M1GOterms.info)
M1GOfunctions=c()
for(ii in 1:length(M1GOterms.info)){
M1GOfunctions=c(M1GOfunctions,M1GOterms.info[[ii]]@Term)
}
xx=table(M1GOfunctions)
par(mar=c(15,3,3,3),cex.lab=.4,cex.axis=.5)
barplot(sort(xx[which(xx>5)],decreasing=TRUE),las=2)
dev.off()
par(mar=c(15,3,3,3),cex.lab=.4,cex.axis=.2)
barplot(sort(table(M1GOterms$`rep(gene.names$symbol[ii], length(xx))`),decreasing=TRUE),las=2)


### How to choose which evidence is the strongest
sort(table(M1GOterms$Evidence),decreasing=TRUE)
# IEA TAS IDA IMP ISS IBA NAS  IC IEP IGI  ND EXP 
# 513 379 287 143 122  67  39  29  12  11   6   1 
rank.evidence=M1GOterms$Evidence
rank.evidence=mapvalues(rank.evidence,from = c(unique(M1GOterms$Evidence)), to=c(12,4,9,1,10,6,5,11,8,2,3,7))
rank.evidence=as.numeric(rank.evidence)
gene.GO=as.character(M1GOterms$`rep(gene.names$symbol[ii], length(xx))`)
bestGO=mapvalues(rank.evidence,from = c(12,4,9,1,10,6,5,11,8,2,3,7), to=unique(M1GOterms$Evidence))

## Picks the GO function with the stronger evidence
bestGO=c()
yy=c()
for(ii in 1:length(unique(gene.GO))){
  xx=grep(paste0('^',unique(gene.GO)[ii],'$'),gene.GO)
  xx=M1GOterms$go_id[xx[which(rank.evidence[xx]==min(rank.evidence[xx]))]]
  bestGO=c(bestGO,xx)
  yy=as.vector(c(yy,rep(unique(gene.GO)[ii],length(xx))))
}
names(bestGO)=yy

##### Finds the description of the GO functions of the geneset
bestGO.info=GOterms.info[match(bestGO,names(GOterms.info))]
bestGOfunctions=c()
names(bestGO.info)=yy
#yy=c()
for(ii in 1:length(bestGO.info)){
  bestGOfunctions=c(bestGOfunctions,bestGO.info[[ii]]@Term)
  #yy=c(yy,bestGO.info[[ii]]@GOID)
}
#names(bestGOfunctions)=yy
names(bestGOfunctions)=yy

par(mar=c(12,3,3,3),cex.lab=.4,cex.axis=.5)
barplot(sort(table(bestGOfunctions),decreasing=TRUE)[1:25],las=2)
dev.off()


### Takes the most common GO functions for each gene
GOcommon=mapvalues(bestGOfunctions,from = sort(unique(bestGOfunctions)), to=order(table(bestGOfunctions),decreasing=TRUE))
bestGO2=c()
yy=c()
for(ii in 1:length(unique(names(bestGOfunctions)))){
  xx=grep(paste0('^',unique(names(bestGOfunctions))[ii],'$'),names(bestGOfunctions))
  if(length(xx)==1){
    xx=bestGOfunctions[xx]
    bestGO2=c(bestGO2,xx)
  }else{
    xx=bestGOfunctions[xx[which(GOcommon[xx]==min(GOcommon[xx]))]]
    bestGO2=c(bestGO2,xx)
  }
  yy=c(yy,(unique(names(bestGO.info)))[ii])
}
names(bestGO2)=yy
sort(table(bestGO2))

par(mar=c(12,3,3,3),cex.lab=.4,cex.axis=.7)
barplot(sort(table(bestGO2),decreasing=TRUE)[1:15],las=2)
