allen.cell.type=as.matrix(read.csv("cell_classification.csv"))
allen.rpkms=as.matrix(read.csv("genes_rpkm.csv",row.names=1))
allen.meta=as.data.frame(read.csv("cell_metadata.csv"))

## checking that the ids are in the same order 
match(allen.meta[,1],dimnames(allen.rpkms)[[2]][1:1679])==1:1679


summary(allen.meta$sub_class)
###                      Astrocyte           Chodl     Endothelial            L2/3      L2/3_Syt10 
#43              50               2              21              94              12 
#L4              L5       L5_Chrna6             L5a            L5a1            L5a2 
#285               6               3              43              10              59 
#L5b            L5b1            L5b2             L6a            L6a1            L6a2 
#35              13              12              12              80              27 
#L6b       Microglia            Ndnf Oligodendrocyte             OPC           Pvalb 
#39              22              78              42              24             270 
#Sst       Sst_Cbln4       Sst_Chodl      Sst_Clstn2        Sst_Nmbr       Sst_Nr2f2 
#7              38              40              14              10              69 
#Th         Unknown             Vip 
#14              10             195 

summary(allen.meta$major_class)

##Endothelial  Excitatory        Glia  Inhibitory     Unknown 
##21         756         138         761           3 

neuron.avg.allen=apply(allen.rpkms[,allen.meta[,11] %in% c("Excitatory","Inhibitory")],1,mean)
microglia.avg.allen=apply(allen.rpkms[,allen.meta[,12] %in% c("Microglia")],1,mean)
astrocyte.avg.allen=apply(allen.rpkms[,allen.meta[,12] %in% c("Astrocyte")],1,mean)
oligo.avg.allen=apply(allen.rpkms[,allen.meta[,12] %in% c("Oligodendrocyte")],1,mean)

fold.change.microglia.allen=microglia.avg.allen/(neuron.avg.allen+astrocyte.avg.allen+oligo.avg.allen)

allen.cells=cbind(allen.rpkms[,allen.meta[,12] %in% c("Microglia")],allen.rpkms[,allen.meta[,11] %in% c("Excitatory","Inhibitory")],allen.rpkms[,allen.meta[,12] %in% c("Astrocyte")],allen.rpkms[,allen.meta[,12] %in% c("Oligodendrocyte")])

##Remove low abundance in Microglia genes 
allen.cells=allen.cells[which(microglia.avg.allen>1),]
microglia.avg.allen=microglia.avg.allen[which(microglia.avg.allen>1)]
fold.change.microglia.allen=fold.change.microglia.allen[which(microglia.avg.allen>1)]

##Doing t tests with the original dataset
library(genefilter) # we load this library for fastT to do many t-tests in parallel
allen.t.scores <- fastT(allen.cells, 1:22, 23:1631)$z # generate many t-scores in parallel  
allen.pvalues=2 * pt( -abs(allen.t.scores), df = length(allen.cells[1,]-2))
sum(!is.na(which(allen.pvalues<.05))) #without correction 6934 probes are significant
allen.pvalues=p.adjust(allen.pvalues,method="fdr")
sum(!is.na(which(allen.pvalues<.05))) #with correction 5524 probes are significant 

## Taking out the biggest 500 fold changes 
allen.rpkms.microglia=allen.cells[order(fold.change.microglia.allen,decreasing=TRUE),][1:500,]
genes.fc.allen=cbind(dimnames(allen.cells)[[1]][1:500],fold.change.microglia.allen[order(fold.change.microglia.allen,decreasing=TRUE)][1:500])



##Compare to darmanis and zhang
 
darmanis=cbind(microglia.genes,microglia.signal,nonmicroglia.signal,microglia.fc)
zhang=cbind(microglia.genes.zhang,microglia.signal.zhang,nonmicroglia.signal.zhang,microglia.fc.zhang)
allen=cbind(toupper(dimnames(allen.cells)[[1]][order(fold.change.microglia.allen,decreasing=TRUE)],microglia.avg.allen[order(fold.change.microglia.allen,decreasing=TRUE)],fold.change.microglia.allen[order(fold.change.microglia.allen,decreasing=TRUE)])
allen[which(allen[,3]=="Inf")]=max(allen[,3])

darmanis=darmanis[which(darmanis[,2]>1),] 
zhang=zhang[which(zhang[,2]>1),] 

zhang=zhang[1:500,] 
zhang[,1]=toupper(zhang[,1])
### Let's get rid of genes in the zhang list that don't appear to have any human equivalent
library(org.Hs.eg.db)
EGIDs <- keys(org.Hs.eg.db)
genesymbol.key=unlist(mget(EGIDs, org.Hs.egSYMBOL))# N=56340 genes in this homo sapiens database

###go back and get the full list from Zhang 
#Get rid of zhang genes without human matches 604 genes not matching (by using symbols anyway)
zhang=zhang[which(!is.na(match(zhang[,1],genesymbol.key))),] 

## Now lets make a top 500 genes and cross reference them to darmanis
overlap.allen=allen[!is.na(match(allen[,1],darmanis[,1])),]
overlap.darmanis=darmanis[!is.na(match(darmanis[,1],allen[,1])),]

overlap.allen=overlap.allen[order(overlap.allen[,1]),]
overlap.darmanis=overlap.darmanis[order(overlap.darmanis[,1]),]

match(overlap.darmanis[,1],overlap.allen[,1])

plot(log2(as.numeric(overlap.darmanis[,4])),log2(as.numeric(overlap.allen[,3])),xlab="Darmanis Fold Change",ylab="Allen Fold Change")

cor(log2(as.numeric(overlap.darmanis[,4])),log2(as.numeric(overlap.allen[,3])))

plot(log2(as.numeric(overlap.darmanis[,2])),log2(as.numeric(overlap.allen[,2])),xlab="Darmanis Microglia signal",ylab="Allen Microglia signal",main="R=0.095")

cor(log2(as.numeric(overlap.darmanis[,2])),log2(as.numeric(overlap.allen[,2])))
