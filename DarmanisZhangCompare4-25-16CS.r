## Goals: compare the mouse microglia distinctive genes list (From Zhang et al paper) to the human microglia distinctive genes list (From Darmanis et al paper)
## See how much they do or do not overlap with each other and with the m state gene lists. 

# Step 1) Load both lists
# Step 2) See what number of genes show up in both
# Step 3) See how the fold changes correlate
# Step 4) Look into the genes that don't show up on the human list, see why that might be. 
# Step 5) Load the M state lists and find our which genes are in all the lists

#1) Load both microglia distinctive gene lists. Both are sorted in order of the highest fold changes to the lowest fold changes.

### search and replace "microglia" with another "neuron" or "astrocyte" if you want to look at those cells instead.

darmanis=cbind(neuron.genes,neuron.signal,nonneuron.signal,neuron.fc)
zhang=cbind(neuron.genes.zhang,neuron.signal.zhang,nonneuron.signal.zhang,neuron.fc.zhang)

darmanis=darmanis[which(darmanis[,2]>10),] 
zhang=zhang[which(zhang[,2]>10),] 

zhang=zhang[1:500,] 
zhang[,1]=toupper(zhang[,1])
### Let's get rid of genes in the zhang list that don't appear to have any human equivalent
library(org.Hs.eg.db)
EGIDs <- keys(org.Hs.eg.db)
genesymbol.key=unlist(mget(EGIDs, org.Hs.egSYMBOL))# N=56340 genes in this homo sapiens database

###go back and get the full list from Zhang 
#Get rid of zhang genes without human matches 604 genes not matching (by using symbols anyway)
zhang=zhang[which(!is.na(match(zhang[,1],genesymbol.key))),] 

## Now lets make a top 500 genes and cross reference them to zhang
overlap.darmanis=darmanis[!is.na(match(darmanis[,1],zhang[,1])),]
overlap.zhang=zhang[!is.na(match(zhang[,1],darmanis[,1])),]

overlap.darmanis=overlap.darmanis[order(overlap.darmanis[,1]),]
overlap.zhang=overlap.zhang[order(overlap.zhang[,1]),]

match(overlap.zhang[,1],overlap.darmanis[,1])

plot(log2(as.numeric(overlap.zhang[,4])),log2(as.numeric(overlap.darmanis[,4])),xlab="Zhang Fold Change",ylab="Darmanis Fold Change")

plot(log2(as.numeric(overlap.zhang[,2])),log2(as.numeric(overlap.darmanis[,2])),xlab="Zhang signal",ylab="Darmanis signal")

cor(log2(as.numeric(overlap.zhang[,2])),log2(as.numeric(overlap.darmanis[,2])))

### Only using the cell types' expression 
### While we're at it, let's look at how many of these genes are in the M state lists. 
 
m0m1=read.table("M1-M0genelistInfo3-25-16CS.txt",header=TRUE)
m2m0=read.table("M2-M0genelistInfo3-25-16CS.txt",header=TRUE)
m1m2=read.table("M1-M2genelistInfo3-25-16CS.txt",header=TRUE)

mstate=rbind(m0m1,m2m0,m1m2)
mstate$gene.id=as.character(mstate$gene.id)
mstate=mstate[order(mstate$gene.id),]

m0=tapply(as.numeric(mstate$M0.avg.signal),mstate$gene.id,FUN=mean)
m1=tapply(as.numeric(mstate$M1.avg.signal),mstate$gene.id,FUN=mean)
m2=tapply(as.numeric(mstate$M2.avg.signal),mstate$gene.id,FUN=mean)

m.genes=unique(mstate$gene.id)

##How does zhang match to mstate genes?
sum(!is.na(match(m.genes,zhang[,1]))) #N=116 matches in the mouse microglial distinctive list 
zhang.mstate=as.matrix(zhang[match(m.genes,zhang[,1]),])
zhang.mstate=zhang.mstate[!is.na(zhang.mstate[,1]),]

##log2 scale for zhang or otherwise change the comparison

##Is there a particular m state profile that zhang data correlates to better? 
cor(m0[!is.na(match(m.genes,zhang[,1]))],log2(as.numeric(zhang.mstate[,2]))) #R=.15
cor(m1[!is.na(match(m.genes,zhang[,1]))],log2(as.numeric(zhang.mstate[,2]))) #R=.190
cor(m2[!is.na(match(m.genes,zhang[,1]))],log2(as.numeric(zhang.mstate[,2]))) #R=.310

jpeg("M0vsAllZhang4-29-16CS.jpeg")
plot(m0[!is.na(match(m.genes,zhang[,1]))],log2(as.numeric(zhang.mstate[,2])),xlab="Profile M0 gene signal",ylab="Zhang M0 gene signal",sub="R = .15")
dev.off()
jpeg("M1vsAllZhang4-29-16CS.jpeg")
plot(m1[!is.na(match(m.genes,zhang[,1]))],log2(as.numeric(zhang.mstate[,2])),xlab="Profile M1 gene signal",ylab="Zhang M0 gene signal",sub="R = 0.190") 
dev.off()
jpeg("M2vsAllZhang4-29-16CS.jpeg")
plot(m2[!is.na(match(m.genes,zhang[,1]))],log2(as.numeric(zhang.mstate[,2])),xlab="Profile M2 gene signal",ylab="Zhang M0 gene signal",sub="R = 0.310") 
dev.off()

##How does darmanis match to m state genes? 
sum(!is.na(match(m.genes,darmanis[,1]))) #N=218 matches in the mouse microglial distinctive list 
darmanis.mstate=as.matrix(darmanis[match(m.genes,darmanis[,1]),])
darmanis.mstate=darmanis.mstate[!is.na(darmanis.mstate[,1]),]

##Is there a particular m state profile that darmanis data correlates to better? 
cor(m0[!is.na(match(m.genes,darmanis[,1]))],log2(as.numeric(darmanis.mstate[,2]))) #R=.164
cor(m1[!is.na(match(m.genes,darmanis[,1]))],log2(as.numeric(darmanis.mstate[,2]))) #R=-.086
cor(m2[!is.na(match(m.genes,darmanis[,1]))],log2(as.numeric(darmanis.mstate[,2]))) #R=.071

jpeg("M0vsDarmanis4-29-16CS.jpeg")
plot(m0[!is.na(match(m.genes,darmanis[,1]))],log2(as.numeric(darmanis.mstate[,2])),xlab="Profile M0 gene signal",ylab="Darmanis M0 gene signal",sub="R = 0.164") 
dev.off()
jpeg("M1vsDarmanis4-29-16CS.jpeg")
plot(m1[!is.na(match(m.genes,darmanis[,1]))],log2(as.numeric(darmanis.mstate[,2])),xlab="Profile M1 gene signal",ylab="Darmanis M0 gene signal",sub="R = 0-.086") 
dev.off()
jpeg("M2vsDarmanis4-29-16CS.jpeg")
plot(m2[!is.na(match(m.genes,darmanis[,1]))],log2(as.numeric(darmanis.mstate[,2])),xlab="Profile M2 gene signal",ylab="Darmanis M0 gene signal",sub="R = .071") 
dev.off()


###These don't really match that well 
### The human lists do match up better (suppose that's not a suprise) but does this indicate that the 
## darmanis list is more representative of microglial phenotype or is this more just a byproduct of gene name annotation matching better between human data. 
overlap.mstate=zhang.mstate[match(darmanis.mstate[,1],zhang.mstate[,1]),1]
overlap.mstate=overlap.mstate[!is.na(overlap.mstate)]#N=95 genes that are in all lists this is with no expression cutoffs besides there actually being at least some sort of signal

overlap.mstate.m0=m0[!is.na(match(m.genes,overlap.mstate))]
overlap.mstate.m1=m1[!is.na(match(m.genes,overlap.mstate))]
overlap.mstate.m2=m2[!is.na(match(m.genes,overlap.mstate))]

darmanis.mstate=darmanis.mstate[match(overlap.mstate,darmanis.mstate[,1]),c(2,4)]
zhang.mstate=zhang.mstate[match(overlap.mstate,zhang.mstate[,1]),c(2,4)]

names(zhang.mstate)[1:62]=zhang[match(names(zhang.mstate)[1:62],rownames(zhang)),1]
match(names(zhang.mstate)[1:62],names(darmanis.mstate)[1:62])==1:62 #checking that the genes are in the right order

jpeg("DarmanisvsM0Signalw954-29-16CS.jpeg")
plot(overlap.mstate.m0,as.numeric(darmanis.mstate[,1]), ylab="Darmanis RPKM",xlab="M0 Signal Martinez et al",sub="R=.2578") #R=.2578
dev.off()
jpeg("DarmanisvsM1Signalw954-29-16CS.jpeg")
plot(overlap.mstate.m1,as.numeric(darmanis.mstate[,1]),ylab="Darmanis RPKM",xlab="M1 Signal Martinez et al",sub="R=-.19") #R=-.19
dev.off()
jpeg("DarmanisvsM2Signalw954-29-16CS.jpeg")
plot(overlap.mstate.m2,as.numeric(darmanis.mstate[,1]),ylab="Darmanis RPKM",xlab="M2 Signal Martinez et al",sub="R=-.0013") #R=-.0013
dev.off()
jpeg("ZhangvsM0Signalw954-29-16CS.jpeg")
plot(overlap.mstate.m0,as.numeric(zhang.mstate[,1]),ylab="Zhang FPKM",xlab="M0 Signal Martinez et al",sub="R=-.057") #R=-.057
dev.off()
jpeg("ZhangvsM1Signalw954-29-16CS.jpeg")
plot(overlap.mstate.m1,as.numeric(zhang.mstate[,1]),ylab="Zhang FPKM",xlab="M1 Signal Martinez et al",sub="R=.396") #R=.396
dev.off()
jpeg("ZhangvsM2Signalw954-29-16CS.jpeg")
plot(overlap.mstate.m2,as.numeric(zhang.mstate[,1]),ylab="Zhang FPKM",xlab="M2 Signal Martinez et al",sub="R=.330") #R=.330
dev.off()

### Printing this out into a list so I can take a better look at these genes. 
all=cbind(overlap.mstate,darmanis.mstate,zhang.mstate,overlap.mstate.m0,overlap.mstate.m1,overlap.mstate.m2)

write.table(all,"m.state.genes.microglia.distinct.overlap.signals4-27-16CS.csv",sep=",",quote=FALSE,row.names=FALSE,col.names=c("Gene Symbol","Darmanis RPKM","Darmanis Fold Change","Zhang FPKM","Zhang Fold Change","M0 macrophage avg signal","M1 macrophage avg signal","M2 macrophage avg signal"))
