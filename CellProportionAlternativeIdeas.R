
### Looking to figure out #######

setwd('~/Desktop/Neurogenomics Data/CellProportionsProject/')

### Data to Estimate
load("BGI.working.Feb.RData")
obsv.exprs = rpkm
#################################

profiles <- read.csv('barreslab_rnaseq.csv', header = TRUE)#Reads in profiles from Zhang et al data
amounts = 200 #Amount of genes to be considered for each cell-type.
tmp=match(toupper(gene.symbols),toupper(profiles$Gene.symbol)) ### 11941 genes match

## Zhang has 14454 genes
## Stanley data has 
### Which oligo cells should be used here? 
rel.profiles <- cbind(profiles$Astrocytes,profiles$Neuron,profiles$Oligodendrocyte.Precursor.Cell,profiles$Microglia)[!is.na(tmp),] #Hard-coded numbers that should be modified to extract profiles of cell-types of interest. 
rel.obsv.exprs <- obsv.exprs[tmp[!is.na(tmp)],] #Hard-coded numbers that should be modified to extract expression data back from merged datasets. 

## Check that the genes match
## match(toupper(gene.symbols)[!is.na(tmp)],toupper(profiles$Gene.symbol)[tmp[!is.na(tmp)]])==1:11941


###### Taking out the houskpng genes
houskpng.genes=read.table("HousekeepingGeneList.txt",sep="\t")
houskpng.genes=houskpng.genes[c(1,4,7,10,13,16,19,22,25,28,31),1]

tmp=match(houskpng.genes,toupper(profiles$Gene.symbol))
profiles.houskpng <- cbind(profiles$Astrocytes,profiles$Neuron,profiles$Oligodendrocyte.Precursor.Cell,profiles$Microglia)[tmp[!is.na(tmp)],] #Hard-coded numbers that should be modified to extract profiles of cell-types of interest. 
rownames(profiles.houskpng)=houskpng.genes[!is.na(tmp)]

tmp=match(houskpng.genes,toupper(gene.symbols))
obsv.exprs.houskpng <- obsv.exprs[tmp[!is.na(tmp)],] #Hard-coded numbers that should be modified to extract expression data back from merged datasets. 
rownames(obsv.exprs.houskpng)=houskpng.genes[!is.na(tmp)]

tmp=match(rownames(obsv.exprs.houskpng),rownames(profiles.houskpng))

obsv.exprs.houskpng=obsv.exprs.houskpng[!is.na(tmp),]
profiles.houskpng=profiles.houskpng[tmp[!is.na(tmp)],] 


houskpng.genes.ratio=apply(obsv.exprs.houskpng[,-1],1,mean)/profiles.houskpng%*%guessprops[,nrow(profiles.houskpng)]
barplot(t(houskpng.genes.ratio),main="Ratio of Housekeeping Genes")


cell.type.avg=c(mean(profiles$Astrocytes),mean(profiles$Neuron),mean(profiles$Oligodendrocyte.Precursor.Cell),mean(profiles$Microglia))
Header = c('Astrocytes', 'Neurons', 'Oligodendrocytes', 'Microglia')
barplot(c(mean(profiles$Astrocytes),mean(profiles$Neuron),mean(profiles$Oligodendrocyte.Precursor.Cell),mean(profiles$Microglia)),names=Header,main="Avg expression across genes for each cell type")
Header = c('Astrocytes','Astrocytes', 'Neurons','Neurons', 'Oligodendrocytes','Oligodendrocytes', 'Microglia','Microglia')
barplot(c(range(profiles$Astrocytes),range(profiles$Neuron),range(profiles$Oligodendrocyte.Precursor.Cell),range(profiles$Microglia)),names=Header,main="Range expression across cell types")
Header = c('Astrocytes', 'Neurons', 'Oligodendrocytes', 'Microglia')
barplot(c(sd(profiles$Astrocytes),sd(profiles$Neuron),sd(profiles$Oligodendrocyte.Precursor.Cell),sd(profiles$Microglia)),names=Header,main="Sd expression across cell types")

barplot(profiles.houskpng,names=Header)
barplot(profiles.houskpng[2,],names=Header)
barplot(profiles.houskpng[3,],names=Header)
barplot(profiles.houskpng[4,],names=Header)
barplot(profiles.houskpng[5,],names=Header)
barplot(profiles.houskpng[6,],names=Header)



gene.sd=apply(profiles[,c(3,4,5,8)], 1, sd)
gene.avg=apply(profiles[,c(3,4,5,8)], 1, mean)
gene.sd.prop=gene.avg/gene.min
gene.min=apply(profiles[,c(3,4,5,8)], 1, min)

gene.sd.prop=cbind(as.character(profiles$Gene.symbol[order(gene.sd.prop,decreasing=TRUE)]),sort(gene.sd.prop,decreasing=TRUE))
gene.sd.prop.sort=gene.sd.prop[-(grep("Rik",as.character(gene.sd.prop[,1]))),]

cell.type.sd=cbind(profiles[order(cell.type.sd),1],sort(cell.type.sd))
which(is.na(profiles$Gene.symbol))
cor(profiles.houskpng)


### Looking at genes with the least expression in the major cell types 

gene.max=apply(profiles[,c(3,4,5,8)], 1, max)
gene.max=cbind(as.character(profiles$Gene.symbol[order(gene.max)]),sort(gene.max))



upregs <- data.frame(NA) #Data frame that will store the degree of upregulation of each gene, for each cell-type in comparison to others.
genesofi.up <- c() #List that will store index of each gene that is significantly upregulated. 
genesofi.down=c()
sds.other.types=c()
range.other.types=matrix(nrow=nrow(rel.profiles),ncol=2)
for(i in 1:ncol(rel.profiles)) {
  for(j in 1:nrow(rel.profiles)){
    upregs[j,i] = rel.profiles[j,i]/mean(rel.profiles[j,-i]) #Calculates degree of upregulation of gene j in cell-type i.
  sds.other.types[j]=sd(rel.profiles[j,-i])
  range.other.types[j,]=range(rel.profiles[j,-i])
        }
  genesofi.up = c(genesofi.up, order(upregs[,i])[(length(upregs[,i])-amounts+1):length(upregs[,i])])
  genesofi.down=c(genesofi.down,order(upregs[,i],decreasing=TRUE)[(length(upregs[,i])-amounts+1):length(upregs[,i])]) 
    }

####### Selecting only abundant genes ###########################
range.other.types=rbind(range.other.types[genesofi.up,],range.other.types[genesofi.down,])
sds.other.types=rbind(sds.other.types[genesofi.up],sds.other.types[genesofi.down])

max.sd.ratio=c(range.other.types[1:800,2]/sds.other.types[1:800],range.other.types[801:1600,2]/sds.other.types[801:1600])
profiles$Gene.symbol[c(genesofi.up,genesofi.down)]

sd.ratio=cbind(max.sd.ratio,as.character(profiles$Gene.symbol[c(genesofi.up,genesofi.down)]))
 

#Reduces profiles down to only of those significantly upregulated. 
ind.profiles = rbind(rel.profiles[genesofi.up,],rel.profiles[genesofi.down,])
ind.obsv.exprs = as.matrix(rel.obsv.exprs[genesofi,2:ncol(rel.obsv.exprs)]+0.1) # + .Machine[[1]] #Reduces expression data down to only of those significantly upregulated. 0.1 is added to prevent errors related to near-0 values. 

max.expr = apply(ind.obsv.exprs, 1, max, na.rm=TRUE) #Calculates max observed expression value for each gene across all areas.
ind.obsv.exprs=ind.obsv.exprs[which(max.expr>2),]
ind.profiles=ind.profiles[which(max.expr>2),] 


######################################################################
guessprops = diag(1, ncol = ncol(ind.obsv.exprs), nrow = ncol(ind.profiles))  #Creates matrix to store first-guess proportions for each cell-type.
guessprops[,] = c(0.2, 0.5, 0.2, 0.1) #Assigns guess values to aforementioned matrix.
guessexprs = ind.profiles%*%guessprops + 0.1  #Calculates guess expression values predicted by guess cell-type proportions. 0.1 added to prevent errors related to near-0 values. 
guessratios = ind.obsv.exprs/guessexprs #Calculates ratio of true expression values to the guess expression values as first-guess of proportionality constants.
ratio = apply(guessratios,1, mean, na.rm=TRUE)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-1.75000 -0.41850 -0.06879 -0.10420  0.25260  1.47500


### Keeps only ratios within the middle quartiles
tmp=which(log(ratio,10)<summary(log(ratio,10))[5] & log(ratio,10)>summary(log(ratio,10))[2])
ratiodiag=diag(ratio[tmp])
ind.obsv.exprs=ind.obsv.exprs[tmp,]
ind.profiles=ind.profiles[tmp,]
