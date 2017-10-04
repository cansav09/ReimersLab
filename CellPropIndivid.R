setwd('~/Desktop/Neurogenomics Data/CellProportionsProject/Original Data')

### Data to Estimate
load("BGI.working.Feb.RData")
obsv.exprs = rpkm
#################################
setwd('~/Desktop/Neurogenomics Data/CellProportionsProject/')
profiles <- read.csv('barreslab_rnaseq.csv', header = TRUE)#Reads in profiles from Zhang et al data
tmp=match(toupper(gene.symbols),toupper(profiles$Gene.symbol)) ### 11941 genes match

rel.profiles <- cbind(profiles$Astrocytes,profiles$Neuron,profiles$Oligodendrocyte.Precursor.Cell,profiles$Microglia)[tmp[!is.na(tmp)],]#Hard-coded numbers that should be modified to extract profiles of cell-types of interest. 
rel.obsv.exprs <- obsv.exprs[!is.na(tmp),] #Hard-coded numbers that should be modified to extract expression data back from merged datasets. 

match(toupper(gene.symbols)[!is.na(tmp)],toupper(profiles$Gene.symbol)[tmp[!is.na(tmp)]])==1:11941
gene.names.oi=toupper(profiles$Gene.symbol)[tmp[!is.na(tmp)]]

cell.fc=rel.profiles[,4]/apply(rel.profiles,1,mean)
cell.fc[which(cell.fc<1)]=-1/cell.fc[which(cell.fc<1)]


cell.type=c("Astrocytes","Neuron","Oligodendrocytes","Microglia")
for(ii in 1:4){

cell.sum=summary(rel.profiles[,ii])  
non.cell.sum=summary(as.vector(rel.profiles[,-ii]))


### Let's identify genes that are high in microglia but nothing else
upgenes=intersect(which(apply(rel.profiles[,-ii],1,max)<non.cell.sum[2]),which(rel.profiles[,ii]>cell.sum[5]))
gene.names.up=gene.names.oi[upgenes]

## Let's identify genes that are low in microglia but high in everything else.
downgenes=intersect(which(apply(rel.profiles[,-ii],1,min)>non.cell.sum[2]),which(rel.profiles[,ii]<cell.sum[2]))
gene.names.down=gene.names.oi[downgenes]

ind.profiles = rel.profiles[c(upgenes,downgenes),] #Reduces profiles down to only of those significantly upregulated. 
ind.obsv.exprs = as.matrix(rel.obsv.exprs[c(upgenes,downgenes),]+0.1) # + .Machine[[1]] #Reduces expression data down to only of those significantly upregulated. 0.1 is added to prevent errors related to near-0 values. 

guessprops = diag(1, ncol = ncol(ind.obsv.exprs), nrow = ncol(ind.profiles))  #Creates matrix to store first-guess proportions for each cell-type.
guessprops[,] = c(0.2, 0.5, 0.2, 0.1) #Assigns guess values to aforementioned matrix.
guessexprs = ind.profiles%*%guessprops + 0.1  #Calculates guess expression values predicted by guess cell-type proportions. 0.1 added to prevent errors related to near-0 values. 
guessratios = ind.obsv.exprs/guessexprs #Calculates ratio of true expression values to the guess exoression values as first-guess of proportionality constants.
ratiodiag = diag(apply(guessratios,1, mean, na.rm=TRUE))  

library(broom) #Library that assists in cleaner data management of linear model objects. To install: install.packages("broom")
ind.profiles.cell=cbind(apply(ind.profiles[,-ii],1,mean),ind.profiles[,ii])
gene.means <- apply(ratiodiag%*%ind.profiles,1,mean, na.rm = TRUE) #Calculates average expression value for each gene in dataset. 
lll <- lm( ind.obsv.exprs ~ -1 + ind.profiles.cell, weights = 1/(.5 + gene.means) ) #Executes weighted linear fitting of data set without diagonal of ratios.
EstPropslll = data.matrix(tidy(coef(lll))[,-1]) #Obtains proportions from model fit. 
EstExprslll = ind.profiles.cell%*%EstPropslll+0.1 #Determines expression values predicted by linear model. 
logRatioslll = log(ind.obsv.exprs/EstExprslll, 10)    #Obtains base-10 logarithm of ratios.

ranked.data=ind.obsv.exprs
for(jj in 1:nrow(ranked.data)){
  ranked.data[jj,]=order(ind.obsv.exprs[jj,])
}
write.csv(ranked.data,file=paste0(cell.type[ii]," Specific Genes Rankings Stanley Samples.csv"))

assign(cell.type[ii],EstPropslll,envir=.GlobalEnv)
assign(paste0(cell.type[ii],".genes"),upgenes,envir=.GlobalEnv)

write.csv(cbind(c(gene.names.up,gene.names.down),ind.profiles.cell),file=paste0(cell.type[ii]," Specific Genes.csv"))


}



plot(apply(ind.obsv.exprs,2,mean),apply(ranked.data,2,mean))
hist(apply(ind.obsv.exprs,2,mean))
hist(apply(ranked.data,2,mean))
