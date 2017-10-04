#Assure in correct directory must contain 'BGI.working.Feb.RData' and 'barreslab_rnaseq.csv' -
setwd('~/Desktop/Neurogenomics Data/CellProportionsProject/Original Data')

### Data to Estimate
load("BGI.working.Feb.RData")
obsv.exprs = rpkm
#################################

profiles <- read.csv('barreslab_rnaseq.csv', header = TRUE)#Reads in profiles from Zhang et al data
amounts = 100 #Amount of genes to be considered for each cell-type.
tmp=match(toupper(gene.symbols),toupper(profiles$Gene.symbol)) ### 11941 genes match
genesofi.bid=c(as.character(profiles$Gene.symbol[genesofi.down]),as.character(profiles$Gene.symbol[genesofi.up]))


## Zhang has 14454 genes
## Stanley data has 
### Which oligo cells should be used here? 
rel.profiles <- cbind(profiles$Astrocytes,profiles$Neuron,profiles$Oligodendrocyte.Precursor.Cell,profiles$Microglia)[tmp[!is.na(tmp)],]#Hard-coded numbers that should be modified to extract profiles of cell-types of interest. 
rel.obsv.exprs <- obsv.exprs[!is.na(tmp),] #Hard-coded numbers that should be modified to extract expression data back from merged datasets. 

plot(apply(rel.obsv.exprs,1,mean),apply(rel.obsv.exprs,1,mean)-apply(rel.profiles,1,mean)) #19.41

## Check that the genes match
match(toupper(gene.symbols)[!is.na(tmp)],toupper(profiles$Gene.symbol)[tmp[!is.na(tmp)]])==1:11941

upregs <- data.frame(NA) #Data frame that will store the degree of upregulation of each gene, for each cell-type in comparison to others.
genesofi <- c() #List that will store index of each gene that is significantly upregulated. 
for(i in 1:ncol(rel.profiles)) {
  for(j in 1:nrow(rel.profiles)){
    upregs[j,i] = rel.profiles[j,i]/mean(rel.profiles[j,-i]) #Calculates degree of upregulation of gene j in cell-type i.
  }
  genesofi = c(genesofi, order(upregs[,i])[(length(upregs[,i])-amounts+1):length(upregs[,i])])  #Appends the indexes of top AMOUNTS (as previously defined) amount of genes specific to the cell-type of interest.
}


####### Selecting only abundant genes ###########################
ind.profiles = rel.profiles[genesofi,] #Reduces profiles down to only of those significantly upregulated. 
ind.obsv.exprs = as.matrix(rel.obsv.exprs[genesofi,1:ncol(rel.obsv.exprs)]+0.1) # + .Machine[[1]] #Reduces expression data down to only of those significantly upregulated. 0.1 is added to prevent errors related to near-0 values. 
gene.names=profiles$Gene.symbol[tmp[!is.na(tmp)]][genesofi]
cell.type=c(rep("Astrocyte",100),rep("Neuron",100),rep("Oligodendrocyte",100),rep("Microglia",100))

max.expr = apply(ind.obsv.exprs, 1, max, na.rm=TRUE) #Calculates max observed expression value for each gene across all areas.
ind.obsv.exprs=ind.obsv.exprs[which(max.expr>2),]
ind.profiles=ind.profiles[which(max.expr>2),] 
cell.type=cell.type[which(max.expr>2)]

#############################################################################                    
astrocyte.ratio=matrix(nrow=length(which(cell.type=="Astrocyte")),ncol=ncol(ind.obsv.exprs))
neuron.ratio=matrix(nrow=length(which(cell.type=="Neuron")),ncol=ncol(ind.obsv.exprs))
oligodendro.ratio=matrix(nrow=length(which(cell.type=="Oligodendrocyte")),ncol=ncol(ind.obsv.exprs))
microglia.ratio=matrix(nrow=length(which(cell.type=="Microglia")),ncol=ncol(ind.obsv.exprs))

for (ii in 1:ncol(ind.obsv.exprs)) {
astrocyte.ratio[,ii]=ind.profiles[which(cell.type=="Astrocyte"),1]/ind.obsv.exprs[which(cell.type=="Astrocyte"),ii]
neuron.ratio[,ii]=ind.profiles[which(cell.type=="Neuron"),2]/ind.obsv.exprs[which(cell.type=="Neuron"),ii]
oligodendro.ratio[,ii]=ind.profiles[which(cell.type=="Oligodendrocyte"),3]/ind.obsv.exprs[which(cell.type=="Oligodendrocyte"),ii]
microglia.ratio[,ii]=ind.profiles[which(cell.type=="Microglia"),4]/ind.obsv.exprs[which(cell.type=="Microglia"),ii]
}

apply(astrocyte.ratio,1,mean)/apply(astrocyte.ratio,1,sd)
apply(neuron.ratio,1,mean)/apply(neuron.ratio,1,sd)
apply(oligodendro.ratio,1,mean)/apply(oligodendro.ratio,1,sd)
apply(microglia.ratio,1,mean)/apply(microglia.ratio,1,sd)

plot((log2(ind.profiles[,2])*log2(apply(ind.obsv.exprs,1,mean)))/2,log2((ind.profiles[,2]/apply(ind.obsv.exprs,1,mean))))

plot(apply(ind.obsv.exprs,1,mean),ind.profiles[,1]) #R=.23
plot(apply(ind.obsv.exprs,1,mean),ind.profiles[,2]) #R=.48
plot(apply(ind.obsv.exprs,1,mean),ind.profiles[,3]) #R=.42
plot(apply(ind.obsv.exprs,1,mean),ind.profiles[,4]) #R=.17
#For just the genes of interest: .11 .33 .06 .04

astrocyte.pca=prcomp(t(ind.obsv.exprs[which(cell.type=="Astrocyte"),]))
plot(astrocyte.pca$x)

neuron.pca=prcomp(t(ind.obsv.exprs[which(cell.type=="Neuron"),]))
plot(neuron.pca$x)

oligodendro.pca=prcomp(t(ind.obsv.exprs[which(cell.type=="Oligodendrocyte"),]))
plot(oligodendro.pca$x)

microglia.pca=prcomp(t(ind.obsv.exprs[which(cell.type=="Microglia"),]))
plot(microglia.pca$x)

apply(astrocyte.ratio,2,mean)
plot(apply(neuron.ratio,1,mean))
plot(apply(oligodendro.ratio,1,mean))
plot(apply(microglia.ratio,1,mean))

sampl.pca=prcomp(t(ind.obsv.exprs))
cbind(sampl.pca$rotation[,1][order(sampl.pca$rotation[,1])],as.character(gene.names[order(sampl.pca$rotation[,1])]))
plot(sampl.pca$rotation,cex.lab=1.5)

astrocyte.pca=prcomp(t(astrocyte.ratio))
cbind(astrocyte.pca$rotation[,1][order(astrocyte.pca$rotation[,1])],as.character(gene.names[order(astrocyte.pca$rotation[,1])]))
plot(astrocyte.pca$rotation,cex.lab=1.5)

neuron.pca=prcomp(t(neuron.ratio))
cbind(neuron.pca$rotation[,1][order(neuron.pca$rotation[,1])],as.character(gene.names[order(neuron.pca$rotation[,1])]))
plot(neuron.pca$rotation,cex.lab=1.5)

oligodendrocyte.pca=prcomp(t(oligodendrocyte.ratio))
cbind(oligodendrocyte.pca$rotation[,1][order(oligodendrocyte.pca$rotation[,1])],as.character(gene.names[order(oligodendrocyte.pca$rotation[,1])]))
plot(oligodendro.pca$rotation,cex.lab=1.5)

microglia.pca=prcomp(t(microglia.ratio))
cbind(microglia.pca$rotation[,1][order(microglia.pca$rotation[,1])],as.character(gene.names[order(microglia.pca$rotation[,1])]))
plot(microglia.pca$rotation,cex.lab=1.5)
######################################################################

######################################################################
guessprops = diag(1, ncol = ncol(ind.obsv.exprs), nrow = ncol(ind.profiles))  #Creates matrix to store first-guess proportions for each cell-type.
guessprops[,] = c(0.2, 0.5, 0.2, 0.1) #Assigns guess values to aforementioned matrix.
guessexprs = ind.profiles%*%guessprops + 0.1  #Calculates guess expression values predicted by guess cell-type proportions. 0.1 added to prevent errors related to near-0 values. 
guessratios = ind.obsv.exprs/guessexprs #Calculates ratio of true expression values to the guess expression values as first-guess of proportionality constants.
ratio = apply(guessratios,1, mean, na.rm=TRUE)
ratiodiag=diag(ratio)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-1.75000 -0.41850 -0.06879 -0.10420  0.25260  1.47500


### Keeps only ratios within the middle quartiles
tmp=which(log(ratio,10)<summary(log(ratio,10))[5] & log(ratio,10)>summary(log(ratio,10))[2])
ratiodiag=diag(ratio[tmp])
ind.obsv.exprs=ind.obsv.exprs[tmp,]
ind.profiles=ind.profiles[tmp,]

#Averages over all guessed ratios for each gene to obtain an estimated scaling factor for the amplification of each gene. 
jpeg("Base-10 logarithm of ratio of observed expression data to projected predicted expression.jpeg")
hist(log(diag(ratio),10), freq=FALSE, xlab='Base-10 logarithm of ratio of observed expression data to projected predicted expression', main='Distribution of estimated ratios for Stanley Institute Data', col='darkgreen')
dev.off()

library(broom) #Library that assists in cleaner data management of linear model objects. To install: install.packages("broom")

gene.means <- apply(ind.obsv.exprs,1,mean, na.rm = TRUE) #Calculates average expression value for each gene in dataset. 
lll <- lm( ind.obsv.exprs ~ -1 + ind.profiles, weights = 1/(.5 + gene.means) ) #Executes weighted linear fitting of data set without diagonal of ratios.
EstPropslll = data.matrix(tidy(coef(lll))[,-1]) #Obtains proportions from model fit. 
EstExprslll = ind.profiles%*%EstPropslll+0.1 #Determines expression values predicted by linear model. 
logRatioslll = log(ind.obsv.exprs/EstExprslll, 10)    #Obtains base-10 logarithm of ratios.

jpeg("Correlation of model fit obsv exprs data.jpeg")
hist(diag(cor(EstExprslll, ind.obsv.exprs, use='complete.obs'))) #Plots correlation of modelfit with observed expression data.
dev.off()

gene.means <- apply(ind.obsv.exprs,1,mean, na.rm = TRUE) #Calculates average expression value for each gene in dataset. 
lll2 <- lm( ind.obsv.exprs ~ -1 + ratiodiag%*%ind.profiles, weights = 1/(.5 + gene.means) ) #Executes weighted linear fitting of data set.
EstPropslll2 = data.matrix(tidy(coef(lll2))[,-1]) #Obtains proportions from model fit. 
EstExprslll2 = ratiodiag%*%ind.profiles%*%EstPropslll2 #Determines expression values predicted by weighted linear model. 
logRatioslll2 = log(ind.obsv.exprs/EstExprslll2, 10) #Determines base-10 logarithm of ratio of estimated to observed expression values. 
hist(diag(cor(EstExprslll2, ind.obsv.exprs, use='complete.obs'))) #Plots correlation of modelfit with observed expression data.
meanvals = apply(ind.obsv.exprs, 1, mean, na.rm = TRUE) #Obtains gene-wise average value of expression data.
Residual_Proportions = apply(abs(lll2[["residuals"]])/meanvals, 1, mean, na.rm = TRUE) #Calculates proportion of residuals with respect to average observed expression values.

jpeg("Residual Proportions Hist.jpeg")
hist(Residual_Proportions) #Plots proportion of residuals with respect to observed value.
dev.off()


##################################
#Excludes genes with "high" proportional residuals to re-do model fitting. 
ind.obsv.exprs = ind.obsv.exprs[-which(Residual_Proportions>0.7),]
ind.profiles = ind.profiles[-which(Residual_Proportions>0.7),]


guessprops = diag(1, ncol = ncol(ind.obsv.exprs), nrow = ncol(ind.profiles))  #Creates matrix to store first-guess proportions for each cell-type.
guessprops[,] = c(0.2, 0.5, 0.2, 0.1) #Assigns guess values to aforementioned matrix.
guessexprs = ind.profiles%*%guessprops + 0.1  #Calculates guess expression values predicted by guess cell-type proportions. 0.1 added to prevent errors related to near-0 values. 
guessratios = ind.obsv.exprs/guessexprs #Calculates ratio of true expression values to the guess exoression values as first-guess of proportionality constants.
ratiodiag = diag(apply(guessratios,1, mean, na.rm=TRUE))  #Averages over all guessed ratios for each gene to obtain an estimated scaling factor for the amplification of each gene. 
gene.means <- apply(ind.obsv.exprs,1,mean, na.rm = TRUE) #Calculates average expression value for each gene in dataset. 
lll2 <- lm( ind.obsv.exprs ~ -1 + ratiodiag%*%ind.profiles, weights = 1/(.5 + gene.means) ) #Executes weighted linear fitting of data set.
EstPropslll2 = data.matrix(tidy(coef(lll2))[,-1]) #Obtains proportions from model fit. 
EstExprslll2 = ratiodiag%*%ind.profiles%*%EstPropslll2 #Determines expression values predicted by weighted linear model. 
logRatioslll2 = log(ind.obsv.exprs/EstExprslll2, 10) #Determines base-10 logarithm of ratio of estimated to observed expression values. 

jpeg("Distribution of estimated ratios.jpeg")
hist(log(diag(ratiodiag),10), freq=FALSE, xlab='Base-10 logarithm of ratio of observed expression data to projected predicted expression.', main='Distribution of estimated ratios for Stanley Institute Data', col='darkgreen')
dev.off()

jpeg("Distribution of corr of fitted data to estimated ratios.jpeg")
hist(diag(cor(EstExprslll2, ind.obsv.exprs, use='complete.obs')), freq=FALSE, xlab='Correlation', main='Distribution of correlations of fitted data to the estimated ratios for Stanley Institute Data', col='darkgreen') #Plots correlation of modelfit with observed expression data.
dev.off()

hist(cor(EstExprslll2, ind.obsv.exprs, use='complete.obs')[which(cor(EstExprslll2, ind.obsv.exprs, use='complete.obs')<1.0)], freq=FALSE, xlab='Correlation', main='Distribution of correlations of fitted data to the estimated ratios for Stanley Institute Data', col='darkgreen') #Plots correlation of modelfit with observed expression data.


meanvals = apply(ind.obsv.exprs, 1, mean, na.rm = TRUE) #Obtains gene-wise average value of expression data.
Residual_Proportions = apply(abs(ind.obsv.exprs-EstExprslll2)/meanvals, 1, mean, na.rm = TRUE) #Calculates proportion of residuals with respect to average observed expression values
hist(Residual_Proportions) #Calculates proportion of residuals with respect to observed value.

ResidualsPerSamples = apply(t(abs(lll2[["residuals"]])), 1, mean, na.rm = TRUE) #Calculates proportion of residuals with respect to average observed expression values by sample.
hist(ResidualsPerSamples) #Plots proportion of residuals with respect to observed value.

Header = c('Astrocytes', 'Neurons', 'Oligodendrocytes', 'Microglia')

jpeg("Fitted Proportions Box Whisker Plot.jpeg")
boxplot(t(EstPropslll2)[which(Group=="SCZ"),], names=Header, main='Fitted Proportions of cell-types for Stanley Brain Institute Data')
dev.off()


apply(EstPropslll2,2,sum)
prop.sum=apply(EstPropslll2,2,sum)
est.prop.crctdz=matrix(nrow=4,ncol=81)
for (i in 1:4){
  est.prop.crctdz[i,]=EstPropslll2[i,]/prop.sum
}

jpeg("Sum of Proportions Stanley Data.jpeg")
hist(apply(EstPropslll2,2,sum),main="Sum of Estimated Proportions for All Four Cell Types",xlab="Sum of the Estimated Proportions")
dev.off()

####Checking for differences between groups 
Group=as.character(Group)
Group[which(Group=="A")]="SCZ"
Group[which(Group=="B")]="CRTL"
Group[which(Group=="C")]="BPD"
Group=as.factor(Group)
colnames(EstPropslll2)=Group
microglia=EstPropslll2[4,c(which(Group=="CRTL"),which(Group=="SCZ"))]
group.microglia=as.factor(as.character(Group[c(which(Group=="CRTL"),which(Group=="SCZ"))]))

jpeg("Fitted Proportions Box Whisker Plot SCZ vs CRTL.jpeg")
boxplot(microglia~group.microglia, main='Fitted Proportions of cell-types for Stanley Brain Institute Data',ylim=c(0,.25))
dev.off()

t.test(microglia~group.microglia)


write.csv(EstPropslll2,"Estimated Cell Proportions of Stanley Samples.csv",quote=FALSE)

