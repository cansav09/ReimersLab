## Step 1) Load original GEO data files
library(Biobase)
library(GEOquery)
library(beepr)

input='~/Desktop/Neurogenomics Data/CellProportionsProject/Original Data/GSE68559_RAW'
output='~/Desktop/Neurogenomics Data/CellProportionsProject/WebbBrainRegionData'

setwd(input)
fnames=dir()
fnames=grep("expr.txt.gz",fnames,value=TRUE)

tmp= as.vector(read.table(fnames[1],fill=TRUE,nrows=1))
gene.name.webb=read.table(fnames[1])[-1,5]
brain.region = matrix(nr=length(gene.name.webb),nc=length(fnames))

rownames(brain.region)=as.vector(read.table(fnames[1])[-1,1])
for (ii in 1:length(fnames)) {
  xx= read.table(fnames[ii],fill=TRUE,skip=1)
  brain.region[,ii] <- as.numeric(xx[,10])
}
beep()
colnames(brain.region)=unlist(gsub("_isoforms_expr.txt.gz","",fnames))

###### Summing the FPKMS for all transcripts of a particular gene ########################
brain.region.sum=matrix(nrow=length(unique(gene.name.webb)),ncol=ncol(brain.region))
for(ii in 27892:length(unique(gene.name.webb))){
  xx=which(gene.name.webb==unique(gene.name.webb)[ii])
  if(length(xx)>1){
    brain.region.sum[ii,]=apply(brain.region[xx,],2,sum)
  }else{
    brain.region.sum[ii,]=brain.region[xx,]
  }
}
beepr()

saveRDS(brain.region.sum,file="Brain Region Webb Isoform Sum.RDS")
saveRDS(brain.region,file="Brain Region Webb Individ Isoforms.RDS")

########## Set the names for all the samples ############################################

setwd(output)

############### Download SOFT file from GEO #########################
setwd(input)

gse="GSE68559"
if(file.exists(paste0(gse,".soft"))==FALSE){
  getGEOfile(gse,destdir=input)
  gunzip("GSE68559.soft.gz")
}
geo.soft <- getGEO(filename=paste0(gse,".soft"))
## This part may take some time depending on how many samples and how large the file is. 

############### Get sample info from GEO soft files##################
geo.sample.info=matrix(nrow=length(geo.soft@header$sample_id),ncol=7)
for (ii in 1:length(geo.soft@header$sample_id)){
  geo.sample.info[ii,1:2]=c(geo.soft@header$sample_id[ii],geo.soft@gsms[[1]]@header$title[ii])
}
geo.groups=1:length(geo.sample.info[,2])

sample.info=matrix(nrow=length(geo.soft@header$sample_id),ncol=2)
for (ii in 1:length(geo.soft@header$sample_id)){
  sample.info[ii,1]=geo.soft@header$sample_id[ii]
}
sample.info[,2]=geo.soft@gsms[[1]]@header$source_name_ch1

sampleID=unlist(strsplit(dimnames(brain.region)[[2]],"_"))
sampleID=sampleID[grep("GSM",sampleID)]

participant=unlist(strsplit(dimnames(brain.region)[[2]],"_"))
participant=participant[grep("^MB*",participant)]
summary(as.factor(participant))
#MB011 MB052 MB059 MB100 MB147 MB148 MB151 MB160 MB197 MB202 
#10    10     9    10    10    10    10     9    10    10 

summary(as.factor(sample.info[,2]))
#amygdala BA10(prefrontal cort)  BA22(temporal/wernicke's) BA24(Ant cingulate)  BA46(dorsolateral prefrontal cortex) cerebellum 
#10                10                   10                          10                      9                           10 
#hippocampus    insula posterior    putamen     raphae nuclei 
#10                10                10                 9 




####Estimate Cell Proportion
###############################################################################################
#Assure in correct directory must contain 'BGI.working.Feb.RData' and 'barreslab_rnaseq.csv' -
setwd('~/Desktop/Neurogenomics Data/CellProportionsProject/')

### Data to Estimate
obsv.exprs = brain.region.sum
#################################

profiles <- read.csv('barreslab_rnaseq.csv', header = TRUE)#Reads in profiles from Zhang et al data
amounts = 100 #Amount of genes to be considered for each cell-type.
tmp=match(toupper(profiles$Gene.symbol),unique(toupper(gene.name.webb)))### 11941 genes match

## Zhang has 14454 genes
## Stanley data has 
### Which oligo cells should be used here? 
rel.profiles <- cbind(profiles$Astrocytes,profiles$Neuron,profiles$Oligodendrocyte.Precursor.Cell,profiles$Microglia)[!is.na(tmp),] #Hard-coded numbers that should be modified to extract profiles of cell-types of interest. 
rel.obsv.exprs <- obsv.exprs[tmp[!is.na(tmp)],] #Hard-coded numbers that should be modified to extract expression data back from merged datasets. 

## Check that the genes match
match(unique(toupper(gene.name.webb))[tmp[!is.na(tmp)]],toupper(profiles$Gene.symbol)[!is.na(tmp)])==1:11941
genes.names.br=unique(toupper(gene.name.webb))[tmp[!is.na(tmp)]]

#upregs <- data.frame(NA) #Data frame that will store the degree of upregulation of each gene, for each cell-type in comparison to others.
##genesofi <- c() #List that will store index of each gene that is significantly upregulated. 
#for(i in 1:ncol(rel.profiles)) {
#  for(j in 1:nrow(rel.profiles)){
#   upregs[j,i] = rel.profiles[j,i]/mean(rel.profiles[j,-i]) #Calculates degree of upregulation of gene j in cell-type i.
#  }
#  genesofi = c(genesofi, order(upregs[,i])[(length(upregs[,i])-amounts+1):length(upregs[,i])])  #Appends the indexes of top AMOUNTS (as previously defined) amount of genes specific to the cell-type of interest.
#}

####### Selecting only abundant genes ###########################

ind.profiles = rel.profiles[genesofi,] #Reduces profiles down to only of those significantly upregulated. 
ind.obsv.exprs = as.matrix(rel.obsv.exprs[genesofi,1:ncol(rel.obsv.exprs)]+0.1) # + .Machine[[1]] #Reduces expression data down to only of those significantly upregulated. 0.1 is added to prevent errors related to near-0 values. 
gene.names.br=profiles$Gene.symbol[!is.na(tmp)][genesofi]

max.expr = apply(ind.obsv.exprs, 1, max, na.rm=TRUE) #Calculates max observed expression value for each gene across all areas.
ind.obsv.exprs=ind.obsv.exprs[which(max.expr>2),]
ind.profiles=ind.profiles[which(max.expr>2),] ## Only gets rid of 8 genes


######################################################################
guessprops = diag(1, ncol = ncol(ind.obsv.exprs), nrow = ncol(ind.profiles))  #Creates matrix to store first-guess proportions for each cell-type.
guessprops[,] = c(0.2, 0.5, 0.2, 0.1) #Assigns guess values to aforementioned matrix.
guessexprs = ind.profiles%*%guessprops + 0.1  #Calculates guess expression values predicted by guess cell-type proportions. 0.1 added to prevent errors related to near-0 values. 
guessratios = ind.obsv.exprs/guessexprs #Calculates ratio of true expression values to the guess expression values as first-guess of proportionality constants.
ratio = apply(guessratios,1, mean, na.rm=TRUE)

### Ratios from Stanley data that was originally used
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-1.75000 -0.41850 -0.06879 -0.10420  0.25260  1.47500

### Ratios from Brain Region Webb et al data
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.00100  0.04316  0.14210  0.51030  0.42720 12.73000 


### Keeps only ratios within the middle quartiles
tmp=which(log(ratio,10)<summary(log(ratio,10))[5] & log(ratio,10)>summary(log(ratio,10))[2])
ratiodiag=diag(ratio[tmp])
ind.obsv.exprs=ind.obsv.exprs[tmp,]
ind.profiles=ind.profiles[tmp,]
gene.names=gene.names[tmp]
#Averages over all guessed ratios for each gene to obtain an estimated scaling factor for the amplification of each gene. 
jpeg("Base-10 logarithm of ratio of observed expression data to projected predicted expression.jpeg")
hist(log(diag(ratio),10), freq=FALSE, xlab='Base-10 logarithm of ratio of observed expression data to projected predicted expression', main='Distribution of estimated ratios for Brain Region data', col='darkgreen')
dev.off()

library(broom) #Library that assists in cleaner data management of linear model objects. To install: install.packages("broom")

gene.means <- apply(ind.obsv.exprs,1,mean, na.rm = TRUE) #Calculates average expression value for each gene in dataset. 
lll <- lm( ind.obsv.exprs ~ -1 + ind.profiles) #Executes weighted linear fitting of data set without diagonal of ratios.
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
lll2.nowght <- lm( ind.obsv.exprs ~ -1 + ratiodiag%*%ind.profiles ) #Executes weighted linear fitting of data set.
EstPropslll2.nowght = data.matrix(tidy(coef(lll2))[,-1]) #Obtains proportions from model fit. 
EstExprslll2 = ratiodiag%*%ind.profiles%*%EstPropslll2 #Determines expression values predicted by weighted linear model. 
logRatioslll2 = log(ind.obsv.exprs/EstExprslll2, 10) #Determines base-10 logarithm of ratio of estimated to observed expression values. 

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
boxplot(t(EstPropslll2), names=Header, main='Fitted Proportions of cell-types for Stanley Brain Institute Data',ylim=c(0,1))
dev.off()

wacky.samples=c()
for(ii in 1:4){
 wacky.samples= c(wacky.samples,which(EstPropslll2[ii,]>.9))
}
EstProps.sort=EstPropslll2[,-wacky.samples]
sample.info[-wacky.samples,]

prop.sum=apply(EstPropslll2,2,sum)
est.prop.crctd=matrix(nrow=4,ncol=98)
for (i in 1:4){
  est.prop.crctd[i,]=EstPropslll2[i,]/prop.sum
  }
###apply(est.prop.crctd,2,sum)
boxplot(t(est.prop.crctd), names=Header, main='Fitted Proportions of cell-types with correction',ylim=c(0,1))

jpeg("Sum of Proportions Brain Region Data.jpeg")
hist(apply(EstPropslll2,2,sum),main="Sum of Estimated Proportions for All Four Cell Types",xlab="Sum of the Estimated Proportions",breaks=50)
dev.off()


#[,1]         [,2]         
#[1,] "GSM1675509" "BA46"       
#[2,] "GSM1675517" "amygdala"   
#[3,] "GSM1675523" "BA22"       
#[4,] "GSM1675536" "hippocampus"
#[5,] "GSM1675559" "BA22"       
#[6,] "GSM1675568" "BA22"       
#[7,] "GSM1675586" "BA24" 

#amygdala BA10(prefrontal cort)  BA22(temporal/wernicke's) BA24(Ant cingulate)  BA46(dorsolateral prefrontal cortex) cerebellum 
#10                10                   10                          10                      9                           10 
#hippocampus    insula posterior    putamen     raphae nuclei 
#10                10                10                 9 

### There isn't showing any differences between brain regions
anova.props=lm(t(EstProps.sort)~participant[-wacky.samples])
anova.props=aov(lm(t(est.prop.crctd)~sample.info[,2]))
summary(anova.props)

### There are significant differences between participants
anova.props=lm(t(EstProps.sort)~participant[-wacky.samples])
anova.props=aov(lm(t(est.prop.crctd)[-wacky.samples,]~participant[-wacky.samples]))
summary(anova.props)

region=as.factor(sample.info[,2])[-wacky.samples]
subject=as.factor(participant)[-wacky.samples]
df=data.frame(apply(est.prop.crctd[,-wacky.samples]),1,as.numeric),subject,region)
colnames(df)=c(Header,"Subject", "BrainRegion")

anova.props=lm(df$Astrocytes~df$BrainRegion + df$Subject +df$BrainRegion:df$Subject,data=df)

##options(contrasts=c("contr.sum","contr.poly

library(nlme)
library(lme4)
rep.anova=lmer(EstPropslll[1,] ~ region*subject + (1|subject),data = df)

summary(rep.anova)

brain.area.means=matrix(nrow=4,ncol=length(unique(sample.info[-wacky.samples,2])))
brain.area.ses=matrix(nrow=4,ncol=length(unique(sample.info[-wacky.samples,2])))
for( ii in 1:4){
  brain.area.means[ii,]=tapply(EstProps.sort[ii,],sample.info[-wacky.samples,2],mean)
  x=tapply(EstProps.sort[ii,],sample.info[-wacky.samples,2],sd)
  brain.area.ses[ii,]=x/sqrt(table(sample.info[-wacky.samples,2]))
  }

colnames(brain.area.means)=c("dorsolateral prefrontal","prefrontal","temporal/wernicke's","anterior cingulate",unique(sample.info[,2])[-c(1:4)])
colnames(brain.area.ses)=unique(sample.info[-wacky.samples,2])

participant.area.means=matrix(nrow=4,ncol=length(unique(participant[-wacky.samples])))
participant.area.ses=matrix(nrow=4,ncol=length(unique(participant[-wacky.samples])))
for( ii in 1:4){
  participant.area.means[ii,]=tapply(EstProps.sort[ii,],participant[-wacky.samples],mean)
  x=tapply(EstProps.sort[ii,],participant[-wacky.samples],sd)
  participant.area.ses[ii,]=x/sqrt(table(participant[-wacky.samples]))
  }
colnames(participant.area.means)=unique(participant[-wacky.samples])
colnames(participant.area.ses)=unique(participant[-wacky.samples])


jpeg("Proportions Astrocytes Across Brain Regions.jpeg")
par(mai=c(2.8,1,1,1))
xx=barplot(brain.area.means[1,],main= "Proportion of Astrocytes",las=2,ylim=c(0,1))
segments(xx, brain.area.means[1,] - brain.area.ses[1,] * 2, xx,brain.area.means[1,] + brain.area.ses[1,]* 2, lwd = 1.5)
arrows(xx, brain.area.means[1,] - brain.area.ses[1,] * 2, xx,brain.area.means[1,] + brain.area.ses[1,] * 2 * 2, lwd = 1.5, angle = 90,code = 3, length = 0.05)
dev.off()
jpeg("Proportions Neurons Across Brain Regions.jpeg")
par(mai=c(2.8,1,1,1))
barplot(brain.area.means[2,],main= "Proportion of Neurons",las=2,ylim=c(0,1))
segments(xx, brain.area.means[2,] - brain.area.ses[2,] * 2, xx,brain.area.means[2,] + brain.area.ses[2,]* 2, lwd = 1.5)
arrows(xx, brain.area.means[2,] - brain.area.ses[2,] * 2, xx,brain.area.means[2,] + brain.area.ses[2,] * 2 * 2, lwd = 1.5, angle = 90,code = 3, length = 0.05)
dev.off()
jpeg("Proportions Oligodendrocytes Across Brain Regions.jpeg")
par(mai=c(2.8,1,1,1))
barplot(brain.area.means[3,],main= "Proportion of Oligodendrocytes",las=2,ylim=c(0,1))
segments(xx, brain.area.means[3,] - brain.area.ses[3,] * 2, xx,brain.area.means[3,] + brain.area.ses[3,]* 2, lwd = 1.5)
arrows(xx, brain.area.means[3,] - brain.area.ses[3,] * 2, xx,brain.area.means[3,] + brain.area.ses[3,] * 2 * 2, lwd = 1.5, angle = 90,code = 3, length = 0.05)
dev.off()
jpeg("Proportions Microglia Across Brain Regions.jpeg")
par(mai=c(2.8,1,1,1))
barplot(brain.area.means[4,],main= "Proportion of Microglia",las=2,ylim=c(0,1))
segments(xx, brain.area.means[4,] - brain.area.ses[4,] * 2, xx,brain.area.means[4,] + brain.area.ses[4,]* 2, lwd = 1.5)
arrows(xx, brain.area.means[4,] - brain.area.ses[4,] * 2, xx,brain.area.means[4,] + brain.area.ses[4,] * 2 * 2, lwd = 1.5, angle = 90,code = 3, length = 0.05)
dev.off()


jpeg("Fitted Proportions Brain Region Data.jpeg")
par(mai=c(1,1,1,1))
boxplot(t(EstProps.sort), names=Header, main='Fitted Proportions of cell-types for Webb et al Brain Region Data',ylim=c(0,1))
dev.off()


jpeg("Proportions Across Participant Brain Regions.jpeg")
par(mai=c(2.8,1,1,1))
barplot(participant.area.means[4,],main= "Proportions Across Participant Brain Regions",las=2,ylim=c(0,1))
segments(xx, participant.area.means[4,] - participant.area.ses[4,] * 2, xx,participant.area.means[4,] + participant.area.ses[4,]* 2, lwd = 1.5)
arrows(xx, participant.area.means[4,] - participant.area.ses[4,] * 2, xx,participant.area.means[4,] + participant.area.ses[4,] * 2 * 2, lwd = 1.5, angle = 90,code = 3, length = 0.05)
dev.off()



sample.info[-wacky.samples,2]
barplot(sample.info[-wacky.samples,1][1:5],EstProps.sort[,1:5])

participant

df=df[-wacky.samples,]
library(ggplot2)
jpeg("Microglia across brain regions.jpeg")
ggplot(data = df, aes(x = factor(df$BrainRegion), y = as.numeric(df$Microglia), colour = df$Subject)) +       
  geom_line(aes(group = df$Subject)) + geom_point()+ labs(x="Brain Regions",y="Microglia Proportion",colour="Subject")+theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
dev.off()

jpeg("Astrocyte across brain regions.jpeg")
ggplot(data = df, aes(x = factor(df$BrainRegion), y = as.numeric(df$Astrocytes), colour = df$Subject)) +       
  geom_line(aes(group = df$Subject)) + geom_point()+ labs(x="Brain Regions",y="Astrocyte Proportion",colour="Subject")+theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
dev.off()

jpeg("Neuron across brain regions.jpeg")
ggplot(data = df, aes(x = factor(df$BrainRegion), y = as.numeric(df$Neurons), colour = df$Subject)) +       
  geom_line(aes(group = df$Subject)) + geom_point()+ labs(x="Brain Regions",y="Neuron Proportion",colour="Subject")+theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
dev.off()

jpeg("Oligodendrocyte across brain regions.jpeg")
ggplot(data = df, aes(x = factor(df$BrainRegion), y = as.numeric(df$Oligodendrocytes), colour = df$Subject)) +       
  geom_line(aes(group = df$Subject)) + geom_point()+ labs(x="Brain Regions",y="Oligodendrocyte Proportion",colour="Subject")+theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
dev.off()



