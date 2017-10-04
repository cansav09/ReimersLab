## Step 1) Load original GEO data files
library(Biobase)
library(GEOquery)
library(beepr)
library(e1071)
library(DescTools)
library(broom) 

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



###################################################################################################
##############################Estimate Cell Proportions ###################################

### Data to Estimate
setwd(input)
obsv.exprs = brain.region.sum
gene.symbols=unique(gene.name.webb)
#################################
# Import list of validated cell-specific genes 

profiles <- read.csv('SupplTable1_CellTypeSpecificGenes_Master3Sorted.csv', header = TRUE)#Reads in profiles from Zhang et al data
profiles=profiles[-which(as.character(profiles$CellType_Primary)=="Mural"),]

cell.type.edit.stan=cell.type[-unique(c(which(out.pr.gene>0),which(abs(gene.skew)>.7),which(abs(gene.kurt)>3)))]

# Match the genes in the Hagenauer list to the genes in the Stanley Data set
#tmp=match(toupper(gene.symbols),toupper(profiles$Gene.Symbol..Human.)) 
tmp=match(toupper(gene.symbols),toupper(cell.type.genes)) 

rel.profiles <- profiles[tmp[!is.na(tmp)],]
rel.obsv.exprs <- obsv.exprs[!is.na(tmp),] 
cell.type.edit.stan=cell.type.edit.stan[tmp[!is.na(tmp)]]


######################################################################

gene.means <- apply(rel.obsv.exprs,1,mean, na.rm = TRUE) 
gene.medians <- apply(rel.obsv.exprs,1,median, na.rm = TRUE)  

# How centered are these genes? 
plot(gene.means/gene.medians)

xx=which(gene.means/gene.medians>10)
rel.profiles <- rel.profiles[-xx,]
rel.obsv.exprs <- rel.obsv.exprs[-xx,] 
cell.type.edit.stan=cell.type.edit.stan[-xx]


# Calculate the average expression value for each cell type geneset for each sample
samp.gene.means=matrix(ncol=5,nrow=ncol(rel.obsv.exprs)) 
for(ii in 1:ncol(rel.obsv.exprs)){
  samp.gene.means[ii,] <- tapply(rel.obsv.exprs[,ii],cell.type.edit.stan,mean) 
}

# Plot the average cell type geneset expression values by cell type
setwd(output)
jpeg("Boxplot of sample means for cell spec genesets.jpeg")
boxplot(samp.gene.means,names=sort(unique(cell.type)),las=2,cex.names=.7)
dev.off()

# For each cell type attempt to find "housekeeping genes" based on their distributions
# This means that they would be genes' whose expression is fairly normally distributed
# They also would not be too skewed or kurtotic

out.pts=matrix(ncol=2,nrow=0)
gene.skew=1:nrow(rel.obsv.exprs)
gene.kurt=gene.skew=1:nrow(rel.obsv.exprs)
for(ii in 1:nrow(rel.obsv.exprs)){
  gene.skew[ii]=skewness(rel.obsv.exprs[ii,yy])
  gene.kurt[ii]=kurtosis(rel.obsv.exprs[ii,yy],type=3)
  xx=which(rel.obsv.exprs[ii,] %in% Outlier(rel.obsv.exprs[ii,], na.rm=TRUE))
  xx=cbind(rep(ii,length(xx)),xx)
  out.pts=rbind(out.pts,xx)
  }
# This prints out a table of the the number of outliers in each gene of all the genes that have at least one outlier
out.pr.gene=table(profiles$Gene.Symbol..Human.[out.pts[,1]])
hist(out.pr.gene)
summary(out.pr.gene)

hist(gene.skew)
summary(gene.skew)

hist(gene.kurt)
summary(gene.kurt)

# Let's rid ourselves of the genes whose data 1) has any outliers, 2) has a skew, and 3) Is too platykurtotic, and maybe even ones that are too leptokurtotic
#rel.obsv.exprs.edit=rel.obsv.exprs[-unique(c(which(out.pr.gene>3),which(abs(gene.skew)>.7),which(abs(gene.kurt)>1))),]
rel.obsv.exprs.edit=rel.obsv.exprs[-unique(c(which(abs(gene.skew)>1.5),which(abs(gene.kurt)>1.5))),]
#cell.type.edit=cell.type[-unique(c(which(out.pr.gene>3),which(abs(gene.skew)>.5),which(abs(gene.kurt)>1)))]
cell.type.edit.stan=cell.type.edit.stan[-unique(c(which(abs(gene.skew)>1.5),which(abs(gene.kurt)>1.5)))]

# Let's recalculate the average expression value for each cell type geneset for each sample
# after having removed the funky distribution genes
samp.gene.means.edit=matrix(ncol=5,nrow=ncol(rel.obsv.exprs.edit))
for(ii in 1:ncol(rel.obsv.exprs.edit)){
  samp.gene.means.edit[ii,] <- tapply(rel.obsv.exprs.edit[,ii],cell.type.edit.stan,mean) 
}


gene.means.edit <- apply(rel.obsv.exprs.edit,1,mean, na.rm = TRUE)  
gene.medians.edit <- apply(rel.obsv.exprs.edit,1,median, na.rm = TRUE) 

# Let's re-plot the average cell type geneset expression values by cell type
jpeg("Boxplot of sample means for cell spec genesets edited gene lists.jpeg")
boxplot(samp.gene.means.edit,names=sort(unique(cell.type.edit)),las=2,cex.names=.7)
dev.off()

# We'll convert these to proportions for each sample
samp.cell.type.mean=apply(samp.gene.means.edit,1,sum)
samp.gene.means.edit2=matrix(ncol=5,nrow=ncol(rel.obsv.exprs.edit))
for(ii in 1:ncol(rel.obsv.exprs.edit)){
  samp.gene.means.edit2[ii,] <- samp.gene.means.edit[ii,]/samp.cell.type.mean[ii]
}

# Make a boxplot of the proportions
jpeg("Boxplot of proportions for each sample before regression.jpeg")
boxplot(samp.gene.means.edit2,names=sort(unique(cell.type.edit)),las=2,cex.names=.7)
dev.off()

## we'll use the average for each gene as the "consensus or comparison profile"
rel.profiles=apply(rel.obsv.exprs.edit,1,mean)

# use Regression to identify to what extent each sample fits to the consensus cell type profile
adjust.props=samp.gene.means.edit2 # object to store the adjusted proportions in 
props.resid.lll=1:length(gene.means.edit) # object to store the residuals of the proportions in. 

for (ii in 1:length(unique(cell.type.edit))){ # This set up only uses the genes specific to the cell type being estimated to estimate the proportion
  xx=which(cell.type.edit==unique(cell.type.edit)[ii]) # index of the data of the cell type being evaluated
  lll <- lm( rel.obsv.exprs.edit[xx,] ~ -1 + rel.profiles[xx], weights = 1/(.5 + gene.means.edit[xx]) ) # weighted by gene means
  adjust.props[,ii]=coef(lll)*apply(samp.gene.means.edit2,2,mean)[ii] 
  props.resid.lll[xx]=apply(abs(lll[["residuals"]])/gene.means.edit[xx], 1, mean, na.rm = TRUE) # Adjust the residuals in proportion to their gene means
}

# Make a boxplot of the proportions of hte adjusted proportions
jpeg("Fitted Proportions for Cell Types - Hagenauer List.jpeg")
boxplot(adjust.props,names=sort(unique(cell.type.edit)),las=2,cex.names=.7)
dev.off()

table(cell.type.edit)
#Identifier genes for each cell type
#Astrocyte     Endothelial       Microglia          Neuron Oligodendrocyte 
#30              28              25             138              83 

jpeg("Residual Proportions Hist -Hagenauer List.jpeg") 
hist(props.resid.lll) #Plots proportion of residuals with respect to observed value.
dev.off()

hist(ResidualsPerSamples) #Plots proportion of residuals with respect to observed value.


####Checking for differences between SCZ and CTL groups 
Group=as.character(Group)
Group[which(Group=="A")]="SCZ"
Group[which(Group=="B")]="CRTL"
Group[which(Group=="C")]="BPD"
Group2=Group[-which(Group=="BPD")]
Group2=as.factor(Group2)

length(Group2)
length(adjust.props[,2])

setwd(output)
for(ii in 1:length(unique(cell.type.edit))){
  jpeg(paste0("Fitted Proportions ",sort(unique(cell.type.edit))[ii]," Box Whisker Plot SCZ vs CRTL.jpeg"))
  xx=t.test(adjust.props[,ii][-which(Group=="BPD")]~Group2)
  boxplot(adjust.props[,ii][-which(Group=="BPD")]~Group2, main=paste0("Proportions of ",sort(unique(cell.type.edit))[ii]," t=",round(xx$statistic,3)," p=",round(xx$p.value,3)))
  dev.off()
}

setwd(output)
write.csv(adjust.props,"Estimated Cell Proportions of Stanley Samples.csv",quote=FALSE)


