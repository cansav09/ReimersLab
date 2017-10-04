### Estimating Cell Proportions using Megan's List
library(e1071)
library(DescTools)
library(broom) 

input='~/Desktop/Neurogenomics Data/CellProportionsProject/Original Data'
output='~/Desktop/Neurogenomics Data/CellProportionsProject/CellProportionHagenauerList'

### Data to Estimate
setwd(input)
load("BGI.working.Feb.RData")
Group=as.character(Group)
Group[which(Group=="A")]="SCZ"
Group[which(Group=="B")]="CTL"
Group[which(Group=="C")]="BPD"
Group2=as.factor(Group2)

obsv.exprs = rpkm
#################################
# Import list of validated cell-specific genes 

profiles <- read.csv('SupplTable1_CellTypeSpecificGenes_Master3Sorted.csv', header = TRUE)#Reads in profiles from Zhang et al data
profiles=profiles[-which(as.character(profiles$CellType_Primary)=="Mural"),]


# Match the genes in the Hagenauer list to the genes in the Stanley Data set
tmp=match(toupper(gene.symbols),toupper(profiles$Gene.Symbol..Human.)) 

rel.profiles <- profiles[tmp[!is.na(tmp)],]
rel.obsv.exprs <- obsv.exprs[!is.na(tmp),] 
gene.names=gene.symbols[!is.na(tmp)] 
######################################################################
cell.type=as.factor(as.character(rel.profiles$Umbrella.Cell.Type))
gene.means <- apply(rel.obsv.exprs,1,mean, na.rm = TRUE) 
gene.medians <- apply(rel.obsv.exprs,1,median, na.rm = TRUE)  

# How centered are these genes? 
plot(gene.means,gene.medians)

cell.type.genes=c()
for(ii in 1:length(unique(cell.type))){
  
cell.gene.avg=gene.means[which(cell.type==unique(cell.type)[ii])]
assign(paste0(unique(cell.type)[ii],".means"),cell.gene.avg,envir=.GlobalEnv)
cell.gene.med=gene.medians[which(cell.type==unique(cell.type)[ii])]
assign(paste0(unique(cell.type)[ii],".means"),cell.gene.avg,envir=.GlobalEnv)
cell.genes=gene.names[which(cell.type==unique(cell.type)[ii])]

xx=which((cell.gene.avg/cell.gene.med)>1.1)
hist(cell.gene.avg[-xx],breaks=40)

yy=which(cell.gene.avg<as.numeric(summary(cell.gene.avg)[2]))
hist(cell.gene.avg[-yy],breaks=40)

#zz=which(cell.gene.avg>min(Outlier(cell.gene.avg, na.rm=TRUE)))
zz=which(cell.gene.avg>as.numeric(summary(cell.gene.avg)[5]))
hist(cell.gene.avg[-zz],breaks=40)

rid=unique(c(xx,yy,zz))
assign(paste0(unique(cell.type)[ii],".genes"),cell.genes[-rid],envir=.GlobalEnv)
cell.type.genes=c(cell.type.genes,cell.genes[-rid])

}
tmp=match(toupper(gene.symbols),toupper(cell.type.genes))

rel.profiles <- profiles[tmp[!is.na(tmp)],]
rel.obsv.exprs <- obsv.exprs[!is.na(tmp),] 
gene.names=gene.symbols[!is.na(tmp)] 

tmp=match(toupper(gene.names),toupper(cell.type.genes))
cell.type.edit=as.factor(as.character(rel.profiles$Umbrella.Cell.Type))[tmp]
  
out.pts=matrix(ncol=2,nrow=0)
gene.skew=1:nrow(rel.obsv.exprs)
gene.kurt=gene.skew=1:nrow(rel.obsv.exprs)
for(ii in 1:nrow(rel.obsv.exprs)){
  gene.skew[ii]=skewness(rel.obsv.exprs[ii,])
  gene.kurt[ii]=kurtosis(rel.obsv.exprs[ii,],type=3)
  xx=match(names(Outlier(rel.obsv.exprs[ii,])),colnames(rel.obsv.exprs))
  xx=cbind(rep(ii,length(xx)),xx)
  out.pts=rbind(out.pts,xx)
}
# This prints out a table of the the number of outliers in each gene of all the genes that have at least one outlier
out.pr.gene=table(gene.names[out.pts[,1]])
rid2=unique(c(which(abs(gene.skew)>1),which(gene.kurt>3),which(out.pr.gene>1)))

table(cell.type.edit[-rid2])

rel.profiles.edit = rel.profiles[-rid2,]
rel.obsv.exprs.edit = rel.obsv.exprs[-rid2,] 
gene.names.edit=gene.names[-rid2] 
cell.type.edit=cell.type.edit[-rid2]


# Calculate the average expression value for each cell type geneset for each sample
samp.gene.means=matrix(ncol=5,nrow=ncol(rel.obsv.exprs.edit)) 
for(ii in 1:ncol(rel.obsv.exprs.edit)){
  samp.gene.means[ii,] <- tapply(rel.obsv.exprs.edit[,ii],cell.type.edit,mean) 
}

# Plot the average cell type geneset expression values by cell type
setwd(output)
jpeg("Boxplot of sample means for cell spec genesets.jpeg")
boxplot(samp.gene.means,names=sort(unique(cell.type)),las=2,cex.names=.7)
dev.off()

# For each cell type attempt to find "housekeeping genes" based on their distributions
# This means that they would be genes' whose expression is fairly normally distributed
# They also would not be too skewed or kurtotic

jpeg("Outliers Per Cell Specific Gene List Histogram Stanley-CTL based.jpeg")
hist(out.pr.gene)
dev.off()
summary(out.pr.gene)

jpeg("Skewness of Cell Specific Gene List Histogram Stanley-CTL based.jpeg")
hist(gene.skew)
dev.off()
summary(gene.skew)

jpeg("Kurtosis of Cell Specific Gene List Histogram Stanley-CTL based.jpeg")
hist(gene.kurt)
dev.off()
summary(gene.kurt)

# Let's rid ourselves of the genes whose data 1) has any outliers, 2) has a skew, and 3) Is too platykurtotic, and maybe even ones that are too leptokurtotic

xx=which(out.pr.gene>2)
xx=match(names(xx),gene.symbols)
xx=xx[!is.na(xx)]
rel.obsv.exprs.edit=rel.obsv.exprs[-unique(c(xx,which(abs(gene.skew)>.7),which(abs(gene.kurt)>3))),]
cell.type.edit=cell.type[-unique(c(xx,which(abs(gene.skew)>.7),which(abs(gene.kurt)>3)))]
cell.type.genes=gene.symbols[-unique(xx,which(abs(gene.skew)>.7),which(abs(gene.kurt)>3))]

# Let's recalculate the average expression value for each cell type geneset for each sample
# after having removed the funky distribution genes
gene.means.edit <- apply(rel.obsv.exprs.edit,1,mean, na.rm = TRUE)  
gene.medians.edit <- apply(rel.obsv.exprs.edit,1,median, na.rm = TRUE) 
mean.med.ratio=gene.means.edit/gene.medians.edit

xx=unique(ratio.keep,mean.keep,median.keep)
rel.obsv.exprs.edit=rel.obsv.exprs.edit[xx,]
cell.type.edit=cell.type.edit[xx]
cell.type.genes=cell.type.genes[xx]

samp.gene.means.edit=matrix(ncol=5,nrow=ncol(rel.obsv.exprs.edit))
for(ii in 1:ncol(rel.obsv.exprs.edit)){
  samp.gene.means.edit[ii,] <- tapply(rel.obsv.exprs.edit[,ii],cell.type.edit,mean) 
}


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

