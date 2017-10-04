
####### Script for One way ANOVA and Tukey's ###########

data="put a matrix with your data here that has the rows=probes/genes and columns=samples "
cutoff="set your abundance cutoff"
group="vector with the groups of the different conditions. If you have multiple, put group1*group2"
csvfile="put the name you want the csv file to be called.csv"
pvalcutoff=.05



### Import your gene expression data matrix with dimnames of genes and sample names set
read.table("")
data=data[which(apply(data,1,mean)>cutoff),]
data=t(data)
data=as.data.frame(data)

## ANOVA p values
data.anova=matrix(nrow=ncol(data),ncol=(length(unique(group)))+3)
for (ii in 1:ncol(data)){
  data.anova[ii,1:2]= summary(aov(data[,ii]~group,data=data))[[1]][["Pr(>F)"]]
}

## FDR
data.anova[,2]=p.adjust(data.anova[,1],method ="hochberg")

## Get rid of p values below cutoff
signif=data[,which(data.anova[,2]<.pvalcutoff)]

### Take average by group
signif.avg=matrix(ncol=ncol(signif),nrow=length(unique(group)))
for (ii in 1:ncol(signif)){
  signif.avg[,ii]=tapply(signif[,ii],group,mean)
}
colnames(signif.avg)=dimnames(signif)[[2]]
rownames(signif.avg)=sort(unique(group))

##Post hoc of Tukey's HSD
tukey.pval=matrix(nrow=ncol(signif),ncol=length(TukeyHSD(x=aov(signif[,1]~group,data=signif), c('group'), conf.level=0.95)$group[,1]))
for (ii in 1:ncol(signif)){
  tmp=TukeyHSD(x=aov(signif[,ii]~group,data=signif), c('group'), conf.level=0.95)
  tukey.pval[ii,]=t(tmp$group[,4])
}
colnames(tukey.pval)=dimnames(tmp$group)[[1]]
rownames(tukey.pval)=dimnames(signif)[[2]]


#### Calculate fold changes for each combination of cohorts
fold.change=matrix(ncol=nrow(signif.avg)*(nrow(signif.avg)-1),nrow=ncol(signif.avg))
dimnames(fold.change)[[2]]=rep(NA,nrow(signif.avg)*(nrow(signif.avg)-1))
dimnames(fold.change)[[1]]=dimnames(signif)[[2]]

for(ii in 1:nrow(signif.avg)){
  xx=rep(c(1:nrow(signif.avg)),nrow(signif.avg))[c(ii+1:(nrow(signif.avg)-1))]
    fold.change[,(((ii-1)*5)+1):(ii*5)]=t(signif.avg[ii,]/signif.avg[c(xx),])
    dimnames(fold.change)[[2]][(((ii-1)*5)+1):(ii*5)]=c(paste(sort(unique(group))[ii],"/",sort(unique(group))[xx]))
}

##### Put it all together in a CSV file
for(ii in 1:ncol(tukey.pval)){
  xx=order(tukey.pval[,ii])
  combo1=paste(strsplit(colnames(tukey.pval)[ii],"-")[[1]][1],"/",strsplit(colnames(tukey.pval)[ii],"-")[[1]][2])
  combo2=paste(strsplit(colnames(tukey.pval)[ii],"-")[[1]][2],"/",strsplit(colnames(tukey.pval)[ii],"-")[[1]][1])
  write.csv(cbind(tukey.pval[xx,ii],t(signif.avg)[xx,],fold.change[xx,which(dimnames(fold.change)[[2]]==combo1|dimnames(fold.change)[[2]]==combo2)]),file=paste(colnames(tukey.pval)[ii],csvfile),quote=FALSE)
}



