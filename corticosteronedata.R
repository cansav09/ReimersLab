## urigen et al data sorting 
data.name="Urigen Corticosterone"
urigen.data=as.matrix(read.table("GSE9797_series_matrix.txt",sep="\t",row.names=1))

urigen.affy.ids=rownames(urigen.data)[-1]
sample.IDs=urigen.data[1,]
urigen.data=urigen.data[2:nrow(urigen.data),]
urigen.data=apply(urigen.data,2,as.numeric)

### Matching Annotation for ST 1.0 array to dataset
library(rat2302.db)

xx=rat2302SYMBOL
probeids <- mappedkeys(xx)
yy=as.data.frame(xx[probeids])
gene.name=yy$symbol[match(urigen.affy.ids,yy$probe_id)] #This tells you in order of the queried entries,what the coordinates are of the second group of entries only  the first match is given 


### Get rid of the data that don't have gene names that match to them
urigen.data=urigen.data[!is.na(gene.name),]
sum(!is.na(gene.name)) #N=20995 out of 31100 probes have gene names
gene.name=gene.name[!is.na(gene.name)]


## Get rid of genes that aren't very abundant. 
probe.avg=apply(urigen.data,1,mean)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.7763  0.9946  1.0050  1.0090  1.0170  1.8190 

urigen.abund=urigen.data[which(probe.avg> summary(probe.avg)[2]),]
urigen.gene.names.abund=gene.name[which(probe.avg> summary(probe.avg)[2])]

## Get sample info from SOFT files
library(Biobase)
library(GEOquery)
urigen <- getGEO(filename='GSE9797_family.soft')
urigen.sample.info=matrix(nrow=length(urigen.data[1,]),ncol=2)
for (ii in 1:length(urigen.data[1,])){
  urigen.sample.info[ii,1:2]=c(urigen@header$sample_id[ii],urigen@gsms[[ii]]@header$title,)
}
urigen.groups=as.vector(nrow(urigen.sample.info))
urigen.groups[grep("Control",urigen.sample.info[,2])]="Control"
urigen.groups[grep("corticosterone",urigen.sample.info[,2])]="Corticosterone"



############ quick ttest to check for differences between groups
ttest.pval=c()
for(ii in 1:nrow(urigen.abund)){
ttest.pval[ii]=t.test(urigen.abund[ii,which(urigen.groups=="Control")],urigen.abund[ii,which(urigen.groups=="Corticosterone")])$p.value
}
hist(-log10(ttest.pval))
abline(v=1.30,col="red")
dev.off()

library(psych)
######### Hotelling MANOVA #####################################
fstat=c()#### A bunch of empty vectors to store stats in. 
hotelling.trace=c()
hotelling.pval.test=c()
num.probes=c() # number of probes per gene
df1=c()
df2=c()
for(ii in 1:length(unique(urigen.gene.names.abund))){ #For each gene name, 
  xx=which(urigen.gene.names.abund==unique(urigen.gene.names.abund)[ii]) 
  ### Finds row number index for all the probes for a particular gene
  num.probes[ii]=p=length(xx) ### Takes note of how many probes are in each set. Sets "p" as the dimensions
  xx=as.matrix(urigen.abund[xx,])#Takes the data from the row indices indicating previously indicated 
  n=ncol(xx) ## number of samples which is always 12 for this dataset.
  if(p<2){ #When there is only 1 probe for a gene, This set of code does normal ANOVA and by passes the rest of the code to the next loop
    xx=summary(aov(xx~as.factor(urigen.groups),data=as.data.frame(xx)))## summarizes aov object into the report stats.
    # note: The only reason this is set to groups  is that we removed the necrotic groups
    fstat[ii]=xx[[1]][["F value"]][1] ### Stores F statistic in the same vector as genes with multiple probes
    hotelling.pval.test[ii]=xx[[1]][["Pr(>F)"]][1] ## Stores p value for 1 probe sets with the other pvalues
    df1[ii]=NA # These are just place holders because the last ten genes only have one probe, and we need to still have the same number of rows
    df2[ii]=NA 
    hotelling.trace[ii]=NA
    next
  }else{ # For genes with more than one probe: 
    xx=manova(t(xx)~as.factor(urigen.groups))
    ss.resid =crossprod(xx$residuals) ## Finds the crossproduct of residuals matrix
    ss.effect=crossprod(xx$effects[-1,]) ## Finds the crossproduct of effect matrix, leaves out the intercept, keeps only the effects for each group
    cc=try(solve(ss.resid), silent=T) # This tests whether or not it can find the inverse of the residual matrix
    if(is(cc,"try-error")) {## If the inverse of the matrix cannot be found (which is what happens to all the p>>n scenarios; 
      # when there are so many probes per gene)
      A1=ss.effect%*%condreg(ss.resid,8)$invS # Those gene sets regularized by the function condreg from Won et al, 2013
    }else{ # if the inverse of the matrix can be found, the variance regularization step is skipped
      A1=ss.effect%*%solve(ss.resid) ### 
    }
    m=length(xx$assign) ### This comes out to number of groups +1 for intercept
    u=(n-m-p-1)/2 ### Here is the traditional u is calculated. 
    u=ifelse(u<0,0,u) ## However, when p>>n u becomes negative, which throws off the df calculation,
    #so here if it turns out that u is negative, we use 0 instead. (I don't know if this is a good way to do it.)
    t=(abs(p-m+1)-1)/2 ## t gets abnormally large when p>>n, don't know if this is something that needs to be accounted for. 
    s=min(c(p,m-1))### will usually be 3, however, for probesets where there are only two probes, then s=2
    df1[ii]=s*(2*t+s+1)### will be between 6 and 24 for "normal" p<n data. Our data however has a mean df1 of 28 and a max of 321. 
    ###Don't know how these super large df1's should be accounted for. 
    df2[ii]=2*(s*u+1) #### will be between 2 and 14 for our data. 
    ### with a really large p, this "should" become negative, however, that doesn't make sense as a df.
    ## So here I've limited it to be no lower than 2. But, again, don't know how this should really be dealt with. 
    hotelling.trace[ii]=tr(A1) ### Hotelling stat is based on the trace of the matrix found from ss.effects*ss.residuals^-1
    fstat[ii]=(tr(A1)/s)*(df2[ii]/df1[ii])## Adjust F statistic to be proportional to the degrees of freedom for effects and residuals
    hotelling.pval.test[ii]=1-pf(fstat[ii], df1[ii], df2[ii]) ## calculates the p value based on the f distribution and df parameters
  }
}


length(which(hotelling.pval.test<.05)) #Before FDR correction, N=9801 genes are found significant
pval.fdr=p.adjust(hotelling.pval.test,method="hochberg") ##After FDR, 15 significant genes p<.05
length(which(pval.fdr<.05))

################ AVERAGE by groups#####################
urigen.avg=matrix(nrow=nrow(urigen.abund),ncol=length(unique(urigen.groups)))
for (ii in 1:length(urigen.abund[,1])){
  tmp=tapply(urigen.abund[ii,],urigen.groups,FUN=mean)
  urigen.avg[ii,]=tmp
}

##################### Fold Change ##########################################
dimnames(urigen.avg)[[2]]=unique(urigen.groups)
dimnames(urigen.avg)[[1]]=urigen.gene.names.abund

urigen.fc=matrix(ncol=5,nrow=nrow(gene.names.abund))
for(ii in 2:length(unique(urigen.groups))){
  urigen.fc[,ii-1]=urigen.avg[,ii]/urigen.avg[,1]
}

urigen.fc=urigen.avg[,1]/urigen.avg[,2]
attributes(urigen.fc)[[1]]=dimnames(urigen.avg)[[1]]

hist(log2(urigen.fc))

####################### Other descriptive statistics for the probesets for each gene#########################
sd.probes=c()
avg.group.signal=matrix(nrow=length(unique(urigen.gene.names.abund)),ncol=2)
for(ii in 1:length(unique(urigen.gene.names.abund))){
  xx=which(urigen.gene.names.abund==unique(urigen.gene.names.abund)[ii])
  if(length(xx)<2){
    avg.group.signal[ii,]=urigen.avg[xx,]
    sd.probes[ii]=NA
    next
  }else{
    avg.group.signal[ii,]=apply(urigen.avg[xx,],2,mean)
    sd.probes[ii]=sd(apply(urigen.abund[xx,],1,mean))
  }
}

signif=cbind(hotelling.trace,fstat,hotelling.pval.test,pval.fdr,num.probes,sd.probes,avg.group.signal)
rownames(signif)=unique(urigen.gene.names.abund)
colnames(signif)=c("Hotelling Trace","F-stat","p value","FDR adj p val","Num of Probesets","SD of probsets",unique(urigen.groups))

write.csv(signif[order(signif[,4]),],file=paste(data.name,"Genes.csv"))
plot(num.probes,-log10(pval.fdr))

