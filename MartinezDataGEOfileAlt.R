#### Alternate Differential expression tests. 
name="MartinezDataGEOfileAlt"
output="/Users/cansav091/Desktop/Neurogenomics Data/Martinez M0M1M2 Lists/MartinezDataGEOfileAlt"
############ Directory Commands ##################
# Type this before creating any output files
setwd(output)

### Need these packages: 
library(Biobase)
library(GEOquery)
library(psych)
library(CondReg)
library(Biostrings)
library(calibrate)
library(MASS)

######### Run the Martinez Data GEO file code first ###############

geo.groups.alt=geo.groups
geo.data.alt=geo.data


colnames(geo.data.alt)=geo.groups.alt
####### Get rid of genes/probes that aren't very abundant ##########
####### Get rid of bottom quartile genes ###########################
# (renamed gene.avg to be called probe.avg)
probe.avg=apply(geo.data.alt,1,mean)
cutoff=.10
cutoff=sort(probe.avg)[nrow(geo.data)*cutoff]

geo.abund=geo.data.alt[which(probe.avg> cutoff),]
gene.abund.alt=gene.name[which(probe.avg> cutoff)]


#### If there are multiple probes per gene, use Hotelling's to find DE genes #########
if(length(unique(gene.name))!=length(gene.name)){
  ######### Hotelling's MANOVA #####################################
  fstat=c()#### A bunch of empty vectors to store stats in. 
  hotelling.trace=c()
  pval=c()
  num.probes=c() # number of probes per gene
  df1=c()
  df2=c()
  ctl.sd=c()
  for(ii in 1:length(unique(gene.abund.alt))){ 
    xx=which(gene.abund.alt==unique(gene.abund.alt)[ii])  ### Finds row number index for all the probes for a particular gene
    if(length(xx)<2){
      ctl.sd[ii]=sd(geo.abund[xx,which(geo.groups=="Control")])
    }else{
      ctl.sd[ii]=sd(apply(geo.abund[xx,],2,mean)[which(geo.groups=="Control")])
    }
    num.probes[ii]=p=length(xx) ### Takes note of how many probes are in each set. Sets "p" as the dimensions
    xx=as.matrix(geo.abund[xx,]) ## Takes the data from the row indices indicating previously indicated 
    
    if(p<2){ ## When there is only 1 probe for a gene, This set of code does normal ANOVA and by passes the rest of the code to the next loop
      n=length(xx) 
      xx=summary(aov(xx~as.factor(geo.groups.alt),data=as.data.frame(xx)))## Summarizes aov object into the report stats.
      fstat[ii]=xx[[1]][["F value"]][1] ### Stores F statistic in the same vector as genes with multiple probes
      pval[ii]=xx[[1]][["Pr(>F)"]][1] ## Stores p value for 1 probe sets with the other pvalues
      df1[ii]=NA 
      df2[ii]=NA 
      hotelling.trace[ii]=NA
      next
    }else{ ## For genes with more than one probe: 
      n=ncol(xx) 
      xx=manova(t(xx)~as.factor(geo.groups.alt))
      ss.resid =crossprod(xx$residuals) ## Finds the crossproduct of residuals matrix
      ss.effect=crossprod(xx$effects[-1,]) ## Finds the crossproduct of effect matrix, leaves out the intercept, keeps only the effects for each group
      
      if(1/kappa(ss.resid)<1e-10) {## If estimated condition number for the inverse matrix is too small (which is what happens to all the p>>n scenarios; When there are so many probes per gene)
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
      pval[ii]=1-pf(fstat[ii], df1[ii], df2[ii]) ## calculates the p value based on the f distribution and df parameters
    }
  }
  length(which(pval<.05)) #Before FDR correction
  pval.fdr=p.adjust(pval,method="hochberg") 
  length(which(pval.fdr<.05))
}
#################################################################
######Put basic info in a summary file ##########################
basic.info=as.vector(c(geo.soft@header$geo_accession,geo.soft@header$summary[1],"hgu133A and hgu133b",table(geo.groups.alt),length(unique(geo.groups.alt)),summary(probe.avg)[1:6],nrow(geo.data),sum(!is.na(gene.name)),length(unique(gene.name)),length(which(pval<.05)),length(which(pval.fdr<.05))))
names(basic.info)=c("GSE","Experiment Description:","Data/Microarray Type:",sort(unique(geo.groups.alt)),"Number of Total Samples:","Min:","1st Qu:","Median:","Mean:", "3rd Qu:","Max:","Number of probes:","Number of Probes with Gene Names:","Number of Genes:","Number of p<.05 Genes before FDR","Number of p<.05 Genes FDR")
setwd(outDir)
write.table(basic.info,file=paste(data.name,"Basic Info File Alt.txt"),quote=FALSE)

############################################################################
############## Find average for each probe/gene by groups###################
geo.avg=matrix(nrow=nrow(geo.abund),ncol=length(unique(geo.groups.alt)))
for (ii in 1:length(geo.abund[,1])){
  tmp=tapply(geo.abund[ii,],geo.groups.alt,FUN=mean)
  geo.avg[ii,]=tmp
}

############################################################################
########### Calculate fold changes across groups ###########################
dimnames(geo.avg)[[2]]=unique(geo.groups.alt)
dimnames(geo.avg)[[1]]=gene.abund.alt

nn=length(unique(geo.groups.alt))
columns.names=c()
comp=matrix(nrow=2,ncol=0)
for(ii in 1:nn){
  tmp=matrix(ncol=nn,nrow=length(gene.abund.alt))
  for(jj in 1:nn){
    tmp[,jj]=geo.avg[,ii]/geo.avg[,jj]
    columns.names=c(columns.names,paste0(unique(geo.groups.alt)[ii],"\\",unique(geo.groups.alt)[jj]))
    comp=cbind(comp,rbind(ii,jj))
  }
  if(ii==1){
    geo.fc=tmp
  }else{
   geo.fc=cbind(geo.fc,tmp)
  }
}
geo.fc=geo.fc[,-which(comp[1,]==comp[2,])]
colnames(geo.fc)=columns.names[-which(comp[1,]==comp[2,])]
comp=comp[,-which(comp[1,]==comp[2,])]
rownames(geo.fc)=gene.abund.alt

########### Put fold changes into a .csv file ##############################
setwd(outDir)
write.csv(geo.fc,paste(data.name,"Fold Changes Alt.csv"),quote=FALSE)

#### Descriptive statistics for the probesets for each gene ###########
if(exists("hotelling.trace")==TRUE){
  sd.probes=c(1:length(unique(gene.abund.alt)))
  avg.group.signal=matrix(nrow=length(unique(gene.abund.alt)),ncol=length(unique(geo.groups.alt)))
  for(ii in 1:length(unique(gene.abund.alt))){
    xx=which(gene.abund.alt==unique(gene.abund.alt)[ii])
    if(length(xx)<2){
      avg.group.signal[ii,]=geo.avg[xx,]
      sd.probes[ii]=NA
      next
    }else{
      avg.group.signal[ii,]=apply(geo.avg[xx,],2,mean)
      sd.probes[ii]=sd(apply(geo.abund[xx,],1,mean))
    }
  }
}

  signif=cbind(hotelling.trace,fstat,pval,pval.fdr,num.probes,sd.probes,avg.group.signal)
  rownames(signif)=unique(gene.abund.alt)
  colnames(signif)=c("Hotelling Trace","F-stat","p value","FDR adj p val","Num of Probesets","SD of probsets",unique(geo.groups.alt))


########Put the significance info into a csv file #######################
setwd(outDir)
write.csv(signif[order(signif[,4]),],file=paste(data.name,"Genes Alt.csv"))



#################### Plots #############################################

setwd(output)
jpeg("FC vs ANOVA -log10pval Martinez Data.jpeg",width=1500,height=1500)
xx=cor(abs(geo.fc[]),-log10(pval))
plot(abs(overall.fc),-log10(pval),pch="",xlab="Absolute Value for High Dominance/Low Dominance FC",cex.lab=2 )
textxy(abs(overall.fc),-log10(pval),gene.abund,offset=0,cex=1.3)
abline(h=1.30,col="red")
abline(v=1.30,col="blue")
dev.off()
