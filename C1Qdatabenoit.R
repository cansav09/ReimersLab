## Benoit et al data sorting 

benoit.data=as.matrix(read.table("C1qMacrophageGSE30177.txt",row.names=1))

benoit.genes=rownames(benoit.data)[2:253319]
sample.IDs=benoit.data[1,]
benoit.data=benoit.data[2:253319,]
benoit.data=apply(benoit.data,2,as.numeric)

library(Biobase)
library(GEOquery)
library(hugene10stprobeset.db)

### Matching Annotation for ST 1.0 array to dataset

xx=hugene10stprobesetSYMBOL
probeids <- mappedkeys(xx)
yy=as.data.frame(xx[probeids])
gene.name.benoit=yy$symbol[match(benoit.genes,yy$probe_id)] #This tells you in order of the queried entries,what the coordinates are of the second group of entries only  the first match is given 

### Get rid of the data that don't have gene names that match to them
benoit.data=benoit.data[!is.na(gene.name.benoit),]
sum(!is.na(gene.name.benoit))
## N= 232492 probes with gene names to them.
gene.name.benoit=gene.name.benoit[!is.na(gene.name.benoit)]

## Get rid of genes that aren't very abundant. 
gene.avg=apply(benoit.data,1,mean)

# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -4.4090 -1.2030  0.7090  0.6993  2.4970  8.0480 

benoit.abund=benoit.data[which(gene.avg> 0),]
benoit.gene.names.abund=gene.name.benoit[which(gene.avg> 0)]


# This leaves N= 196038 probes that are more abundant than an average signal of -2

## Get sample info from SOFT files
Benoit <- getGEO(filename='GSE30177_family.soft')
benoit.sample.info=matrix(nrow=length(benoit.data[1,]),ncol=2)
for (ii in 1:length(benoit.data[1,])){
  benoit.sample.info[ii,1:2]=c(Benoit@header$sample_id[ii],Benoit@gsms[[ii]]@header$title)
}
#When combined with propidium iodide (PI) the double labeling procedure
# necrotic (annexin V+/PI+) cells = LAL
# apoptotic (annexin V+/PI-) cells = EAL

benoit.groups= c(rep("control",3),rep("LPS",3),rep("Apoptotic LPS",3),rep("C1q Apoptotic LPS",3),rep("Necrotic LPS",3),rep("C1q necrotic LPS",3))

###### MANOVA and Hotelling-Lawley trace with flexibility for p>>n cases ##############################
non.necrotic.groups=signif[,1:12]

fstat=c()
hotelling.trace=c()
hotelling.pval.test=c()
for(ii in 1:length(unique(benoit.gene.names.abund))){
  xx=t(as.matrix(non.necrotic.groups[which(benoit.gene.names.abund==unique(benoit.gene.names.abund[2])),]))
  n=ncol(xx)## Number of observations
  p=nrow(xx)
  if(p<3){
    next
  }else{
    xx=manova(xx~benoit.groups[1:12])
    ss.resid =crossprod(xx$residuals)
    ss.effect=crossprod(xx$effects[2:4,])
    if(p>8){
      A1=ss.effect%*%condreg(ss.resid,8)$invS
      u=.5
      t=1
    }else{
      m=4 ### Number of groups +1 for intercept
      A1=ss.effect%*%solve(ss.resid)
      t=(abs(p-m+1)-1)/2### will be between 1 and -.5 t increases with more dimensions, p. When p=4, t=0
      u=(n-m-p-1)/2 ##will be between .5 and 3 as the dimensions increase, u becomes more negative. When p=7, u=0, if p>7, u is negative
    }### four groups
    s=min(c(p,m-1))### s will always be 2. 
    df1=s*(2*t+s+1)### will be between 9 and 18 This will be really big with a really big p
    df2=2*(s*u+1) #### will be between 5 and 20 This will become negative and really small with a large p
    hotelling.trace[ii]=tr(A1)
    fstat[ii]=(tr(A1)/s)*(df2/df1)
    hotelling.pval.test[ii]=1-pf(fstat[ii], df1, df2)
  }
}

################ Looking at correlations between probes of the same genes################
library(calibrate)
library(MASS)
probe.per.gene=table(as.factor(gene.name.benoit))
sd(gene.avg[which(dimnames(probe.per.gene)[[1]][1]==gene.name.benoit)])

probeset.corr=matrix(nrow=length(dimnames(probe.per.gene)[[1]]),ncol=4)
probe.set.sd=c()
loadings.per.probeset=matrix(ncol=3)
sample.gene.composite=matrix(nrow=length(dimnames(probe.per.gene)[[1]]),ncol=18)
for (ii in 1: length(unique(dimnames(probe.per.gene)[[1]]))){
  if(probe.per.gene[ii]<3){
    next
  }else{ 
  probe.set.sd[ii]=sd(gene.avg[which(dimnames(probe.per.gene)[[1]][ii]==gene.name.benoit)])
  xx=benoit.data[which(dimnames(probe.per.gene)[[1]][ii]==gene.name.benoit),]
  yy=as.vector(cor(t(xx)))[which(as.vector(cor(t(xx)))<1)]
  probeset.corr[ii,]=c(median(yy),sd(yy),range(yy))
  xx=prcomp(t(xx))
  loadings.per.probeset=rbind(loadings.per.probeset,cbind(gene.name.benoit[which(dimnames(probe.per.gene)[[1]][ii]==gene.name.benoit)],probeids[which(dimnames(probe.per.gene)[[1]][ii]==gene.name.benoit)],xx$rotation[,1]))
  jpeg(paste(dimnames(probe.per.gene)[[1]][ii],"PCA plot of probeset per gene.jpeg "))
  plot(xx$rotation)
  textxy(xx$rotation[,1],xx$rotation[,2],probeids[which(dimnames(probe.per.gene)[[1]][ii]==gene.name.benoit)],cex=.75,offset=.6)
  dev.off()
  }
}

jpeg(paste(dimnames(probe.per.gene)[[1]][ii]," PCA plot of probeset per gene.jpeg "))
plot(xx$rotation)
textxy(xx$rotation[,1],xx$rotation[,2],probeids[which(dimnames(probe.per.gene)[[1]][ii]==gene.name.benoit)],cex=.75,offset=.6)
dev.off()

################ AVERAGE by groups#####################
benoit.avg=matrix(nrow=length(benoit.abund[,1]),ncol=6)
for (ii in 1:length(benoit.abund[,1])){
  tmp=tapply(benoit.abund[ii,],benoit.groups,FUN=mean)
  benoit.avg[ii,]=tmp
}


##################### Fold Change ##########################################

dimnames(benoit.avg)[[2]]=unique(benoit.groups)
dimnames(benoit.avg)[[1]]=benoit.gene.names.abund

benoit.fc=matrix(ncol=5,nrow=196038)
for(ii in 2:6){
  benoit.fc[,ii-1]=benoit.avg[,ii]/benoit.avg[,1]
}
dimnames(benoit.fc)[[1]]=dimnames(benoit.avg)[[1]]

hist(log2(benoit.fc[,1]))

head(cbind(dimnames(benoit.fc)[[1]][order(abs(benoit.fc[,1]),decreasing=TRUE)],benoit.fc[order(abs(benoit.fc[,1]),decreasing=TRUE)]),20)
head(cbind(dimnames(benoit.fc)[[1]][order(abs(benoit.fc[,1]),decreasing=TRUE)],benoit.fc[order(abs(benoit.fc[,1]),decreasing=TRUE)]),20)
head(cbind(dimnames(benoit.fc)[[1]][order(abs(benoit.fc[,1]),decreasing=TRUE)],benoit.fc[order(abs(benoit.fc[,1]),decreasing=TRUE)]),20)
head(cbind(dimnames(benoit.fc)[[1]][order(abs(benoit.fc[,1]),decreasing=TRUE)],benoit.fc[order(abs(benoit.fc[,1]),decreasing=TRUE)]),20)



############### HCL of samples #####################
library(dplyr)
library(ggplot2)
library(plyr)
library(dendextend)

color.group=mapvalues(benoit.groups,from = c(unique(benoit.groups)), to = c("blue","red","brown","green","black","purple"))

jpeg("Cluster Analysis of Samples w IDs.jpeg",width=700,height=500)
par(mar=c(10,3,3,3))
dend <- as.data.frame(t(benoit.abund)) %>%  scale %>% dist %>% hclust %>% as.dendrogram
dend %>% plot
labels(benoit.groups)
dend %>% set("labels_col", color.group[order(unlist(dend))]) %>% set("labels_cex", 1)  %>% plot(main = "Cluster of Cell Types/Treatments")

legend(900,c(unique(benoit.groups)),lty=c(1,1), lwd=c(2.5,2.5),col=c(unique(color.group)),cex=.75) 
dev.off()


######################## ANOVA ###############################################################
benoit.abund=t(benoit.abund)
benoit.abund=as.data.frame(benoit.abund)
colnames(benoit.abund)=benoit.gene.names.abund

benoit.data.anova=matrix(nrow=ncol(benoit.abund),ncol=2)
for (ii in 1:ncol(benoit.abund)){
 benoit.data.anova[ii,]= summary(aov(benoit.abund[,ii]~benoit.groups,data=benoit.abund))[[1]][["Pr(>F)"]]
}


### Adjust the p values
benoit.data.anova[,2]=p.adjust(benoit.data.anova[,1],method = "hochberg")
signif=benoit.abund[which(benoit.data.anova[,1]<.05),]
signif.avg=benoit.avg[which(benoit.data.anova[,2]<.05),]

length(unique(gsub(".[0-9]{1,2}$","",dimnames(signif)[[2]])))## 280 genes are significant after p value correction 
probe.per.signif.gene=table(as.factor(gsub("\\.[0-9]{1,2}$","",dimnames(signif)[[2]])))

hist(probe.per.signif.gene/probe.per.gene[match(dimnames(probe.per.signif.gene)[[1]],dimnames(probe.per.gene)[[1]])],main="Proportion of Probes that Show Up as Signif",xlab="Ratio")
xx=probe.per.signif.gene/probe.per.gene[match(dimnames(probe.per.signif.gene)[[1]],dimnames(probe.per.gene)[[1]])]

plot(as.numeric(probe.per.signif.gene),as.numeric(xx),pch=20,ylab="Proportion of probes that showed as signif",xlab="Number of probes total")
textxy(as.numeric(probe.per.signif.gene),as.numeric(xx),dimnames(xx)[[1]])

# 878 probes are significant after p value correction 

### post hoc on probes that are significant 
benoit.tukey.diff=matrix(ncol=15,nrow=ncol(signif))
benoit.tukey.pval=matrix(ncol=15,nrow=ncol(signif))
for (ii in 1:ncol(signif)){
  xx=aov(signif[,ii]~benoit.groups,data=benoit.abund)
  tmp=TukeyHSD(x=xx, 'benoit.groups', conf.level=0.95)
  benoit.tukey.pval[ii,]=t(tmp$benoit.groups[,4])
  benoit.tukey.diff[ii,]=t(tmp$benoit.groups[,1])
}
colnames(benoit.tukey.pval)=dimnames(tmp$benoit.groups)[[1]]
colnames(benoit.tukey.diff)=dimnames(tmp$benoit.groups)[[1]]

rownames(benoit.tukey.pval)=dimnames(signif)[[2]]
rownames(benoit.tukey.diff)=dimnames(signif)[[2]]

###########
#[1] "C1q Apoptotic LPS-Apoptotic LPS"    
#[2] "C1q necrotic LPS-Apoptotic LPS"    
#[3] "control-Apoptotic LPS"              
#[4] "LPS-Apoptotic LPS"                 
#[5] "Necrotic LPS-Apoptotic LPS"        
#[6] "C1q necrotic LPS-C1q Apoptotic LPS"
#[7] "control-C1q Apoptotic LPS"          
#[8] "LPS-C1q Apoptotic LPS"             
#[9] "Necrotic LPS-C1q Apoptotic LPS"     
#[10] "control-C1q necrotic LPS"          
#[11] "LPS-C1q necrotic LPS"             
#[12] "Necrotic LPS-C1q necrotic LPS"     
#[13] "LPS-control"                       
#[14] "Necrotic LPS-control"              
#[15] "Necrotic LPS-LPS"   

grep("control",dimnames(tmp$benoit.groups)[[1]])


######### Sorting the pvals and output #################################
colnames(signif.avg)
length(unique(rownames(signif.avg)))

## QQ plot

apop.v.control.fc=signif.avg[,3]/signif.avg[,1]
c1qapop.v.control.fc=signif.avg[,4]/signif.avg[,1]
necro.v.control.fc=signif.avg[,5]/signif.avg[,1]
c1qnecro.v.control.fc=signif.avg[,6]/signif.avg[,1]
lps.control.fc=signif.avg[,2]/signif.avg[,1]

index=c(3,7,14,10,13)
index2=c(3,4,5,6,2)
name=c("Apoptotic","C1q Apoptotic","Necrotic","C1q Necrotic","LPS")

### Write into a csv file and make QQ plots

for (ii in 1:5){
  xx=order(benoit.tukey.pval[,index[ii]])
  fc=signif.avg[,index2[ii]]/signif.avg[,1]
  write.csv(cbind(sort(benoit.tukey.pval[,index[ii]]),fc[xx],signif.avg[xx,]),file=paste(name[ii], "p vals.csv"))
  jpeg(paste(name[ii],"Q-Q plot.jpeg"))
  qqnorm(benoit.tukey.pval[,index[ii]])
  qqline(benoit.tukey.pval[,index[ii]])
  dev.off()
}


############### evaluating probesets ########################################
genes=gsub(".[0-9]{1,2}$","",dimnames(signif)[[2]])

probeset.mean=matrix(nrow=nrow(signif),ncol=length(unique(genes)))
probeset.sd=matrix(nrow=nrow(signif),ncol=length(unique(genes)))
for(ii in 1:nrow(signif)){
  probeset.mean[ii,]=tapply(as.matrix(signif)[ii,],genes,mean)
  probeset.sd[ii,]=tapply(as.matrix(signif)[ii,],genes,sd)
}
colnames(probeset.mean)=sort(unique(genes))
colnames(probeset.sd)=sort(unique(genes))



########## Comparing this list to the C1q list from Galvan #########################

galvan=read.csv("Control vs C1q Treatment p vals.csv")

apop.v.control=sort(benoit.tukey.pval[,3])
c1qapop.v.control=sort(benoit.tukey.pval[,7])
necro.v.control=sort(benoit.tukey.pval[,14])
c1qnecro.v.control=sort(benoit.tukey.pval[,10])
lps.control=sort(benoit.tukey.pval[,13])

######################################################################

cor.abund.lps=cor(lps.abund[1:28,])
install.packages("GGally")
library(GGally)

for (ii in 1:6){
  xx=cor(lps.abund[which(lps.avg[,ii]<200),which(group==unique(group)[1])])
  jpeg(paste("Correlation matrix low abundance genes",unique(group)[ii],".jpeg"))
  pairs(lps.abund[which(lps.avg[,ii]<200),which(group==unique(group)[ii])],main=paste("Correlations of",unique(group)[ii],"R=",substr(mean(xx),1,5)),pch=20,xlim=c(0,10),ylim=c(0,10))
  dev.off()
}

library(manipulate)
manipulate(plot(1:x), x = slider(5, 10))
