##### Alternative ANOVA for Benoit et al data 


benoit.abund=benoit.data[which(gene.avg> 0),]
benoit.gene.names.abund=gene.name.benoit[which(gene.avg> 0)]
benoit.groups= c(rep("control",3),rep("LPS",3),rep("Apoptotic LPS",3),rep("C1q Apoptotic LPS",3),rep("Necrotic LPS",3),rep("C1q necrotic LPS",3))


####### Let's see if the necrotic samples are indeed more variable #############
####### Compare the two variances ##############################################

non.necrotic.groups=benoit.abund[,1:12]
necrotic.groups=benoit.abund[,13:18]

# F test to compare two variances
f.test=matrix(ncol=2,nrow=length(benoit.gene.names.abund))
for(ii in 1:length(benoit.gene.names.abund)){
  f.test[ii,1]=var.test(non.necrotic.groups[ii,],necrotic.groups[ii,])$statistic
  f.test[ii,2]=var.test(non.necrotic.groups[ii,],necrotic.groups[ii,])$p.value
}
length(which(f.test[,2]<.05))/length(benoit.gene.names.abund)### 14% are unequal variance

f.test.eg=matrix(ncol=2,nrow=length(benoit.gene.names.abund))
for(ii in 1:length(benoit.gene.names.abund)){
  xx=sample(1:18,9)
  yy=c(1:18)[is.na(match(1:18,xx))]
  zz=var.test(benoit.abund[ii,xx],benoit.abund[ii,yy])
  f.test.eg[ii,1]=zz$statistic
  f.test.eg[ii,2]=zz$p.value
}
length(which(f.test.eg[,2]<.05))/length(benoit.gene.names.abund)### .07 are unequal variance


###### MANOVA and Hotelling-Lawley trace with flexibility for p>>n cases ##############################
fstat=c()
hotelling.trace=c()
hotelling.pval.test=c()
num.probes=c()
df1=c()
df2=c()
for(ii in 1:length(unique(benoit.gene.names.abund))){
  xx=which(benoit.gene.names.abund==unique(benoit.gene.names.abund)[ii])
  num.probes[ii]=p=length(xx)
  xx=as.matrix(non.necrotic.groups[xx,])
  n=ncol(xx)
  if(p<2){
    xx=summary(aov(xx~as.factor(benoit.groups[1:12]),data=as.data.frame(xx)))
    fstat[ii]=xx[[1]][["F value"]][1]
    hotelling.pval.test[ii]=xx[[1]][["Pr(>F)"]][1]
    df1[ii]=NA### will be between 9 and 18 This will be really big with a really big p
    df2[ii]=NA #### will be between 5 and 20 This will become negative and really small with a large p
    hotelling.trace[ii]=NA
    next
  }else{
    xx=manova(t(xx)~as.factor(benoit.groups[1:12]))
    ss.resid =crossprod(xx$residuals)
    ss.effect=crossprod(xx$effects[2:4,])
    cc=try(solve(ss.resid), silent=T) 
    if(is(cc,"try-error")) {
      A1=ss.effect%*%condreg(ss.resid,8)$invS
    }else{
      m=4 ### Number of groups +1 for intercept
      A1=ss.effect%*%solve(ss.resid)
    }### four groups
    u=(n-m-p-1)/2
    u=ifelse(u<0,0,u)
    t=(abs(p-m+1)-1)/2
    s=min(c(p,m-1))### s will always be 3. 
    df1[ii]=s*(2*t+s+1)### will be between 9 and 18 This will be really big with a really big p
    df2[ii]=2*(s*u+1) #### will be between 5 and 20 This will become negative and really small with a large p
    hotelling.trace[ii]=tr(A1)
    fstat[ii]=(tr(A1)/s)*(df2[ii]/df1[ii])
    hotelling.pval.test[ii]=1-pf(fstat[ii], df1[ii], df2[ii])
  }
}


length(which(hotelling.pval.test<.05)) #N=9801
pval.fdr=p.adjust(hotelling.pval.test,method="hochberg") ##15 significant genes p<.05


sd.probes=c()
avg.group.signal=matrix(nrow=length(unique(benoit.gene.names.abund)),ncol=4)
for(ii in 1:length(unique(benoit.gene.names.abund))){
  xx=which(benoit.gene.names.abund==unique(benoit.gene.names.abund)[ii])
if(length(xx)<2){
  avg.group.signal[ii,]=benoit.avg[xx,1:4]
  sd.probes[ii]=NA
  next
}else{
avg.group.signal[ii,]=apply(benoit.avg[xx,1:4],2,mean)
sd.probes[ii]=sd(apply(benoit.abund[xx,1:4],1,mean))
}
}

signif=cbind(hotelling.trace,fstat,hotelling.pval.test,pval.fdr,num.probes,sd.probes,avg.group.signal)
rownames(signif)=unique(benoit.gene.names.abund)
colnames(signif)=c("Hotelling Trace","F-stat","p value","FDR adj p val","Num of Probesets","SD of probsets","Controls","LPS","Apoptotic LPS","C1q Apoptotic")

write.csv(signif[order(signif[,4]),],file="Benoit C1Q Genes.csv")
plot(num.probes,-log10(pval.fdr))

xx=signif[which(signif[,4]<.05),]
xx=prcomp(t(xx[,8:10]))
plot(xx$rotation)






