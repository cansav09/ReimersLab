###### MANOVA and Hotelling-Lawley trace with flexibility for p>>n cases ##############################

### Overall these data are representing 3 different treatment groups with 4 samples each. 
### Only probes with expression >0 are included. 
### non.necrotic.groups is a dataset that has probe x sample matrix. It has the necrotic group samples removed.
### benoit.groups is a vector that holds the group that each respective sample is in. It is subsetted to removed the necrotic groups 
### benoit.gene.names.abund is the vector that holds the corresponding gene names of each row of data of the non.necrotic.groups matrix


### There are 17,229 genes represented with a wide range of how many probes cover them. 
#Descriptive stats for number of probes per gene: 
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1     3         6         8     11      107
### Around half of the genes have p>>n which is why we need to integrate the variance regularization step for those p>>n genes



fstat=c()#### A bunch of empty vectors to store stats in. 
hotelling.trace=c()
hotelling.pval.test=c()
num.probes=c() # number of probes per gene
df1=c()
df2=c()
for(ii in 1:length(unique(benoit.gene.names.abund))){ #For each gene name, 
  xx=which(benoit.gene.names.abund==unique(benoit.gene.names.abund)[2]) 
  ### Finds row number index for all the probes for a particular gene
  num.probes[ii]=p=length(xx) ### Takes note of how many probes are in each set. Sets "p" as the dimensions
  xx=as.matrix(non.necrotic.groups[xx,]) #Takes the data from the row indices indicating previously indicated 
  colnames(xx)=benoit.groups[1:12]

  ## non.necrotic groups is a dataset with all the probes (for which there is an average signal of 2 or greater) 
  #for all 12 samples that we are investigating
  n=ncol(xx) ## number of samples which is always 12 for this dataset.
  if(p<2){ #When there is only 1 probe for a gene, This set of code does normal ANOVA and by passes the rest of the code to the next loop
    xx=summary(aov(xx~as.factor(benoit.groups[1:12]),data=as.data.frame(xx)))## summarizes aov object into the report stats.
    # note: The only reason this is set to groups [1:12] is that we removed the necrotic groups
    fstat[ii]=xx[[1]][["F value"]][1] ### Stores F statistic in the same vector as genes with multiple probes
    hotelling.pval.test[ii]=xx[[1]][["Pr(>F)"]][1] ## Stores p value for 1 probe sets with the other pvalues
    df1[ii]=NA # These are just place holders because the last ten genes only have one probe, and we need to still have the same number of rows
    df2[ii]=NA 
    hotelling.trace[ii]=NA
    next
  }else{ # For genes with more than one probe: 
    xx=manova(t(xx)~as.factor(benoit.groups[1:12]))
    ss.resid =crossprod(xx$residuals) ## Finds the crossproduct of residuals matrix
    ss.effect=crossprod(xx$effects[2:4,]) ## Finds the crossproduct of effect matrix, leaves out the intercept, keeps only the effects for each group
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



####################### Other descriptive statistics for the probesets for each gene#########################
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