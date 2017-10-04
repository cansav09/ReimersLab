## Step 1) Load original GEO data files
setwd('~/Desktop/Neurogenomics Data/CellProportionsProject/GSE68559_RAW/')
fnames=dir()

tmp=read.table(fnames[1])
brain.region = matrix(nr=nrow(read.table(fnames[1])),nc=length(dir()))

rownames(brain.region)<-tmp[,1]
for ( ii in 1:ncol(brain.region)) {
  tmp= read.table(fnames[ii])
  if (length(tmp[[ii+1]]) != nrow(brain.region)) {warning('sample',ii,'has',);next}
  brain.region[,ii] <- tmp[,2] 
}

setwd('..')
dimnames(darmanis)[[2]]<-sub('_.*C','.C',sub('.csv.gz','',fnames))
saveRDS( darmanis, file='brain.region.Webb.RDS')
