
# Removing Outliers from MState Gene List

name="OutlierTestMstateGenes"
outDir=paste0(mainDir,"/",name)
outDirDay=paste0(mainDir,"/",name,"/",as.character(Sys.Date()))
outDirGraphsDay=paste0(outDir,"/Graphs","/",as.character(Sys.Date()))
if(dir.exists(outDir)==FALSE){
  dir.create(outDir)
}
if(dir.exists(outDirDay)==FALSE){
  dir.create(outDirDay)
  dir.create(outDirGraphsDay)
}

### Let's check if there might be outliers that are pulling these data
for(ii in 1:ncol(Mstategib.edit)){
    jpeg(paste0(Mstategib.gene.name.edit[ii],"histogram Gibbs.jpeg"))
    hist(Mstategib.edit[,ii],xlab=Mstategib.gene.name.edit[ii])
    dev.off()
}

# Doing a simple Z score Grubb's test to identify if any of the genes in the M state list 
# show irregular distributions. 
mstateZscore=matrix(nrow=291,ncol=ncol(Mstategib.edit))
for(ii in 1:ncol(Mstategib.edit)){
  pop_sd= sd(Mstategib.edit[,ii])*sqrt((290-1)/(290))
  pop_mean=mean(Mstategib.edit[,ii])
  mstateZscore[,ii]=( Mstategib.edit[,ii]- pop_mean) / pop_sd
}

outliers=which(pnorm(mstateZscore)<.05,arr.ind = TRUE)

outliers.sample=table(outliers[,1])
as.numeric(names(sort(table(outliers[,2]),decreasing=TRUE)))


Mstategib.gene.name.edit[as.numeric(names(sort(table(outliers[,2]),decreasing=TRUE)))]

cor(cbind(gene.fam.pca$x[,1],Mstatepca.edit$x[,1],house.pca$x[,1],rand.pca$x[,1]))
