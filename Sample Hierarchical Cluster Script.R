
### Script makes a heirarchical cluster for samples and keeps the ID


sample.IDs= "Put a vector that has the individual sample IDs you want to show "
groups= "Put a vector that has group that each sample isa part of in the order that the samples are entered in the data set"
data= "put a matrix with the samples in each column and the rows each genes' values"


library(dplyr)
library(ggplot2)
library(plyr)
library(dendextend)

color.group=mapvalues(groups,from = c(unique(groups)), to = c("blue","red","brown","green","black","purple"))

jpeg("Cluster Analysis of Samples w IDs.jpeg",width=700,height=500)
par(mar=c(10,3,3,3))
dend <- t(data) %>%  scale %>% dist %>% hclust %>% as.dendrogram
dend %>% plot
labels(sample.IDs)
dend %>% set("labels_col", color.group) %>% set("labels_cex", 1)  %>% plot(main = "Cluster of Cell Types/Treatments")

legend(175,c(unique(groups)),lty=c(1,1), lwd=c(2.5,2.5),col=c(color.group),cex=.75) 
dev.off()



