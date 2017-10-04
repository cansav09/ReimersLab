### Comparing datasets 


galvan=read.csv("Galvan et al C1q Hotelling data.csv",header=TRUE)
benoit=signif

tmp=match(galvan[,1],rownames(benoit))
galvan=galvan[!is.na(tmp),]
benoit=benoit[tmp[!is.na(tmp)],]

galvan.signif=galvan[which(galvan[,7]<.2),]
benoit.signif=benoit[which(benoit[,7]<.2),]

tmp=match(galvan.signif[,1],rownames(benoit.signif))
galvan.signif=galvan.signif[!is.na(tmp),]
benoit.signif=benoit.signif[tmp[!is.na(tmp)],]

plot(galvan[,7],benoit[,7])

plot(-log10(benoit.signif[,3]),-log10(galvan.signif[,4]))
     