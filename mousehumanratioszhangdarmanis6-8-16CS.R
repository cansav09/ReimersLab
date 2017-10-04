### Still trying to figure out why the Darmanis and Zhang dataset aren't matching
### Let's look at fold changes between datasets and compare them across cell types and see if that tells us anything

### Let's use datasets that are using abundant enough based on the gene avg across all cell types
human.rpkm.abun=human.rpkm.avg[which(apply(human.rpkm.avg,1,mean)>1),]
zhang.fpkm.abun=zhang.fpkm.avg[which(apply(zhang.fpkm.avg,1,mean)>1),]

match(dimnames(human.rpkm.avg)[[1]],toupper(dimnames(zhang.fpkm.avg)[[1]]))## 13134/15500 genes in human data 13134/22462 zhang
## 7573/9294 genes in human data 7573/12440 zhang

zhang.overlap=zhang.fpkm.abun[match(dimnames(human.rpkm.abun)[[1]],toupper(dimnames(zhang.fpkm.abun)[[1]])),]## 7573/9294 genes in human data 7573/12440 zhang
human.overlap=human.rpkm.abun[match(toupper(dimnames(zhang.fpkm.abun)[[1]]),dimnames(human.rpkm.abun)[[1]]),]

human.overlap=human.overlap[!is.na(dimnames(human.overlap)[[1]]),]
zhang.overlap=zhang.overlap[!is.na(dimnames(zhang.overlap)[[1]]),]

### Checking that all the genes match by the order the datasets are in.
sum(!is.na(match(dimnames(human.overlap)[[1]],toupper(dimnames(zhang.overlap)[[1]]))))

neuron.ratio=human.overlap[,6]/zhang.overlap[,7]
astrocyte.ratio=human.overlap[,1]/zhang.overlap[,1]

plot(neuron.ratio,astrocyte.ratio)


