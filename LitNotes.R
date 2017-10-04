
library(pdftools)

setwd( "/Users/cansav091/Desktop/Papers")
folders=grep(".pdf",dir(),invert=TRUE,value=TRUE)
folders=grep(".PDF",folders,invert=TRUE,value=TRUE)

papers=grep(".pdf",dir(),value=TRUE)
papers=c(papers,grep(".PDF",dir(),value=TRUE))
txt <- pdf_text(papers[ii])

category=rep("NA",length(papers))
for(ii in 1:length(folders)){
  setwd(paste0("/Users/cansav091/Desktop/Papers/",folders[ii]))
  xx=grep(".pdf",dir(),value=TRUE)
  papers=c(papers,xx)
  category=c(category,rep(folders[ii],length(xx)))
}

chart=as.data.frame(cbind(as.character(papers),category))
chart$V1=gsub(",","",as.character(chart$V1))
colnames(chart)=c("File Name","Folder")

setwd( "/Users/cansav091/Desktop/Current Projects ")
write.csv(chart,file="PaperDirectory.csv",quote=FALSE)

