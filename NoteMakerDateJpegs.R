library(chron)
library(rtf)
setwd("/Users/cansav091/Desktop/Neurogenomics Data")
projects=dir()
output=paste0("/Users/cansav091/Desktop/Current Projects /Notes")
proj.names=c()
proj.notes=c()
proj.notes[[1]]=paste0("Notes on all Neurogenomics projects as of ",Sys.Date())
for(ii in 1:length(projects)){
  setwd(paste0("/Users/cansav091/Desktop/Neurogenomics Data/",projects[ii]))
  #total.notes=c(total.notes,paste0("Project Folder: ",projects[ii]," Started: ",proj.start))
  files=dir()
  files=grep(".R$",files,value=TRUE)
  if(length(files)>0){
    proj.names=c(proj.names,rep(projects[ii],length(files)))
    proj.notes=c(proj.notes,paste0(projects[ii],": ",files))
  }else{
    next
  }
}
proj.vector=proj.notes
proj.notes=as.list(proj.notes)

mod.dates=c()
all.notes=c()
lineswplots.dir=c()
plot.names=c()
plots.dir=c()
for(ii in 1:length(projects)){
  setwd(paste0("/Users/cansav091/Desktop/Neurogenomics Data/",projects[ii]))
  files=dir()
  plot.files=grep(".jpeg$",files,value=TRUE)
  files=grep(".R$",files,value=TRUE)
  if(length(files)>0){ 
    for(jj in 1:length(files)){
      mod.dates=c(mod.dates,as.character(file.mtime(paste0("/Users/cansav091/Desktop/Neurogenomics Data/",projects[ii],"/",files[jj]))))
      xx=readLines(files[jj],warn=FALSE)
      lineswplots=grep(".jpeg",xx)
      if(length(c(lineswplots,plot.files)>0)){
      plots.dir=c(lineswplots.dir,paste0("/Users/cansav091/Desktop/Neurogenomics Data/",projects[ii],"/",substr(files[jj],1,nchar(files[jj])-2)))
      for(qq in 1:length(plots.dir)){
        if(dir.exists(plots.dir[qq])){
          for(ww in 1:length(lineswplots)){
            setwd(plots.dir[qq])
            plot.names=grep(".jpeg",dir(),value=TRUE)
          }
          }
        }
      if(length(plot.files>0)){
        lineswplots.dir=c(lineswplots.dir,paste0("/Users/cansav091/Desktop/Neurogenomics Data/",projects[ii],"/",plot.files))
      }
      if(length(plot.names)>0){
      lineswplots.dir=c(lineswplots.dir,paste0(plots.dir,"/",plot.names))
      }
      setwd(paste0("/Users/cansav091/Desktop/Neurogenomics Data/",projects[ii]))
      plot.names=c()
      plots.dir=c()
       }
      xx=xx[unique(c(grep("\\#",xx),grep(".jpeg",xx),grep("write.",xx)))]
      xx[unique(grep(".jpeg",xx),grep("write",xx))]=paste0("\t",xx[unique(grep(".jpeg",xx),grep("write.",xx))])
      xx=grep("\\#",xx,value=TRUE)
      xx=unlist(strsplit(xx,"\\#"))
      xx=grep("%",xx,value=TRUE,invert=TRUE)
      xx=grep("\\+",xx,value=TRUE,invert=TRUE)
      xx=grep("=",xx,value=TRUE,invert=TRUE)
      xx=grep("[:alpha:]+",xx,value=TRUE)
      xx=grep("\\{",xx,value=TRUE,invert=TRUE)
      xx=as.data.frame(xx)
      xx=apply(xx,2,function(x) paste("\t \t ", x))
      fmod=as.character(paste0("\t Code:",files[jj]," Notes as of ",file.mtime(paste0("/Users/cansav091/Desktop/Neurogenomics Data/",projects[ii],"/",files[jj])) ))
      xx=rbind(fmod,xx)
      yy=grep(paste0(projects[ii],": ",files[jj],"$"),proj.vector)
      proj.notes[[yy]]=xx
    }
  }else{
    next
  }
}

names(proj.notes)=c("Neurogenomics Data",proj.names)
mod.dates=sort(mod.dates, method = "radix",index.return=TRUE,decreasing=TRUE)
proj.notes.ordered=c(proj.notes[1],proj.notes[2:length(proj.notes)][mod.dates$ix])
lineswplots=grep("jpeg",proj.notes)

setwd(output)
write.csv(unlist(proj.notes.ordered),file=paste0("Notes on Code from All Projects as of ",Sys.Date()),row.names=FALSE,quote=FALSE)
files=dir()
old=dir()
for(ii in 1:length(files)){
  old[ii]=file.mtime(paste0("/Users/cansav091/Desktop/Current Projects /Notes/",files[ii]))
}
new=readLines(paste0("Notes on Code from All Projects as of ",Sys.Date()))

old=readLines(files[order(old,decreasing=TRUE)][2])
old=grep("fs20",old,value=TRUE,invert=TRUE)
old=substr(old,1,nchar(old)-5)
old=old[11:length(old)]

changedlines=which(is.na(match(new,old)))
lineswplots=grep("jpeg",new)
new[lineswplots]

rtf<-RTF(paste0("Notes on Code from All Projects as of ",Sys.Date()),width=8.5,height=11,font.size=10,omi=c(1,1,1,1))
for(ii in 2:length(new)){
  if(!is.na(match(ii,lineswplots))){
    addPng(rtf,)
  }
  if(!is.na(match(ii,changedlines))){
  addText(rtf,new[ii],bold=TRUE)
  }else{
    addParagraph.RTF(rtf,new[ii])
  }
}
done.RTF(rtf)

