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
jpeg.locs=c()
jpeg.proj=c()
for(ii in 1:length(projects)){
  setwd(paste0("/Users/cansav091/Desktop/Neurogenomics Data/",projects[ii]))
  files=dir()
  files=grep(".R$",files,value=TRUE)
  if(length(files)>0){ 
    for(jj in 1:length(files)){
      mod.dates=c(mod.dates,as.character(file.mtime(paste0("/Users/cansav091/Desktop/Neurogenomics Data/",projects[ii],"/",files[jj]))))
      xx=readLines(files[jj],warn=FALSE)
      xx=xx[unique(c(grep("\\#",xx),grep(".jpeg",xx),grep("write.",xx)))]
      xx[unique(grep(".jpeg",xx),grep("write.",xx))]=paste0("\t \t",xx[unique(grep(".jpeg",xx),grep("write.",xx))])
      xx=unlist(strsplit(xx,"\\#"))
      xx=grep("^$",xx,value=TRUE, invert=TRUE)
      xx=grep("%",xx,value=TRUE,invert=TRUE)
      xx=grep("\\+",xx,value=TRUE,invert=TRUE)
      xx=grep("\\[",xx,value=TRUE,invert=TRUE)
      xx=grep("=",xx,value=TRUE,invert=TRUE)
      xx=grep("[:alpha:]+",xx,value=TRUE)
      xx=grep("\\{",xx,value=TRUE,invert=TRUE)
      xx=as.data.frame(xx)
      xx=apply(xx,2,function(x) paste("\t \t ", x))
      fmod=as.character(paste0("\t Code:",files[jj]," Notes as of ",file.mtime(paste0("/Users/cansav091/Desktop/Neurogenomics Data/",projects[ii],"/",files[jj])) ))
      xx=rbind(fmod,xx)
      yy=grep(paste0(projects[ii],": ",files[jj],"$"),proj.vector)
      proj.notes[[yy]]=xx
      
      if(length(grep(".jpeg",dir()))>0){
        jpeg.locs=c(jpeg.locs,paste0(rep(getwd(),nn),"/",grep(".jpeg",dir(),value=TRUE)))
        jpeg.proj=c(jpeg.proj,rep(projects[ii],nn))
      }
      
      jpeg.dir=paste0(getwd(),"/",files)
      for(ll in 1:length(jpeg.dir)){
        if(dir.exists(jpeg.dir[ll])==TRUE){
          setwd(jpeg.dir[ll])
          if(length(grep(".jpeg",dir()))>0){
          nn=length(grep(".jpeg",dir()))
          jpeg.locs=c(jpeg.locs,paste0(rep(jpeg.dir[ll],nn),"/",grep(".jpeg",dir(),value=TRUE)))
          jpeg.proj=c(jpeg.proj,rep(projects[ii],nn))
          }
        }
      }
    }
  }else{
    next
  }
}

names(proj.notes)=c("Neurogenomics Data",proj.names)
mod.dates=sort(mod.dates, method = "radix",index.return=TRUE,decreasing=TRUE)
proj.notes.ordered=c(proj.notes[1],proj.notes[2:length(proj.notes)][mod.dates$ix])
proj.notes.ordered=unlist(proj.notes.ordered)

setwd(output)
rtf<-RTF(paste0("Notes on Code from All Projects as of ",Sys.Date()),width=8.5,height=11,font.size=10,omi=c(1,1,1,1))
addParagraph.RTF(rtf,proj.notes.ordered)
done.RTF(rtf)

new <- readLines(paste0("Notes on Code from All Projects as of ",Sys.Date()))
g = grep("\\{", new,invert=TRUE, value = TRUE)
g = gsub("\\\\par}","", g)
g = grep("\\}", g,invert=TRUE, value = TRUE)
g = grep("\t",g,value=TRUE)
new=g

files=dir()
old=dir()
for(ii in 1:length(files)){
  old[ii]=file.mtime(paste0("/Users/cansav091/Desktop/Current Projects /Notes/",files[ii]))
}
old=readLines(files[order(old,decreasing=TRUE)][2])
g = grep("\\{", old,invert=TRUE, value = TRUE)
g = gsub("\\\\par}","", g)
g = grep("\\}", g,invert=TRUE, value = TRUE)
g = grep("\t",g,value=TRUE)
old=g

changedlines=which(is.na(match(new,old)))
lineswplots=grep("jpeg\\(",new)

jpeg.locs


rtf<-RTF(paste0("Notes on Code from All Projects as of ",Sys.Date()),width=8.5,height=11,font.size=10,omi=c(1,1,1,1))
for(ii in 1:length(new)){
   # if(!is.na(match(ii,lineswplots))){
    #  addPng(rtf,)
    #}
    if(!is.na(match(ii,changedlines))){
      addText(rtf,paste0("\n",new[ii]),bold=TRUE)
    }else{
      addParagraph.RTF(rtf,new[ii])
    }
  }
  done.RTF(rtf)
  
  
