
############################## Master Project file for "       " Data #######################
# Candace Savonen
# Date Created:

# Overall Project Summary:


############################## List of .R files for this directory ###############################
##################################################################################################
# File Name: "MartinezDataGEOfile.R"
# Task:
# Last Update: 
# Problems:

##################################################################################################
# File Name:
# Task:
# Last Update: 
# Problems:

# To keep things organized puts output files in a folder based on the date
name="NAME"
out=setwd(outDir)
home=setwd(mainDir)
mainDir=paste0("/Users/cansav091/Desktop/Neurogenomics Data/",name)
input=paste0(mainDir,"/","Original Data")
outDir=paste0(mainDir,"/","Output/",as.character(Sys.Date()))
outDirGraphs=paste0(outDir,"/","Graphs")
if(dir.exists(file.path(mainDir, "/Output"))==FALSE){
    dir.create(file.path(mainDir, "/Output"))
    dir.create(input)
}
if(dir.exists(outDir)==FALSE){
  dir.create(outDir)
}
if(dir.exists(outDirGraphs)==FALSE){
  dir.create(outDirGraphs)
}

############ Directory Commands ##################
# Type this before creating any output files
out=setwd(outDir)
# Type this before creating any graphs
graph=setwd(outDirGraphs)
# Type this before referencing something in the main directory
home=setwd(mainDir)
# Type this before referencing original data
org=setwd(input)


