
############################## Master Project file for Tung et al Data #######################
# Candace Savonen
# Date Created: Mar 31, 2016
#
# Overall Project Summary: To identify gene expression markers that are indicative of increases in 
# numbers of dendritic spines.  Use both mass spec data from Frese et al and gene expression data. 
#
#
############################## List of .R files for this directory ###############################
##################################################################################################
# File Name: "FreseDataImportClean.R"
# Task: Import data files that I got from Frese et al, 2017 supplement section. 
# Last Update: Mar 31 2017
# Problems: NEED to run this code before using any of the other code in this project. 
#

##################################################################################################
############################ File creation and organization code #################################
name="TungStressData"
mainDir=paste0("/Users/cansav091/Desktop/Neurogenomics Data/",name)
input=paste0(mainDir,"/","Original Data")
if(dir.exists(input)==FALSE){
  dir.create(input)
}

############ Directory Commands ##################
# Type this before referencing something in the main directory
setwd(mainDir)
# Type this before referencing original data folder
setwd(input)
