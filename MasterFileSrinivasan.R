
############################## Master Project file for Srinivasan et al Data #######################
# Candace Savonen
# Date Created: Aug 8, 2016
#
# Overall Project Summary: To identify gene expression signatures for the major brain cell types when
# treated wtih LPS.
#
#
############################## List of .R files for this directory ###############################
##################################################################################################
# File Name: "LPSCellTypes.R"
# Task: Downloads data and annotation from GEO and does ANOVA and Hotelling's test to identify genes
# differentially expressed across M0, M1 and M2 treated macrophages. Makes MA plots of all combinations 
# of experimental groups. Set FC to your desired FC cut off. 
# Last Update: Mar 30 2017
# Problems: NEED to run this code before using any of the other code in this project. 
#
##################################################################################################
############################ File creation and organization code #################################
name="CellTypesAndLPSSrnivasan"
mainDir=paste0("/Users/cansav091/Desktop/Neurogenomics Data/",name)
input=paste0(mainDir,"/","Original Data")
if(dir.exists(input)==FALSE){
  dir.create(input)
}

############ Directory Commands ##################
# Type this before referencing something in the main directory
setwd(mainDir)
# Type this before referencing original data
setwd(input)
