
############################## Master Project file for Martinez et al Data #######################
# Candace Savonen
# Date Created: Feb 11, 2016
#
# Overall Project Summary: To identify gene expression signatures for the traditional macrophage
# treatments of LPS and IFN-gamma from Martinez et al and troubleshoot the data analysis workflow. 
#
#
############################## List of .R files for this directory ###############################
##################################################################################################
# File Name: "MartinezDataGEOfile.R"
# Task: Downloads data and annotation from GEO and does ANOVA and Hotelling's test to identify genes
# differentially expressed across M0, M1 and M2 treated macrophages. Makes MA plots of all combinations 
# of experimental groups. Set FC to your desired FC cut off. 
# Last Update: Mar 30 2017
# Problems: NEED to run this code before using any of the other code in this project. 
#
##################################################################################################
# File Name: "ComparingZhangwithMstateGenes.R"
# Task: Basic file modeled after GEO R code template. Compares the differentially expressed M state 
# genes identified with Zhang RNA-Seq Expression with the different cell types. Edits out genes 
# that are being expressed in other cell 
# types besides microglia 
# Last Update: Mar 21 2017
# Problems:
# Findings: Tried taking out M state genes that are expressed in other cell types and then 
# redoing PCA. 


##################################################################################################
# File Name: "PCAproblem.R"
# Task:Investigating why the PCA score correlations are always so inflated. 
# Last Update: Mar 21 2017
# Problems:
#
##################################################################################################
# File Name: "GOontologydownloaderM1.R"
# Task:Downloads GO ontology and matches it to each gene in the M1 List
# Last Update: Feb 2 2017
# Problems:Depending on the gene list that you try to download from HUGO the URL may not work. 
#
##################################################################################################
# File Name: "GOontologydownloaderSynapse.R"
# Task:Downloads GO ontology and matches it to each gene in the M1 List
# Last Update: Feb 2 2017
# Problems:It isn't completed yet. Needs rankings to be made for the different types of evidence.
#
##################################################################################################
# File Name: "MStateGibbsWorkflow.R"
# Task: Takes out the M state genes' data from Gibbs et al data and downloads gene of interest
# info and does PCA then correlates M1 PCA scores with "gene family of interest" Which you can 
# download the HUGO gene list from it.
# Last Update:  Mar 3 2017
# Problems:Depending on the gene list that you try to download from HUGO the URL may not work. 
#
##################################################################################################
# File Name: "HierarchicalCluster.R"
# Task: Clusters genes of interest expression's based on Gibb's illumina probe data
# Last Update: Jan 20 2017
# Problems: Just need to make sure to set it to the genes you want and change things accordingly
#
##################################################################################################
# File Name: "SynapseDB.R"
# Task: Loads into R the lists of gene that are identified by Syntaptome DB as involved in synapse
# structure. Then does a hierarchical cluster using Gibbs data. 
# Last Update: Mar 27 2017
# Problems: Not completely edited yet. 
#
##################################################################################################
# File Name: "GABAGluGibbs.R"
# Task: Examines the relationship between GABA Glutamate and M1 related genes. Does the same functions
# as the PritzkerM1GABAGlu.R does 
# Last Update: May 8 2017
# Problems: Not completely edited yet. 
#
##################################################################################################
##################################################################################################
############################ File creation and organization code #################################
name="Martinez M0M1M2 Lists"
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
