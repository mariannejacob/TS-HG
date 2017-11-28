###preparations for statistics
####
######################################################################################
#This script was written by Marianne Jacob (orcid.org/0000-0002-9351-8014) 
#to analyse Illumina MiSeq data for an Arctic sediment bacterial dataset.
#
#This is a supplement script to "Rscript for data anlysis of time series.R" and serves for the preparation of subset data
#This script needs to be used with the appropriate files and scripts for analysis.
####14.02.2017###
######################################################################################
#
#location of Input files
setwd("//tux/mjacob/Documents/PhD/1. Projects/time series till 2015/Data analysis/Input files")

#location of scripts
setwd("//tux/mjacob/Documents/PhD/1. Projects/time series till 2015/Data analysis/supporting scripts")

######load/save workspace
#load("//tux/mjacob/Documents/PhD/1. Projects/time series till 2015/Data analysis/time serie 2003-2015.RData")
#load("//tux/mjacob/Documents/PhD/1. Projects/time series till 2015/Data analysis/time serie 2003-2015_new.RData")
load("//tux/mjacob/Documents/PhD/1. Projects/time series till 2015/Data analysis/time serie 2003-2016.RData")

#save.image("//tux/mjacob/Documents/PhD/1. Projects/time series till 2015/Data analysis/time serie 2003-2015.RData")
#save.image("//tux/mjacob/Documents/PhD/1. Projects/time series till 2015/Data analysis/time serie 2003-2015_new.RData")
save.image("//tux/mjacob/Documents/PhD/1. Projects/time series till 2015/Data analysis/time serie 2003-2016.RData")

###R packages needed
require(vegan)
#require(MASS)
#functions
source("Y:\\Documents\\PhD\\14. R\\software\\merger\\replicate_merger_ALk_consensus_RFI 1.0.r")


#read in OTU and Taxonomy file and Metadata file
rawOTU=read.table("OTU x Sample full table.txt",header=T, row.names=1)
rawTaxo=read.table("Taxo.path full table.txt",header=T, row.names=1, sep="\t")
rawTaxo_split=read.table("Taxo sep full table.txt",header=T, row.names=1, sep="\t")
rawMetadata_ENV=read.table("Metadata_ENV time series.txt",header=T, row.names=1, sep="\t") 




#################################################################
###remove Archaea, Eukarya ,etc.
#######remove unclassified sequences
levels(rawTaxo_split$Domain) #check if anything else than Bacteria is in the data
###[1] "Archaea"     "Bacteria"    "No Relative"


##"No Relative"
#how many sequences could not be classified
length(which(rawTaxo_split$Domain == "No Relative"))
#remove unclassified
rows=which(rawTaxo_split$Domain == "No Relative")
cleanOTU=rawOTU[-rows,] # remove all "No Relative"
cleanTaxo=as.data.frame(rawTaxo[-rows,])
cleanTaxo_split=rawTaxo_split[-rows,]

##"Archaea" 
#how many sequences could not be classified
length(which(cleanTaxo_split$Domain == "Archaea" ))
#remove unclassified
rows=which(cleanTaxo_split$Domain == "Archaea" )
cleanOTU=cleanOTU[-rows,] # remove all "Archaea" 
cleanTaxo=as.data.frame(cleanTaxo[-rows,])
cleanTaxo_split=cleanTaxo_split[-rows,]
#######################################################################

###check abundance distribution of taxa per sample
CleanOTU_gr_0=(cleanOTU)>0
hist(as.matrix(cleanOTU[CleanOTU_gr_0]),breaks=100000,xlim=c(1,50)) #full table
hist(apply(as.matrix(cleanOTU[CleanOTU_gr_0]),1,sum),breaks=100000,xlim=c(1,50)) # subset for sum of sample
hist(as.numeric(apply(cleanOTU,1,sum)),breaks=1000000,xlim=c(1,50))  # subset for sum of sample

############remove singletons!!!!######
cleanOTU_woSSOabs=as.data.frame(cleanOTU[apply(cleanOTU,1,sum)>1,])
cleanTaxo_woSSOabs=as.data.frame(cleanTaxo[apply(cleanOTU,1,sum)>1,1])
cleanTaxo_split_woSSOabs=as.data.frame(cleanTaxo_split[apply(cleanOTU,1,sum)>1,])
###245 Taxa were removed!!!->1400 OTU in full dataset (with old DNA extracts from Josi...)


# ################make data relative
# cleanOTUrel=decostand(cleanOTU,"total",MARGIN=2)
# apply(cleanOTUrel,2,sum) #check if samples =1
# cleanOTU_woSSOabs_rel=decostand(cleanOTU_woSSOabs,"total",MARGIN=2)
# apply(cleanOTU_woSSOabs_rel,2,sum) #check if samples =1
# #######################
# 
# ###############subsets####################
# #################################################
# #subset depth transect 2010 and Josi data
# OTU_comparison_raw=cleanOTU_woSSOabs[,c(60:64,133:137)] ##took oTU woSSOabs!!!!
# OTU_comparison=OTU_comparison_raw[apply(OTU_comparison_raw,1,sum)!=0,] #remove empty OTU #926 OTU
# OTU_comparison_rel=decostand(OTU_comparison,"total",MARGIN=2) #make relative
# apply(OTU_comparison_rel,2,sum) #check if samples =1
# ##Taxonomy
# Taxo_comparison=as.data.frame(cleanTaxo_woSSOabs[apply(OTU_comparison_raw,1,sum)!=0,1])
# Taxo_comparison_split=as.data.frame((cleanTaxo_split_woSSOabs[apply(OTU_comparison_raw,1,sum)!=0,]))
# 

###############################################
##new Time series dataset##
Index=  c(1:132) ###list with samples to be used
OTU=cleanOTU_woSSOabs
Taxo=cleanTaxo_woSSOabs
Taxo_split=cleanTaxo_split_woSSOabs
Env=rawMetadata_ENV
c(nrow(OTU),nrow(Taxo),nrow(Taxo_split),nrow(Env)) #check if dimension fit

##OTU
###calculations###  
OTU_subset_raw=OTU[,Index] #make subset #132 samples
OTU_subset=OTU_subset_raw[apply(OTU_subset_raw,1,sum)!=0,] #remove empty OTU #1398 OTU
OTU_subsetrel=decostand(OTU_subset,"total",MARGIN=2) #make relative
apply(OTU_subsetrel,2,sum) #check if samples =1

##Taxonomy
Taxo_subset=as.data.frame(Taxo[apply(OTU_subset_raw,1,sum)!=0,1])
Taxo_subset=droplevels(Taxo_subset)
Taxo_subset_split=as.data.frame((Taxo_split[apply(OTU_subset_raw,1,sum)!=0,]))

##clean SSOabs
OTU_subsetrel_woSSOabs=as.data.frame(OTU_subsetrel[apply(OTU_subset,1,sum)>1,])
Taxo_subset_woSSOabs=as.data.frame(Taxo_subset[apply(OTU_subset,1,sum)>1,1])
Taxo_subset_split_woSSOabs=as.data.frame(Taxo_subset_split[apply(OTU_subset,1,sum)>1,])
OTU_subset_woSSOabs=as.data.frame(OTU_subset[apply(OTU_subset,1,sum)>1,])
###1394 Taxa/ OTU in full dataset 

####ENv data
cleanENV=Env[Index,]

###cleaning####
c(nrow(OTU_subsetrel_woSSOabs),nrow(Taxo_subset_woSSOabs),
    nrow(Taxo_subset_woSSOabs),nrow(OTU_subset_woSSOabs), nrow(cleanENV))#check dim
OTU_0316_woSSOabs=OTU_subset_woSSOabs
OTU_0316_woSSOabs_rel=OTU_subsetrel_woSSOabs
Taxo_0316_woSSOabs=Taxo_subset_woSSOabs
Taxo_0316_woSSOabs_split=Taxo_subset_split_woSSOabs
Metadata_Env_0316=cleanENV
rm(OTU,OTU_subset, OTU_subset_raw,Taxo_subset,Taxo_subset_split,
    Taxo,Taxo_split,OTU_subsetrel_woSSOabs,Taxo_subset_woSSOabs,Taxo_subset_split_woSSOabs
    ,cleanENV, Env, OTU_subset_woSSOabs) #remove matrices for calculations


##make averages per year for the time series dataset
#averages per year
OTU_0316_average_woSSOabs = Merging_ALk(t(rbind(OTU_0316_woSSOabs,Metadata_Env_0316$Year)), k=1)
OTU_0316_average_woSSOabs_rel=decostand(OTU_0316_average_woSSOabs,"total",MARGIN=1) #make relative
colnames(OTU_0316_average_woSSOabs)=Taxo_0316_woSSOabs$`Taxo_subset[apply(OTU_subset, 1, sum) > 1, 1]`
colnames(OTU_0316_average_woSSOabs_rel)=Taxo_0316_woSSOabs$`Taxo_subset[apply(OTU_subset, 1, sum) > 1, 1]`
###############################################
##subset for comparison with old 454 data
##only N and HGIV stations from 2003-2009 without 2005
Index=  c(3,4,9,14,15,16,27,28,29,32,36,37,41,45,46,47,52,55,56,57,58) ###list with samples to be used
OTU=OTU_0316_woSSOabs
Taxo=Taxo_0316_woSSOabs
Taxo_split=Taxo_0316_woSSOabs_split
c(nrow(OTU),nrow(Taxo),nrow(Taxo_split)) #check if dimension fit#1394 rows


###calculations###  
OTU_subset_raw=OTU[,Index] #make subset #21 samples
OTU_subset=OTU_subset_raw[apply(OTU_subset_raw,1,sum)!=0,] #remove empty OTU #1036 OTU
OTU_subsetrel=decostand(OTU_subset,"total",MARGIN=2) #make relative
apply(OTU_subsetrel,2,sum) #check if samples =1

##Taxonomy
Taxo_subset=as.data.frame(Taxo[apply(OTU_subset_raw,1,sum)!=0,1])
Taxo_subset=droplevels(Taxo_subset)
Taxo_subset_split=as.data.frame((Taxo_split[apply(OTU_subset_raw,1,sum)!=0,]))

##clean SSOabs
OTU_subsetrel_woSSOabs=as.data.frame(OTU_subsetrel[apply(OTU_subset,1,sum)>1,])
Taxo_subset_woSSOabs=as.data.frame(Taxo_subset[apply(OTU_subset,1,sum)>1,1])
Taxo_subset_split_woSSOabs=as.data.frame(Taxo_subset_split[apply(OTU_subset,1,sum)>1,])
OTU_subset_woSSOabs=as.data.frame(OTU_subset[apply(OTU_subset,1,sum)>1,])
###909 Taxa/ OTU in full dataset 


###cleaning####
c(nrow(OTU_subset_woSSOabs),nrow(OTU_subsetrel_woSSOabs),nrow(Taxo_subset_woSSOabs),nrow(Taxo_subset_split_woSSOabs))#check dim
OTU_0309_woSSOabs=OTU_subset_woSSOabs
OTU_0309_woSSOabs_rel=OTU_subsetrel_woSSOabs
Taxo_0309_woSSOabs=Taxo_subset_woSSOabs
Taxo_0309_woSSOabs_split=Taxo_subset_split_woSSOabs
rm(OTU,OTU_subset, OTU_subset_raw,Taxo_subset,Taxo_subset_split,
   Taxo,Taxo_split,OTU_subsetrel_woSSOabs,Taxo_subset_woSSOabs,Taxo_subset_split_woSSOabs
   , OTU_subset_woSSOabs) #remove matrices for calculations


###################time series without 2005 and 2011
##
Index=  which(Metadata_Env_0316$Year!= "2005" & Metadata_Env_0316$Year!= "2011" &Metadata_Env_0316$Optional_label<30) ###list with samples to be used

OTU=OTU_0316_woSSOabs
Taxo=Taxo_0316_woSSOabs
Taxo_split=Taxo_0316_woSSOabs_split
Env=rawMetadata_ENV
c(nrow(OTU),nrow(Taxo),nrow(Taxo_split),nrow(Env)) #check if dimension fit


###calculations###  
OTU_subset_raw=OTU[,Index] #make subset #73 samples
OTU_subset=OTU_subset_raw[apply(decostand(OTU_subset_raw,"pa"),1,sum)!=0,] #remove empty OTU #1342 OTU
OTU_subsetrel=decostand(OTU_subset,"total",MARGIN=2) #make relative
apply(OTU_subsetrel,2,sum) #check if samples =1

##Taxonomy
Taxo_subset=as.data.frame(Taxo[apply(decostand(OTU_subset_raw,"pa"),1,sum)!=0,1])
Taxo_subset=droplevels(Taxo_subset)
Taxo_subset_split=as.data.frame((Taxo_split[apply(decostand(OTU_subset_raw,"pa"),1,sum)!=0,]))

##clean SSOabs
OTU_subsetrel_woSSOabs=as.data.frame(OTU_subsetrel[apply(OTU_subset,1,sum)>1,])
Taxo_subset_woSSOabs=as.data.frame(Taxo_subset[apply(OTU_subset,1,sum)>1,1])
Taxo_subset_split_woSSOabs=as.data.frame(Taxo_subset_split[apply(OTU_subset,1,sum)>1,])
OTU_subset_woSSOabs=as.data.frame(OTU_subset[apply(OTU_subset,1,sum)>1,])
###1314Taxa/ OTU in full dataset 

####ENv data
cleanENV=Env[Index,]


###cleaning####
c(nrow(OTU_subsetrel_woSSOabs),nrow(Taxo_subset_woSSOabs),
  nrow(Taxo_subset_woSSOabs),nrow(OTU_subset_woSSOabs), nrow(cleanENV))#check dim

OTU_0316_wo0511_woSSOabs=OTU_subset_woSSOabs
OTU_0316_wo0511_woSSOabs_rel=OTU_subsetrel_woSSOabs
Taxo_0316_wo0511_woSSOabs=Taxo_subset_woSSOabs
Taxo_0316_wo0511_split_woSSOabs=Taxo_subset_split_woSSOabs
Metadata_Env_0316_wo0511=cleanENV
rm(OTU,OTU_subset, OTU_subset_raw,Taxo_subset,Taxo_subset_split,
   Taxo,Taxo_split,OTU_subsetrel_woSSOabs,Taxo_subset_woSSOabs,Taxo_subset_split_woSSOabs
   ,cleanENV, Env, OTU_subset_woSSOabs) #remove matrices for calculations

#############################

#####################################################################
#depth transect
##subset with only stations from the depth transect
Index=  which(Metadata_Env_0316$Optional_label<10) ###list with samples to be used

OTU=OTU_0316_woSSOabs
Taxo=Taxo_0316_woSSOabs
Taxo_split=Taxo_0316_woSSOabs_split
Env=rawMetadata_ENV
c(nrow(OTU),nrow(Taxo),nrow(Taxo_split),nrow(Env)) #check if dimension fit


###calculations###  
OTU_subset_raw=OTU[,Index] #make subset #73 samples
OTU_subset=OTU_subset_raw[apply(decostand(OTU_subset_raw,"pa"),1,sum)!=0,] #remove empty OTU #1342 OTU
OTU_subsetrel=decostand(OTU_subset,"total",MARGIN=2) #make relative
apply(OTU_subsetrel,2,sum) #check if samples =1

##Taxonomy
Taxo_subset=as.data.frame(Taxo[apply(decostand(OTU_subset_raw,"pa"),1,sum)!=0,1])
Taxo_subset=droplevels(Taxo_subset)
Taxo_subset_split=as.data.frame((Taxo_split[apply(decostand(OTU_subset_raw,"pa"),1,sum)!=0,]))

##clean SSOabs
OTU_subsetrel_woSSOabs=as.data.frame(OTU_subsetrel[apply(OTU_subset,1,sum)>1,])
Taxo_subset_woSSOabs=as.data.frame(Taxo_subset[apply(OTU_subset,1,sum)>1,1])
Taxo_subset_split_woSSOabs=as.data.frame(Taxo_subset_split[apply(OTU_subset,1,sum)>1,])
OTU_subset_woSSOabs=as.data.frame(OTU_subset[apply(OTU_subset,1,sum)>1,])
###1272 Taxa/ OTU in full dataset 

####ENv data
cleanENV=Env[Index,]


###cleaning####
c(nrow(OTU_subsetrel_woSSOabs),nrow(Taxo_subset_woSSOabs),
  nrow(Taxo_subset_woSSOabs),nrow(OTU_subset_woSSOabs), nrow(cleanENV))#check dim

OTU_depth_woSSOabs=OTU_subset_woSSOabs
OTU_depth_woSSOabs_rel=OTU_subsetrel_woSSOabs
Taxo_depth_woSSOabs=Taxo_subset_woSSOabs
Taxo_depth_split_woSSOabs=Taxo_subset_split_woSSOabs
Metadata_Env_depth=cleanENV
rm(OTU,OTU_subset, OTU_subset_raw,Taxo_subset,Taxo_subset_split,
   Taxo,Taxo_split,OTU_subsetrel_woSSOabs,Taxo_subset_woSSOabs,Taxo_subset_split_woSSOabs
   ,cleanENV, Env, OTU_subset_woSSOabs) #remove matrices for calculations


  #############################################
#latitudinal transect

##subset with only stations from the latitudinal transect
Index=  which(Metadata_Env_0316$Optional_label>10 & Metadata_Env_0316$Optional_label< 30 | Metadata_Env_0316$Optional_label==4) ###list with samples to be used

OTU=OTU_0316_woSSOabs
Taxo=Taxo_0316_woSSOabs
Taxo_split=Taxo_0316_woSSOabs_split
Env=rawMetadata_ENV
c(nrow(OTU),nrow(Taxo),nrow(Taxo_split),nrow(Env)) #check if dimension fit


###calculations###  
OTU_subset_raw=OTU[,Index] #make subset #73 samples
OTU_subset=OTU_subset_raw[apply(decostand(OTU_subset_raw,"pa"),1,sum)!=0,] #remove empty OTU #1342 OTU
OTU_subsetrel=decostand(OTU_subset,"total",MARGIN=2) #make relative
apply(OTU_subsetrel,2,sum) #check if samples =1

##Taxonomy
Taxo_subset=as.data.frame(Taxo[apply(decostand(OTU_subset_raw,"pa"),1,sum)!=0,1])
Taxo_subset=droplevels(Taxo_subset)
Taxo_subset_split=as.data.frame((Taxo_split[apply(decostand(OTU_subset_raw,"pa"),1,sum)!=0,]))

##clean SSOabs
OTU_subsetrel_woSSOabs=as.data.frame(OTU_subsetrel[apply(OTU_subset,1,sum)>1,])
Taxo_subset_woSSOabs=as.data.frame(Taxo_subset[apply(OTU_subset,1,sum)>1,1])
Taxo_subset_split_woSSOabs=as.data.frame(Taxo_subset_split[apply(OTU_subset,1,sum)>1,])
OTU_subset_woSSOabs=as.data.frame(OTU_subset[apply(OTU_subset,1,sum)>1,])
###1272 Taxa/ OTU in full dataset 

####ENv data
cleanENV=Env[Index,]


###cleaning####
c(nrow(OTU_subsetrel_woSSOabs),nrow(Taxo_subset_woSSOabs),
  nrow(Taxo_subset_woSSOabs),nrow(OTU_subset_woSSOabs), nrow(cleanENV))#check dim

OTU_NS_woSSOabs=OTU_subset_woSSOabs
OTU_NS_woSSOabs_rel=OTU_subsetrel_woSSOabs
Taxo_NS_woSSOabs=Taxo_subset_woSSOabs
Taxo_NS_split_woSSOabs=Taxo_subset_split_woSSOabs
Metadata_Env_NS=cleanENV
rm(OTU,OTU_subset, OTU_subset_raw,Taxo_subset,Taxo_subset_split,
   Taxo,Taxo_split,OTU_subsetrel_woSSOabs,Taxo_subset_woSSOabs,Taxo_subset_split_woSSOabs
   ,cleanENV, Env, OTU_subset_woSSOabs) #remove matrices for calculations

###############
##make averages per year for the NS dataset
#averages per year
OTU_NS_average_woSSOabs = Merging_ALk(t(rbind(OTU_NS_woSSOabs,Metadata_Env_NS$Year)), k=1)
OTU_NS_average_woSSOabs_rel=decostand(OTU_NS_average_woSSOabs,"total",MARGIN=1) #make relative

#############################################
#2000-3000 m waterdepth stations (latitudinal transect + HGIII)

##subset with only stations from the latitudinal transect
Index=  which(Metadata_Env_0316$Elevation_of_event< -1999 & Metadata_Env_0316$Elevation_of_event> -3001 & Metadata_Env_0316$Optional_label<30) ###list with samples to be used

OTU=OTU_0316_woSSOabs
Taxo=Taxo_0316_woSSOabs
Taxo_split=Taxo_0316_woSSOabs_split
Env=rawMetadata_ENV
c(nrow(OTU),nrow(Taxo),nrow(Taxo_split),nrow(Env)) #check if dimension fit


###calculations###  
OTU_subset_raw=OTU[,Index] #make subset #73 samples
OTU_subset=OTU_subset_raw[apply(decostand(OTU_subset_raw,"pa"),1,sum)!=0,] #remove empty OTU #1342 OTU
OTU_subsetrel=decostand(OTU_subset,"total",MARGIN=2) #make relative
apply(OTU_subsetrel,2,sum) #check if samples =1

##Taxonomy
Taxo_subset=as.data.frame(Taxo[apply(decostand(OTU_subset_raw,"pa"),1,sum)!=0,1])
Taxo_subset=droplevels(Taxo_subset)
Taxo_subset_split=as.data.frame((Taxo_split[apply(decostand(OTU_subset_raw,"pa"),1,sum)!=0,]))

##clean SSOabs
OTU_subsetrel_woSSOabs=as.data.frame(OTU_subsetrel[apply(OTU_subset,1,sum)>1,])
Taxo_subset_woSSOabs=as.data.frame(Taxo_subset[apply(OTU_subset,1,sum)>1,1])
Taxo_subset_split_woSSOabs=as.data.frame(Taxo_subset_split[apply(OTU_subset,1,sum)>1,])
OTU_subset_woSSOabs=as.data.frame(OTU_subset[apply(OTU_subset,1,sum)>1,])
###1166Taxa/ OTU in full dataset 

####ENv data
cleanENV=Env[Index,]


###cleaning####
c(nrow(OTU_subsetrel_woSSOabs),nrow(Taxo_subset_woSSOabs),
  nrow(Taxo_subset_woSSOabs),nrow(OTU_subset_woSSOabs), nrow(cleanENV))#check dim

OTU_2000to3000_woSSOabs=OTU_subset_woSSOabs
OTU_2000to3000_woSSOabs_rel=OTU_subsetrel_woSSOabs
Taxo_2000to3000_woSSOabs=Taxo_subset_woSSOabs
Taxo_2000to3000_split_woSSOabs=Taxo_subset_split_woSSOabs
Metadata_Env_2000to3000=cleanENV
rm(OTU,OTU_subset, OTU_subset_raw,Taxo_subset,Taxo_subset_split,
   Taxo,Taxo_split,OTU_subsetrel_woSSOabs,Taxo_subset_woSSOabs,Taxo_subset_split_woSSOabs
   ,cleanENV, Env, OTU_subset_woSSOabs) #remove matrices for calculations



###################
##subset for comparison of Hausgarten and EG stations
##
Index=  c(119,120,121,129,130,131,132) ###list with samples to be used
OTU=OTU_0316_woSSOabs
Taxo=Taxo_0316_woSSOabs
Taxo_split=Taxo_0316_woSSOabs_split
Env=rawMetadata_ENV
c(nrow(OTU),nrow(Taxo),nrow(Taxo_split)) #check if dimension fit#1394 rows


###calculations###  
OTU_subset_raw=OTU[,Index] #make subset #21 samples
OTU_subset=OTU_subset_raw[apply(OTU_subset_raw,1,sum)!=0,] #remove empty OTU #1036 OTU
OTU_subsetrel=decostand(OTU_subset,"total",MARGIN=2) #make relative
apply(OTU_subsetrel,2,sum) #check if samples =1

##Taxonomy
Taxo_subset=as.data.frame(Taxo[apply(OTU_subset_raw,1,sum)!=0,1])
Taxo_subset=droplevels(Taxo_subset)
Taxo_subset_split=as.data.frame((Taxo_split[apply(OTU_subset_raw,1,sum)!=0,]))

##clean SSOabs
OTU_subsetrel_woSSOabs=as.data.frame(OTU_subsetrel[apply(OTU_subset,1,sum)>1,])
Taxo_subset_woSSOabs=as.data.frame(Taxo_subset[apply(OTU_subset,1,sum)>1,1])
Taxo_subset_split_woSSOabs=as.data.frame(Taxo_subset_split[apply(OTU_subset,1,sum)>1,])
OTU_subset_woSSOabs=as.data.frame(OTU_subset[apply(OTU_subset,1,sum)>1,])
###909 Taxa/ OTU in full dataset 

####ENv data
cleanENV=Env[Index,]


###cleaning####
c(nrow(OTU_subset_woSSOabs),nrow(OTU_subsetrel_woSSOabs),nrow(Taxo_subset_woSSOabs),nrow(Taxo_subset_split_woSSOabs))#check dim
OTU_HGEG16_woSSOabs=OTU_subset_woSSOabs
OTU_HGEG16_woSSOabs_rel=OTU_subsetrel_woSSOabs
Taxo_HGEG16_woSSOabs=Taxo_subset_woSSOabs
Taxo_HGEG16_woSSOabs_split=Taxo_subset_split_woSSOabs
Metadata_Env_HGEG16=cleanENV
rm(OTU,OTU_subset, OTU_subset_raw,Taxo_subset,Taxo_subset_split,
   Taxo,Taxo_split,OTU_subsetrel_woSSOabs,Taxo_subset_woSSOabs,Taxo_subset_split_woSSOabs
   , OTU_subset_woSSOabs,cleanENV) #remove matrices for calculations

####
###################
##subset for comparison of Hausgarten and EG stations
##
Index=  c(98,99,100,104,105,106,107,119,120,121,129,130,131,132) ###list with samples to be used
OTU=OTU_0316_woSSOabs
Taxo=Taxo_0316_woSSOabs
Taxo_split=Taxo_0316_woSSOabs_split
Env=rawMetadata_ENV
c(nrow(OTU),nrow(Taxo),nrow(Taxo_split)) #check if dimension fit#1394 rows


###calculations###  
OTU_subset_raw=OTU[,Index] #make subset #21 samples
OTU_subset=OTU_subset_raw[apply(OTU_subset_raw,1,sum)!=0,] #remove empty OTU #1036 OTU
OTU_subsetrel=decostand(OTU_subset,"total",MARGIN=2) #make relative
apply(OTU_subsetrel,2,sum) #check if samples =1

##Taxonomy
Taxo_subset=as.data.frame(Taxo[apply(OTU_subset_raw,1,sum)!=0,1])
Taxo_subset=droplevels(Taxo_subset)
Taxo_subset_split=as.data.frame((Taxo_split[apply(OTU_subset_raw,1,sum)!=0,]))

##clean SSOabs
OTU_subsetrel_woSSOabs=as.data.frame(OTU_subsetrel[apply(OTU_subset,1,sum)>1,])
Taxo_subset_woSSOabs=as.data.frame(Taxo_subset[apply(OTU_subset,1,sum)>1,1])
Taxo_subset_split_woSSOabs=as.data.frame(Taxo_subset_split[apply(OTU_subset,1,sum)>1,])
OTU_subset_woSSOabs=as.data.frame(OTU_subset[apply(OTU_subset,1,sum)>1,])
###909 Taxa/ OTU in full dataset 

####ENv data
cleanENV=Env[Index,]


###cleaning####
c(nrow(OTU_subset_woSSOabs),nrow(OTU_subsetrel_woSSOabs),nrow(Taxo_subset_woSSOabs),nrow(Taxo_subset_split_woSSOabs))#check dim
OTU_HGEG1416_woSSOabs=OTU_subset_woSSOabs
OTU_HGEG1416_woSSOabs_rel=OTU_subsetrel_woSSOabs
Taxo_HGEG1416_woSSOabs=Taxo_subset_woSSOabs
Taxo_HGEG1416_woSSOabs_split=Taxo_subset_split_woSSOabs
Metadata_Env_HGEG1416=cleanENV
rm(OTU,OTU_subset, OTU_subset_raw,Taxo_subset,Taxo_subset_split,
   Taxo,Taxo_split,OTU_subsetrel_woSSOabs,Taxo_subset_woSSOabs,Taxo_subset_split_woSSOabs
   , OTU_subset_woSSOabs,cleanENV) #remove matrices for calculations




#####################################
#normalize, transofrm env data???
