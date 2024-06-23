# missMethyl with C7: immunologic signature gene sets
# https://bioconductor.org/packages/devel/bioc/vignettes/missMethyl/inst/doc/missMethyl.html
# Created by Mannan Luo (l.mannan@erasmusmc.nl)

rm(list=ls())
library(missMethyl)  
library(data.table)

## Gene sets that represent cell states and perturbations within the immune system, 
# i.e., C7: immunologic signature gene sets, from the Molecular Signatures Database (MSigDB), were combined with the R package missMethyl to identify enriched immunologic signatures for each cell state comparison. 

# get C7 file
C7 <- readRDS(url("https://bioinf.wehi.edu.au/MSigDB/v7.1/Hs.c7.all.v7.1.entrez.rds"))

####################################################################
#  load EWAS results _Total score
ewas <- fread("~/revision/merged_meta_analyzed_results_total_infection_model1.csv")
head(ewas)
str(ewas)

# Get all the CpG sites used in the analysis to form the background (extensive filtering of the CpGs has been performed prior to analysi)
all <- ewas$cpg
# Total number of CpG sites tested
length(all)
#393360

## Get hits based on p threshold
sigCpGs <- ewas$cpg[ewas$meta_pval<5e-05]
length(sigCpGs)
#35

# run the function with C7 database
gsa.C7 <- gsameth(sig.cpg=sigCpGs, all.cpg=all, collection=C7, sig.genes = TRUE)
topGSA(gsa.C7)
#SigGenesInSet GSE14769_UNSTIM_VS_360MIN_LPS_BMDM_DN   AMBRA1,TAPBP,VOPP1,CTU1,PJA2
table(gsa.C7$FDR<0.05)

top20<-topGSA(gsa.C7, n=20)
setwd("~/EWAS_Infection_Shadow/")
fwrite(top20, file="C7_TOP20_5E-5_total_Mod1.csv", row.names=T)

rm(ewas,all,sigCpGs, gsa.C7,top20)

####################################################################
#  EWAS_trimester1
ewas <- fread("~/merged_meta_analyzed_results_trimester1_infection_model1.csv")
str(ewas)

# Get all the CpG sites used in the analysis to form the background (extensive filtering of the CpGs has been performed prior to analysi)
all <- ewas$cpg
# Total number of CpG sites tested
length(all)
#393360

## Get hits based on p threshold
sigCpGs <- ewas$cpg[ewas$meta_pval<5e-05]
length(sigCpGs)
#21

# run the function with C7 database
gsa.C7 <- gsameth(sig.cpg=sigCpGs, all.cpg=all, collection=C7, sig.genes = TRUE)
topGSA(gsa.C7)
table(gsa.C7$FDR<0.05)
top20<-topGSA(gsa.C7, n=20)

setwd("~/EWAS_Infection_Shadow/")
fwrite(top20, file="C7_TOP20_5E-5_trimester1_Mod1.csv", row.names=T)

rm(ewas,all,sigCpGs, gsa.C7,top20)

####################################################################
# EWAS_trimester2
ewas <- fread("~/merged_meta_analyzed_results_trimester2_infection_model1.csv")
str(ewas)

# Get all the CpG sites used in the analysis to form the background (extensive filtering of the CpGs has been performed prior to analysi)
all <- ewas$cpg
# Total number of CpG sites tested
length(all)
#393360

## Get hits based on p threshold
sigCpGs <- ewas$cpg[ewas$meta_pval<5e-05]
length(sigCpGs)

# run the function with C7 database
gsa.C7 <- gsameth(sig.cpg=sigCpGs, all.cpg=all, collection=C7, sig.genes = TRUE)
topGSA(gsa.C7)
table(gsa.C7$FDR<0.05)
top20<-topGSA(gsa.C7, n=20)

setwd("~/EWAS_Infection_Shadow/")
fwrite(top20, file="C7_TOP20_5E-5_trimester2_Mod1.csv", row.names=T)

rm(ewas,all,sigCpGs,gsa.C7,top20)

####################################################################
# EWAS_trimester3
ewas <- fread("~/merged_meta_analyzed_results_trimester3_infection_model1.csv")
str(ewas)

# Get all the CpG sites used in the analysis to form the background (extensive filtering of the CpGs has been performed prior to analysi)
all <- ewas$cpg
# Total number of CpG sites tested
length(all)
#393360

## Get hits based on p threshold
sigCpGs <- ewas$cpg[ewas$meta_pval<5e-05]
length(sigCpGs)

# run the function with C7 database
gsa.C7 <- gsameth(sig.cpg=sigCpGs, all.cpg=all, collection=C7, sig.genes = TRUE)
topGSA(gsa.C7)
table(gsa.C7$FDR<0.05)
top20<-topGSA(gsa.C7, n=20)

setwd("path_to_store_results")
fwrite(top20, file="C7_TOP20_5E-5_trimester3_Mod1.csv", row.names=T)

rm(list=ls())

#\ END OF SCRIPT
