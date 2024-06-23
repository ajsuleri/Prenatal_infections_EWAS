#########################################################################
#####################PROJECT: PRENATAL INFECTION EWAS####################
#########################################################################

#Author: Anna Suleri

# Goal of this script: 
# 1. We create a dataframe of CpGs that overlap in 450k and in epic for both 450k and epic
# 2. We create a dataframe for cpgs in epic that are not in the 450k for an extra exploratory analysis
# 3. We run the following EWASES on each dataframe:
######a) total infection - CpG at birth
######b) trimester 1 infection - CpG at birth
######c) trimester 2 infection - CpG at birth
######d) trimester 3 infection - CpG at birth
#4. We run these four EWASES for the following two models:
######1. Model 1: child sex, maternal age at delivery, maternal education, maternal smoking, parity, batch effects, cel type proportions
######2. Model 2: model 1 + additionally gestational age at birth and birth weight 
#4.We check the diagnostics of each EWAS with a Q-Q plot of the pvalues, lambda inflation factor and make a manhattan plot of potential hits  
#5.Then we pool all the results for these 4x EWAS with each 2 models (so we end up with 8 EWAS for 450k and 8 EWAS for EPIC) --> we apply inverse weighted fixed-effects meta analysis to pool the results 
#6.We check the quality of the ewas 
#7.We apply post processing steps such as excluding cross reactive probes and flagging probes containing snps
#8.We create a manhattan plot 
#9. We run a posthoc analysis
#10. We visualize top hits 

### Set working directory and load packages ###
library(foreign) # to load data 
library(haven) # to load data 
library(dplyr) # to manipulate data 
library(pbmcapply) # to run EWAS with parallel processing 

setwd('set_path_to_home_directory')

## source functions file
source('0.Functions_script.R')

## run this script on the GENA server 

###----------------STEP 1: Prepare dataframes----------------###
## load methylation data, first for 450k and then for epic
load("path_450K_methylation_data") # beta matrix is called x
load("path_EPIC_methylation_data") 
epic_df <- GENR_EPICv1METH_Norm_Betas_birth_3IQRwinsorized.data # rename beta matrix of epic 

# check dimensions to see if rownames and colnames are id and cpgs (so identitical for both dataframes, if not, make it identitical)
dim(x) #rows = 458563, columns = 1396
dim(epic_df) #rows = 808183, columns = 1115

# create df of cpgs that are available in both 450k assay and epic assay 
overlapping_betas <- as.data.frame(x[rownames(x) %in% rownames(epic_df),]) #393360 cpgs 
overlapping_betas2 <- as.data.frame(epic_df[rownames(epic_df) %in% rownames(x),])

# create df of cpgs that are only available in epic assay 
epic_only_betas <- as.data.frame(epic_df[!rownames(epic_df) %in% rownames(x),]) #414823 cpgs 

# transpose methylation dataframes so that cpgs are columns and participants are rows
overlap_df <- as.data.frame(t(overlapping_betas))
dim(overlap_df) #check if it went correctly 
overlap_df2 <- as.data.frame(t(overlapping_betas2))
dim(overlap_df2) #check if it went correctly 
epic_only_df <- as.data.frame(t(epic_only_betas)) 
dim(epic_only_df) #check if it went correctly

## load in phenotypic data
df_phenotype_450k <- readRDS('df_final_450k.rds') #Sample_ID
df_phenotype_epic <- readRDS('df_final_epic.rds') #SampleID

# check str of variables and recode if necessary 
str(df_phenotype_450k)
str(df_phenotype_epic)

df_phenotype_450k$Sample_Plate <- as.factor(df_phenotype_450k$Sample_Plate)
df_phenotype_epic$Sample_Plate <- as.factor(df_phenotype_epic$Sample_Plate)

## match methylation data to phenotypic data
overlap_df_450k <- overlap_df[as.character(df_phenotype_450k$Sample_ID),]
table(rownames(overlap_df_450k)) == df_phenotype_450k$Sample_ID

overlap_df_epic <- overlap_df2[as.character(df_phenotype_epic$SampleID),]
table(rownames(overlap_df_epic)) == df_phenotype_epic$SampleID

epic_only_betas_df <- epic_only_df[as.character(df_phenotype_epic$SampleID),]
table(rownames(epic_only_betas_df)) == df_phenotype_epic$SampleID

# save datasets to later use also in regional analyses 
saveRDS(overlap_df_450k, 'overlap_df_450k.rds') 
saveRDS(overlap_df_epic, 'overlap_df_epic.rds') 

#\

###----------------STEP 2: RUN EWAS----------------###
## functions to run ewas models for total infection score
#function for model 1, 450k
regress1 <- function(cpg){
  cpg_phenotype_data <- cbind(cpg, df_phenotype_450k)
  coef(summary(lm(cpg ~ sumscore_inf_tot_standardized + GENDER + SMOKE_ALL + AGE_M_v2 + PARITY + EDUCM_3groups + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC, data = cpg_phenotype_data)))[2,]
}

#function for model 2, 450k
regress2 <- function(cpg){
  cpg_phenotype_data <- cbind(cpg, df_phenotype_450k)
  coef(summary(lm(cpg ~ sumscore_inf_tot_standardized + GENDER + SMOKE_ALL + AGE_M_v2 + PARITY + EDUCM_3groups + GESTBIR + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC, data = cpg_phenotype_data)))[2,]
}

#function for model 1, epic 
regress3 <- function(cpg){
  cpg_phenotype_data <- cbind(cpg, df_phenotype_epic)
  coef(summary(lm(cpg ~ sumscore_inf_tot_standardized + GENDER + SMOKE_ALL + AGE_M_v2 + PARITY + EDUCM_3groups + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC, data = cpg_phenotype_data)))[2,]
}

#function for model 2, epic
regress4 <- function(cpg){
  cpg_phenotype_data <- cbind(cpg, df_phenotype_epic)
  coef(summary(lm(cpg ~ sumscore_inf_tot_standardized + GENDER + SMOKE_ALL + AGE_M_v2 + PARITY + EDUCM_3groups + GESTBIR + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC, data = cpg_phenotype_data)))[2,]
}

## functions for trimester specific infection scores

#function for trimester 1, model 1, 450k
regress5 <- function(cpg){
  cpg_phenotype_data <- cbind(cpg, df_phenotype_450k)
  coef(summary(lm(cpg ~ sumscore_inf_tri1_standardized + GENDER + SMOKE_ALL + AGE_M_v2 + PARITY + EDUCM_3groups + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC, data = cpg_phenotype_data)))[2,]
}

#function for trimester 1, model 2, 450k
regress6 <- function(cpg){
  cpg_phenotype_data <- cbind(cpg, df_phenotype_450k)
  coef(summary(lm(cpg ~ sumscore_inf_tri1_standardized + GENDER + SMOKE_ALL + AGE_M_v2 + PARITY + EDUCM_3groups + GESTBIR + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC, data = cpg_phenotype_data)))[2,]
}

#function for trimester 2, model 1, 450k
regress7 <- function(cpg){
  cpg_phenotype_data <- cbind(cpg, df_phenotype_450k)
  coef(summary(lm(cpg ~ sumscore_inf_tri2_standardized + GENDER + SMOKE_ALL + AGE_M_v2 + PARITY + EDUCM_3groups + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC, data = cpg_phenotype_data)))[2,]
}

#function for trimester 2, model 2, 450k
regress8 <- function(cpg){
  cpg_phenotype_data <- cbind(cpg, df_phenotype_450k)
  coef(summary(lm(cpg ~ sumscore_inf_tri2_standardized + GENDER + SMOKE_ALL + AGE_M_v2 + PARITY + EDUCM_3groups + GESTBIR + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC, data = cpg_phenotype_data)))[2,]
}

#function for trimester 3, model 1, 450k
regress9 <- function(cpg){
  cpg_phenotype_data <- cbind(cpg, df_phenotype_450k)
  coef(summary(lm(cpg ~ sumscore_inf_tri3_standardized + GENDER + SMOKE_ALL + AGE_M_v2 + PARITY + EDUCM_3groups + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC, data = cpg_phenotype_data)))[2,]
}

#function for trimester 3, model 2, 450k
regress10 <- function(cpg){
  cpg_phenotype_data <- cbind(cpg, df_phenotype_450k)
  coef(summary(lm(cpg ~ sumscore_inf_tri3_standardized + GENDER + SMOKE_ALL + AGE_M_v2 + PARITY + EDUCM_3groups + GESTBIR + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC, data = cpg_phenotype_data)))[2,]
}

#function for trimester 1, model 1, epic
regress11 <- function(cpg){
  cpg_phenotype_data <- cbind(cpg, df_phenotype_epic)
  coef(summary(lm(cpg ~ sumscore_inf_tri1_standardized + GENDER + SMOKE_ALL + AGE_M_v2 + PARITY + EDUCM_3groups + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC, data = cpg_phenotype_data)))[2,]
}

#function for trimester 1, model 2, epic
regress12 <- function(cpg){
  cpg_phenotype_data <- cbind(cpg, df_phenotype_epic)
  coef(summary(lm(cpg ~ sumscore_inf_tri1_standardized + GENDER + SMOKE_ALL + AGE_M_v2 + PARITY + EDUCM_3groups + GESTBIR + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC, data = cpg_phenotype_data)))[2,]
}

#function for trimester 2, model 1, epic
regress13 <- function(cpg){
  cpg_phenotype_data <- cbind(cpg, df_phenotype_epic)
  coef(summary(lm(cpg ~ sumscore_inf_tri2_standardized + GENDER + SMOKE_ALL + AGE_M_v2 + PARITY + EDUCM_3groups + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC, data = cpg_phenotype_data)))[2,]
}

#function for trimester 2, model 2, epic
regress14 <- function(cpg){
  cpg_phenotype_data <- cbind(cpg, df_phenotype_epic)
  coef(summary(lm(cpg ~ sumscore_inf_tri2_standardized + GENDER + SMOKE_ALL + AGE_M_v2 + PARITY + EDUCM_3groups + GESTBIR + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC, data = cpg_phenotype_data)))[2,]
}

#function for trimester 3, model 1, epic
regress15 <- function(cpg){
  cpg_phenotype_data <- cbind(cpg, df_phenotype_epic)
  coef(summary(lm(cpg ~ sumscore_inf_tri3_standardized + GENDER + SMOKE_ALL + AGE_M_v2 + PARITY + EDUCM_3groups + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC, data = cpg_phenotype_data)))[2,]
}

#function for trimester 3, model 2, epic
regress16 <- function(cpg){
  cpg_phenotype_data <- cbind(cpg, df_phenotype_epic)
  coef(summary(lm(cpg ~ sumscore_inf_tri3_standardized + GENDER + SMOKE_ALL + AGE_M_v2 + PARITY + EDUCM_3groups + GESTBIR + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC, data = cpg_phenotype_data)))[2,]
}

## run all the ewases and save results 

# first for 450k
#ewas 1: total infection - 450k data - model 1
overlap_ewas_450k_model1.list <- pbmclapply(overlap_df_450k, regress1, mc.cores = 12)
overlap_ewas_450k_model1.data <- as.data.frame(do.call(rbind, overlap_ewas_450k_model1.list))
overlap_ewas_450k_model1.data$cpg <- names(overlap_df_450k)
save(overlap_ewas_450k_model1.data, file = 'overlap_ewas_450k_model1.rdata')

#ewas 2: total infection - 450k data - model 2
overlap_ewas_450k_model2.list <- pbmclapply(overlap_df_450k, regress2, mc.cores = 12)
overlap_ewas_450k_model2.data <- as.data.frame(do.call(rbind, overlap_ewas_450k_model2.list))
overlap_ewas_450k_model2.data$cpg <- names(overlap_df_450k)
save(overlap_ewas_450k_model2.data, file = 'overlap_ewas_450k_model2.rdata')

#ewas 3: trimester 1 - 450k data - model 1
overlap_ewas_450k_trimester1_model1.list <- pbmclapply(overlap_df_450k, regress5, mc.cores = 12)
overlap_ewas_450k_trimester1_model1.data <- as.data.frame(do.call(rbind, overlap_ewas_450k_trimester1_model1.list))
overlap_ewas_450k_trimester1_model1.data$cpg <- names(overlap_df_450k)
save(overlap_ewas_450k_trimester1_model1.data, file = 'overlap_ewas_450k_trimester1_model1.rdata')

#ewas 4: trimester 1 - 450k data - model 2
overlap_ewas_450k_trimester1_model2.list <- pbmclapply(overlap_df_450k, regress6, mc.cores = 12)
overlap_ewas_450k_trimester1_model2.data <- as.data.frame(do.call(rbind, overlap_ewas_450k_trimester1_model2.list))
overlap_ewas_450k_trimester1_model2.data$cpg <- names(overlap_df_450k)
save(overlap_ewas_450k_trimester1_model2.data, file = 'overlap_ewas_450k_trimester1_model2.rdata')

#ewas 5: trimester 2 - 450k data - model 1
overlap_ewas_450k_trimester2_model1.list <- pbmclapply(overlap_df_450k,regress7, mc.cores = 12)
overlap_ewas_450k_trimester2_model1.data <- as.data.frame(do.call(rbind, overlap_ewas_450k_trimester2_model1.list))
overlap_ewas_450k_trimester2_model1.data$cpg <- names(overlap_df_450k)
save(overlap_ewas_450k_trimester2_model1.data, file = 'overlap_ewas_450k_trimester2_model1.rdata')

#ewas 6: trimester 2 - 450k data - model 2
overlap_ewas_450k_trimester2_model2.list <- pbmclapply(overlap_df_450k,regress8, mc.cores = 12)
overlap_ewas_450k_trimester2_model2.data <- as.data.frame(do.call(rbind, overlap_ewas_450k_trimester2_model2.list))
overlap_ewas_450k_trimester2_model2.data$cpg <- names(overlap_df_450k)
save(overlap_ewas_450k_trimester2_model2.data, file = 'overlap_ewas_450k_trimester2_model2.rdata')

#ewas 7: trimester 3 - 450k data - model 1
overlap_ewas_450k_trimester3_model1.list <- pbmclapply(overlap_df_450k,regress9, mc.cores = 12)
overlap_ewas_450k_trimester3_model1.data <- as.data.frame(do.call(rbind, overlap_ewas_450k_trimester3_model1.list))
overlap_ewas_450k_trimester3_model1.data$cpg <- names(overlap_df_450k)
save(overlap_ewas_450k_trimester3_model1.data, file = 'overlap_ewas_450k_trimester3_model1.rdata')

#ewas 8: trimester 3 - 450k data - model 2
overlap_ewas_450k_trimester3_model2.list <- pbmclapply(overlap_df_450k,regress10, mc.cores = 12)
overlap_ewas_450k_trimester3_model2.data <- as.data.frame(do.call(rbind, overlap_ewas_450k_trimester3_model2.list))
overlap_ewas_450k_trimester3_model2.data$cpg <- names(overlap_df_450k)
save(overlap_ewas_450k_trimester3_model2.data, file = 'overlap_ewas_450k_trimester3_model2.rdata')

# then for epic

#ewas 9: total infection - epic data - model 1
overlap_ewas_epic_model1.list <- pbmclapply(overlap_df_epic, regress3, mc.cores = 12)
overlap_ewas_epic_model1.data <- as.data.frame(do.call(rbind, overlap_ewas_epic_model1.list))
overlap_ewas_epic_model1.data$cpg <- names(overlap_df_epic)
save(overlap_ewas_epic_model1.data, file = 'overlap_ewas_epic_model1.rdata')

#ewas 10: total infection - epic data- model 2
overlap_ewas_epic_model2.list <- pbmclapply(overlap_df_epic, regress4, mc.cores = 12)
overlap_ewas_epic_model2.data <- as.data.frame(do.call(rbind, overlap_ewas_epic_model2.list))
overlap_ewas_epic_model2.data$cpg <- names(overlap_df_epic)
save(overlap_ewas_epic_model2.data, file = 'overlap_ewas_epic_model2.rdata')

#ewas 11: total infection - epic cpgs only - model 1
epic_only_ewas_model1.list <- pbmclapply(epic_only_betas_df, regress3, mc.cores = 12)
epic_only_ewas_model1.data <- as.data.frame(do.call(rbind, epic_only_ewas_model1.list))
epic_only_ewas_model1.data$cpg <- names(epic_only_betas_df)
save(epic_only_ewas_model1.data, file = 'epic_only_ewas_model1.rdata')

#ewas 12: total infection - epic cpgs only - model 2
epic_only_ewas_model2.list <- pbmclapply(epic_only_betas_df, regress4, mc.cores = 12)
epic_only_ewas_model2.data <- as.data.frame(do.call(rbind, epic_only_ewas_model2.list))
epic_only_ewas_model2.data$cpg <- names(epic_only_betas_df)
save(epic_only_ewas_model2.data, file = 'epic_only_ewas_model2.rdata')

#ewas 13: trimester 1 - epic data - model 1
overlap_ewas_epic_trimester1_model1.list <- pbmclapply(overlap_df_epic, regress11, mc.cores = 12)
overlap_ewas_epic_trimester1_model1.data <- as.data.frame(do.call(rbind, overlap_ewas_epic_trimester1_model1.list))
overlap_ewas_epic_trimester1_model1.data$cpg <- names(overlap_df_epic)
save(overlap_ewas_epic_trimester1_model1.data, file = 'overlap_ewas_epic_trimester1_model1.rdata')

#ewas 14: trimester 1 - epic data - model 2
overlap_ewas_epic_trimester1_model2.list <- pbmclapply(overlap_df_epic,regress12, mc.cores = 12)
overlap_ewas_epic_trimester1_model2.data <- as.data.frame(do.call(rbind, overlap_ewas_epic_trimester1_model2.list))
overlap_ewas_epic_trimester1_model2.data$cpg <- names(overlap_df_epic)
save(overlap_ewas_epic_trimester1_model2.data, file = 'overlap_ewas_epic_trimester1_model2.rdata')

#ewas 15: trimester 2 - epic data - model 1
overlap_ewas_epic_trimester2_model1.list <- pbmclapply(overlap_df_epic, regress13, mc.cores = 12)
overlap_ewas_epic_trimester2_model1.data <- as.data.frame(do.call(rbind, overlap_ewas_epic_trimester2_model1.list))
overlap_ewas_epic_trimester2_model1.data$cpg <- names(overlap_df_epic)
save(overlap_ewas_epic_trimester2_model1.data, file = 'overlap_ewas_epic_trimester2_model1.rdata')

#ewas 16: trimester 2 - epic data - model 2
overlap_ewas_epic_trimester2_model2.list <- pbmclapply(overlap_df_epic, regress14, mc.cores = 12)
overlap_ewas_epic_trimester2_model2.data <- as.data.frame(do.call(rbind, overlap_ewas_epic_trimester2_model2.list))
overlap_ewas_epic_trimester2_model2.data$cpg <- names(overlap_df_epic)
save(overlap_ewas_epic_trimester2_model2.data, file = 'overlap_ewas_epic_trimester2_model2.rdata')

#ewas 17: trimester 3 - epic data - model 1
overlap_ewas_epic_trimester3_model1.list <- pbmclapply(overlap_df_epic, regress15, mc.cores = 12)
overlap_ewas_epic_trimester3_model1.data <- as.data.frame(do.call(rbind, overlap_ewas_epic_trimester3_model1.list))
overlap_ewas_epic_trimester3_model1.data$cpg <- names(overlap_df_epic)
save(overlap_ewas_epic_trimester3_model1.data, file = 'overlap_ewas_epic_trimester3_model1.rdata')

#ewas 18: trimester 3 - epic data - model 2
overlap_ewas_epic_trimester3_model2.list <- pbmclapply(overlap_df_epic, regress16, mc.cores = 12)
overlap_ewas_epic_trimester3_model2.data <- as.data.frame(do.call(rbind, overlap_ewas_epic_trimester3_model2.list))
overlap_ewas_epic_trimester3_model2.data$cpg <- names(overlap_df_epic)
save(overlap_ewas_epic_trimester3_model2.data, file = 'overlap_ewas_epic_trimester3_model2.rdata')

#ewas 19: trimester 1 - epic cpgs only - model 1
epic_only_ewas_epic_trimester1_model1.list <- pbmclapply(epic_only_betas_df, regress11, mc.cores = 12)
epic_only_ewas_epic_trimester1_model1.data <- as.data.frame(do.call(rbind, epic_only_ewas_epic_trimester1_model1.list))
epic_only_ewas_epic_trimester1_model1.data$cpg <- names(epic_only_betas_df)
save(epic_only_ewas_epic_trimester1_model1.data, file = 'epic_only_ewas_epic_trimester1_model1.rdata')

#ewas 20: trimester 1 - epic cpgs only - model 2
epic_only_ewas_epic_trimester1_model2.list <- pbmclapply(epic_only_betas_df, regress12, mc.cores = 12)
epic_only_ewas_epic_trimester1_model2.data <- as.data.frame(do.call(rbind, epic_only_ewas_epic_trimester1_model2.list))
epic_only_ewas_epic_trimester1_model2.data$cpg <- names(epic_only_betas_df)
save(epic_only_ewas_epic_trimester1_model2.data, file = 'epic_only_ewas_epic_trimester1_model2.rdata')

#ewas 21: trimester 2 - epic cpgs only - model 1
epic_only_ewas_epic_trimester2_model1.list <- pbmclapply(epic_only_betas_df, regress13, mc.cores = 12)
epic_only_ewas_epic_trimester2_model1.data <- as.data.frame(do.call(rbind, epic_only_ewas_epic_trimester2_model1.list))
epic_only_ewas_epic_trimester2_model1.data$cpg <- names(epic_only_betas_df)
save(epic_only_ewas_epic_trimester2_model1.data, file = 'epic_only_ewas_epic_trimester2_model1.rdata')

#ewas 22: trimester 2 - epic cpgs only - model 2
epic_only_ewas_epic_trimester2_model2.list <- pbmclapply(epic_only_betas_df, regress14, mc.cores = 12)
epic_only_ewas_epic_trimester2_model2.data <- as.data.frame(do.call(rbind, epic_only_ewas_epic_trimester2_model2.list))
epic_only_ewas_epic_trimester2_model2.data$cpg <- names(epic_only_betas_df)
save(epic_only_ewas_epic_trimester2_model2.data, file = 'epic_only_ewas_epic_trimester2_model2.rdata')

#ewas 23: trimester 3 - epic cpgs only - model 1
epic_only_ewas_epic_trimester3_model1.list <- pbmclapply(epic_only_betas_df, regress15, mc.cores = 12)
epic_only_ewas_epic_trimester3_model1.data <- as.data.frame(do.call(rbind, epic_only_ewas_epic_trimester3_model1.list))
epic_only_ewas_epic_trimester3_model1.data$cpg <- names(epic_only_betas_df)
save(epic_only_ewas_epic_trimester3_model1.data, file = 'epic_only_ewas_epic_trimester3_model1.rdata')

#ewas 24: trimester 3 - epic cpgs only - model 2
epic_only_ewas_epic_trimester3_model2.list <- pbmclapply(epic_only_betas_df, regress16, mc.cores = 12)
epic_only_ewas_epic_trimester3_model2.data <- as.data.frame(do.call(rbind, epic_only_ewas_epic_trimester3_model2.list))
epic_only_ewas_epic_trimester3_model2.data$cpg <- names(epic_only_betas_df)
save(epic_only_ewas_epic_trimester3_model2.data, file = 'epic_only_ewas_epic_trimester3_model2.rdata')

#\

###----------------STEP 3: Quality checks EWAS----------------###

## QC of EWAS 

# Set the path to the folder containing .rdata files
dataframes_folder <- 'path_to_data_folder'

# List all .rdata files in the folder
rdata_files <- list.files(dataframes_folder, pattern = "\\.rdata$", full.names = TRUE)

# Load each dataframe and assign it to its original name
for (rdata_file in rdata_files) {
  load(rdata_file)
}

# load libraries 
library(bacon)

# obtain data frame names from environment 
dataframe_names <- ls()[sapply(ls(), function(x) is.data.frame(get(x)))]
print(dataframe_names)

# apply function to obtain bacon estimates + plot p values distribution 
performBaconAnalysis(overlap_ewas_450k_model1.data)
performBaconAnalysis(overlap_ewas_450k_model2.data)
performBaconAnalysis(overlap_ewas_450k_trimester1_model1.data)
performBaconAnalysis(overlap_ewas_450k_trimester1_model2.data)
performBaconAnalysis(overlap_ewas_450k_trimester2_model1.data)
performBaconAnalysis(overlap_ewas_450k_trimester2_model2.data)
performBaconAnalysis(overlap_ewas_450k_trimester3_model1.data)
performBaconAnalysis(overlap_ewas_450k_trimester3_model2.data)

performBaconAnalysis(overlap_ewas_epic_model1.data)
performBaconAnalysis(overlap_ewas_epic_model2.data)
performBaconAnalysis(overlap_ewas_epic_trimester1_model1.data)
performBaconAnalysis(overlap_ewas_epic_trimester1_model2.data)
performBaconAnalysis(overlap_ewas_epic_trimester2_model1.data)
performBaconAnalysis(overlap_ewas_epic_trimester2_model2.data)
performBaconAnalysis(overlap_ewas_epic_trimester3_model1.data)
performBaconAnalysis(overlap_ewas_epic_trimester3_model2.data)

performBaconAnalysis(epic_only_ewas_model1.data)
performBaconAnalysis(epic_only_ewas_model2.data)
performBaconAnalysis(epic_only_ewas_epic_trimester1_model1.data)
performBaconAnalysis(epic_only_ewas_epic_trimester1_model2.data)
performBaconAnalysis(epic_only_ewas_epic_trimester2_model1.data)
performBaconAnalysis(epic_only_ewas_epic_trimester2_model2.data)
performBaconAnalysis(epic_only_ewas_epic_trimester3_model1.data)
performBaconAnalysis(epic_only_ewas_epic_trimester3_model2.data)

#\ 

###----------------STEP 4: Metafor meta-analyses----------------###

# load library for meta analysis and to save results
library(metafor)

# we will run an inverse weighted fixed-effects meta-analysis for total infection, trimester 1 infection, trimester 2 infection and trimester 3 infection. 
# For each exposure type we have a model 1 and 2. 
# We will meta-analyze the output of 16 ewases of the 450k and epic array together
# This will result in a final of 8 ewases
# because the vectors in each dataframes are too big, we need to do this with parallel processing 

setwd('set_path_to_results_folder')

# combine datasets for each model and exposure together that we want to meta-analyze together 
total_infection_model1 <- rbind(overlap_ewas_450k_model1.data, overlap_ewas_epic_model1.data)
total_infection_model2 <- rbind(overlap_ewas_450k_model2.data, overlap_ewas_epic_model2.data)

trimester1_infection_model1 <- rbind(overlap_ewas_450k_trimester1_model1.data, overlap_ewas_epic_trimester1_model1.data)
trimester1_infection_model2 <- rbind(overlap_ewas_450k_trimester1_model2.data, overlap_ewas_epic_trimester1_model2.data)

trimester2_infection_model1 <- rbind(overlap_ewas_450k_trimester2_model1.data, overlap_ewas_epic_trimester2_model1.data)
trimester2_infection_model2 <- rbind(overlap_ewas_450k_trimester2_model2.data, overlap_ewas_epic_trimester2_model2.data)

trimester3_infection_model1 <- rbind(overlap_ewas_450k_trimester3_model1.data, overlap_ewas_epic_trimester3_model1.data)
trimester3_infection_model2 <- rbind(overlap_ewas_450k_trimester3_model2.data, overlap_ewas_epic_trimester3_model2.data)

# Use pbmclapply to run the models in parallel
# of note, restart R after every model (and be sure to save the output as datafra,e) otherwise R crashes 
results_total_infection_model1 <- pbmclapply(total_infection_model1$cpg, function(cpg) {
  subset_data <- total_infection_model1[total_infection_model1$cpg == cpg,] 
  result <- meta_analyze(subset_data) 
  c(cpg, result$beta, result$se, result$pval) 
}, mc.cores = 30) 

results_total_infection_model2 <- pbmclapply(total_infection_model2$cpg, function(cpg) {
  subset_data <- total_infection_model2[total_infection_model2$cpg == cpg,] 
  result <- meta_analyze(subset_data) 
  c(cpg, result$beta, result$se, result$pval)
}, mc.cores = 30) 

results_trimester1_infection_model1 <- pbmclapply(trimester1_infection_model1$cpg, function(cpg) {
  subset_data <- trimester1_infection_model1[trimester1_infection_model1$cpg == cpg,] 
  result <- meta_analyze(subset_data) 
  c(cpg, result$beta, result$se, result$pval)
}, mc.cores = 30) 

results_trimester1_infection_model2 <- pbmclapply(trimester1_infection_model2$cpg, function(cpg) {
  subset_data <- trimester1_infection_model2[trimester1_infection_model2$cpg == cpg,] 
  result <- meta_analyze(subset_data) 
  c(cpg, result$beta, result$se, result$pval)
}, mc.cores = 30) 

results_trimester2_infection_model1 <- pbmclapply(trimester2_infection_model1$cpg, function(cpg) {
  subset_data <- trimester2_infection_model1[trimester2_infection_model1$cpg == cpg,] 
  result <- meta_analyze(subset_data) 
  c(cpg, result$beta, result$se, result$pval)
}, mc.cores = 30) 

results_trimester2_infection_model2 <- pbmclapply(trimester2_infection_model2$cpg, function(cpg) {
  subset_data <- trimester2_infection_model2[trimester2_infection_model2$cpg == cpg,] 
  result <- meta_analyze(subset_data) 
  c(cpg, result$beta, result$se, result$pval)
}, mc.cores = 40) 

results_trimester3_infection_model1 <- pbmclapply(trimester3_infection_model1$cpg, function(cpg) {
  subset_data <- trimester3_infection_model1[trimester3_infection_model1$cpg == cpg,] 
  result <- meta_analyze(subset_data) 
  c(cpg, result$beta, result$se, result$pval)
}, mc.cores = 30) 

results_trimester3_infection_model2 <- pbmclapply(trimester3_infection_model2$cpg, function(cpg) {
  subset_data <- trimester3_infection_model2[trimester3_infection_model2$cpg == cpg,] 
  result <- meta_analyze(subset_data) 
  c(cpg, result$beta, result$se, result$pval)
}, mc.cores = 30) 

# Obtain the summary of the meta-analysis results
result_total_infection_model1_df <- as.data.frame(matrix(unlist(results_total_infection_model1), ncol = 4, byrow = TRUE))
colnames(result_total_infection_model1_df) <- c("cpg","meta_beta", "meta_se", "meta_pval")
result_total_infection_model1_df_unique <- subset(result_total_infection_model1_df, !duplicated(cpg))

result_total_infection_model2_df <- as.data.frame(matrix(unlist(results_total_infection_model2), ncol = 4, byrow = TRUE))
colnames(result_total_infection_model2_df) <- c("cpg","meta_beta", "meta_se", "meta_pval")
result_total_infection_model2_df_unique <- subset(result_total_infection_model2_df, !duplicated(cpg))

results_trimester1_infection_model1_df <- as.data.frame(matrix(unlist(results_trimester1_infection_model1), ncol = 4, byrow = TRUE))
colnames(results_trimester1_infection_model1_df) <- c("cpg","meta_beta", "meta_se", "meta_pval")
results_trimester1_infection_model1_df_unique <- subset(results_trimester1_infection_model1_df, !duplicated(cpg))

results_trimester1_infection_model2_df <- as.data.frame(matrix(unlist(results_trimester1_infection_model2), ncol = 4, byrow = TRUE))
colnames(results_trimester1_infection_model2_df) <- c("cpg","meta_beta", "meta_se", "meta_pval")
results_trimester1_infection_model2_df_unique <- subset(results_trimester1_infection_model2_df, !duplicated(cpg))

results_trimester2_infection_model1_df <- as.data.frame(matrix(unlist(results_trimester2_infection_model1), ncol = 4, byrow = TRUE))
colnames(results_trimester2_infection_model1_df) <- c("cpg","meta_beta", "meta_se", "meta_pval")
results_trimester2_infection_model1_df_unique <- subset(results_trimester2_infection_model1_df, !duplicated(cpg))

results_trimester2_infection_model2_df <- as.data.frame(matrix(unlist(results_trimester2_infection_model2), ncol = 4, byrow = TRUE))
colnames(results_trimester2_infection_model2_df) <- c("cpg","meta_beta", "meta_se", "meta_pval")
results_trimester2_infection_model2_df_unique <- subset(results_trimester2_infection_model2_df, !duplicated(cpg))

results_trimester3_infection_model1_df <- as.data.frame(matrix(unlist(results_trimester3_infection_model1), ncol = 4, byrow = TRUE))
colnames(results_trimester3_infection_model1_df) <- c("cpg","meta_beta", "meta_se", "meta_pval")
results_trimester3_infection_model1_df_unique <- subset(results_trimester3_infection_model1_df, !duplicated(cpg))

results_trimester3_infection_model2_df <- as.data.frame(matrix(unlist(results_trimester3_infection_model2), ncol = 4, byrow = TRUE))
colnames(results_trimester3_infection_model2_df) <- c("cpg","meta_beta", "meta_se", "meta_pval")
results_trimester3_infection_model2_df_unique <- subset(results_trimester3_infection_model2_df, !duplicated(cpg))

# Save the results to an Excel file
write.csv(result_total_infection_model1_df_unique, "meta_analyzed_results_total_infection_model1.csv", row.names = FALSE)
write.csv(result_total_infection_model2_df_unique, "meta_analyzed_results_total_infection_model2.csv", row.names = FALSE)

write.csv(results_trimester1_infection_model1_df_unique, "meta_analyzed_results_trimester1_infection_model1.csv", row.names = FALSE)
write.csv(results_trimester1_infection_model2_df_unique, "meta_analyzed_results_trimester1_infection_model2.csv", row.names = FALSE)

write.csv(results_trimester2_infection_model1_df_unique, "meta_analyzed_results_trimester2_infection_model1.csv", row.names = FALSE)
write.csv(results_trimester2_infection_model2_df_unique, "meta_analyzed_results_trimester2_infection_model2.csv", row.names = FALSE)

write.csv(results_trimester3_infection_model1_df_unique, "meta_analyzed_results_trimester3_infection_model1.csv", row.names = FALSE)
write.csv(results_trimester3_infection_model2_df_unique, "meta_analyzed_results_trimester3_infection_model2.csv", row.names = FALSE)

#\ 

###----------------STEP 5: Quality checks meta-analysis EWAS----------------###

## run this part on the gena server 

## QC of EWAS 

# Set the path to the folder containing .csv files
setwd('set_path_to_meta_analysis_results_folder')

# List all CSV files in the current directory
csv_files <- list.files(pattern = "\\.csv$")

# Loop through each CSV file, read it, and assign it to an object with the original name
for (file in csv_files) {
  # Extract the file name without the file extension
  name <- tools::file_path_sans_ext(file)
  
  # Read the CSV file and assign it to an object with the original name
  assign(name, read.csv(file))
}

# Now, each data frame is available in your environment with its original name
# You can access them directly, e.g., `dataframe_name` contains the data from "dataframe_name.csv"

# load libraries 
library(bacon)

# obtain data frame names from environment 
dataframe_names <- ls()[sapply(ls(), function(x) is.data.frame(get(x)))]
print(dataframe_names)

# calculate lambda for meta analyzed results
Lambda<-function(P){
  chisq <- qchisq(1-P,1)
  median(chisq,na.rm=T)/qchisq(0.5,1)
}

Lambda(meta_analyzed_results_total_infection_model1$meta_pval)
Lambda(meta_analyzed_results_trimester1_infection_model1$meta_pval)
Lambda(meta_analyzed_results_trimester2_infection_model1$meta_pval)
Lambda(meta_analyzed_results_trimester3_infection_model1$meta_pval)

#\

###----------------STEP 6: Annotation, creating basic tables----------------###
## run this part of script on gena server 

# set wd 
setwd('set_path_to_meta_analysis_results_folder')

## now we will annotate the files 
library(IlluminaHumanMethylation450kanno.ilmn12.hg19) 

# first annotate dataframes with chromosome and position 
annotation <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
annotation$cpg <- annotation$Name # get cosnistent names 

# Specify the path to your folder containing the CSV files
csv_folder <- getwd()

# Get a list of CSV file names in the folder
csv_files <- list.files(csv_folder, pattern = "\\.csv$", full.names = TRUE)

# Loop through each CSV file and load it into a data frame
data_frames <- list()

for (csv_file in csv_files) {
  # Extract the original name without the .csv extension
  df_name <- tools::file_path_sans_ext(basename(csv_file))
  
  # Load the CSV file into a data frame
  temp_df <- read.csv(csv_file)
  
  # Perform the merge with annotation (adjust this part as needed)
  temp_df <- merge(temp_df, annotation, by = 'cpg', all.x = TRUE)
  
  # Save merged dataframe with 'merged_' prefix and export as CSV
  output_filename <- paste0('merged_', df_name, '.csv')
  write.csv(temp_df, output_filename, row.names = FALSE)
}

# first make tables with only p < 0.05 (nominal significance)

# Loop through each CSV file and process it
for (csv_file in csv_files) {
  # Extract the original name without the .csv extension
  df_name <- tools::file_path_sans_ext(basename(csv_file))
  
  # Load the CSV file into a data frame
  temp_df <- read.csv(csv_file)
  
  # Subset the data based on meta_pval < 0.05
  subset_df <- temp_df[temp_df$meta_pval < 0.05, ]
  
  # Save the subset dataframe with 'sign_results_' prefix and export as CSV
  output_filename <- paste0('sign_results_', df_name, '.csv')
  write.csv(subset_df, output_filename, row.names = FALSE)
}

#\

## now go back to local r studio 

## clean up tables for SI in manuscript 

## first for individual ewases

# Set the directory containing .rdata files
#setwd("set_path_to_results_EWAS_output")
dataframes_folder <- "path_to_individual_EWAS_output"

# List all .rdata files in the folder
rdata_files <- list.files(dataframes_folder, pattern = "\\.rdata$", full.names = TRUE)

# Function to subset significant results and save excel files
apply_filter_write_excel <- function(input_df, output_file) {
  
  # Subset the dataframe based on Pr(>|t|) < 0.05
  sign_df <- subset(input_df, `Pr(>|t|)` < 0.05)
  
  # Write the filtered dataframe to an Excel file
  write_xlsx(sign_df, output_file)
  
  cat("Filtered dataframe has been written to", output_file, "\n")
}

# set wd where we want to store results
setwd("set_wd_to_results_folder")

# load libary to save to excel 
library(writexl)

# Loop through each .rdata file, load, and apply the function
for (rdata_file in rdata_files) {
  
  # Load the dataframe from the .rdata file
  loaded_data <- get(load(rdata_file))
  
  # Extract the dataframe name from the file name
  df_name <- sub(".rdata$", "", basename(rdata_file))
  
  # Construct the output file name
  output_file <- paste0("sign_", df_name, ".xlsx")
  
  # Apply the function to the loaded dataframe
  apply_filter_write_excel(loaded_data, output_file)
}

#\

## then for meta-analyses 
rm(list = ls())

library(dplyr)

#setwd("path_to_annotated_results_folder")
csv_folder <- "path_to_meta_analyzed_results"

setwd("path_to_store_results")

# Get a list of CSV file names in the folder
csv_files <- list.files(csv_folder, pattern = "\\.csv$", full.names = TRUE)

# Loop through each CSV file
for (csv_file in csv_files) {
  
  # Read the CSV file into a dataframe
  original_df <- read.csv(csv_file)
  
  # Process the dataframe (modify as needed)
  processed_df <- dplyr::select(original_df, c('cpg', 'meta_beta', 'meta_se', 'meta_pval', 'chr', 'pos', 'Relation_to_Island'))
  
  cleaned_df <- subset(processed_df, meta_pval < 0.05)
  
  # Get the file name without extension
  file_name <- tools::file_path_sans_ext(basename(csv_file))
  
  # Construct the output Excel file name
  output_file <- paste0("cleaned_sign_", file_name, ".xlsx")
  
  # Write the cleaned dataframe to an Excel file
  write_xlsx(cleaned_df, output_file)
  
  cat("Cleaned dataframe has been written to", output_file, "\n")
}

#\ 

###----------------STEP 7: Post processing----------------###
## run this part on local r studio 

setwd("path_to_store_nominal_significant_annotated_results")

## first load in all csv files for nominal significant results 

# List all CSV files in the current directory
csv_files <- list.files(pattern = "\\.csv$")

# Loop through each CSV file, read it, and assign it to an object with the original name
for (file in csv_files) {
  # Extract the file name without the file extension
  name <- tools::file_path_sans_ext(file)
  
  # Read the CSV file and assign it to an object with the original name
  assign(name, read.csv(file))
}

## load in text files for cross reactive probes and probes containing snps
setwd("path_to_postprocessing_results")

## read crossreactive & flagged probes: combination of Naeem and Chen 
cross_reactive_probes<-read.csv("crossreactiveprobes.csv")
flagged_probes<-read.csv("flaggedprobes.csv")

## exclude cross reactive probes 
cleaned_sign_results_meta_analyzed_results_total_infection_model1 <- subset(sign_results_meta_analyzed_results_total_infection_model1, !(sign_results_meta_analyzed_results_total_infection_model1$cpg %in% cross_reactive_probes$MarkerName))

cleaned_sign_results_meta_analyzed_results_total_infection_model2 <- subset(sign_results_meta_analyzed_results_total_infection_model2, !(sign_results_meta_analyzed_results_total_infection_model2$cpg %in% cross_reactive_probes$MarkerName))

cleaned_sign_results_meta_analyzed_results_trimester1_infection_model1 <- subset(sign_results_meta_analyzed_results_trimester1_infection_model1, !(sign_results_meta_analyzed_results_trimester1_infection_model1$cpg %in% cross_reactive_probes$MarkerName))

cleaned_sign_results_meta_analyzed_results_trimester1_infection_model2 <- subset(sign_results_meta_analyzed_results_trimester1_infection_model2, !(sign_results_meta_analyzed_results_trimester1_infection_model2$cpg %in% cross_reactive_probes$MarkerName))

cleaned_sign_results_meta_analyzed_results_trimester2_infection_model1 <- subset(sign_results_meta_analyzed_results_trimester2_infection_model1, !(sign_results_meta_analyzed_results_trimester2_infection_model1$cpg %in% cross_reactive_probes$MarkerName))

cleaned_sign_results_meta_analyzed_results_trimester2_infection_model2 <- subset(sign_results_meta_analyzed_results_trimester2_infection_model2, !(sign_results_meta_analyzed_results_trimester2_infection_model2$cpg %in% cross_reactive_probes$MarkerName))

cleaned_sign_results_meta_analyzed_results_trimester3_infection_model1 <- subset(sign_results_meta_analyzed_results_trimester3_infection_model1, !(sign_results_meta_analyzed_results_trimester3_infection_model1$cpg %in% cross_reactive_probes$MarkerName))

cleaned_sign_results_meta_analyzed_results_trimester3_infection_model2 <- subset(sign_results_meta_analyzed_results_trimester3_infection_model2, !(sign_results_meta_analyzed_results_trimester3_infection_model2$cpg %in% cross_reactive_probes$MarkerName))

## flag probes containing snps
# first align name with cpg name of our meta analysis results file 
colnames(flagged_probes)[1] <- "cpg"

final_sign_results_meta_analyzed_results_total_infection_model1 <- merge(cleaned_sign_results_meta_analyzed_results_total_infection_model1, flagged_probes, by.x = "cpg", by.y = "cpg", all.x = T, all.y = F)
final_sign_results_meta_analyzed_results_total_infection_model2 <- merge(cleaned_sign_results_meta_analyzed_results_total_infection_model2, flagged_probes, by.x = "cpg", by.y = "cpg", all.x = T, all.y = F)

final_cleaned_sign_results_meta_analyzed_results_trimester1_infection_model1 <- merge(cleaned_sign_results_meta_analyzed_results_trimester1_infection_model1, flagged_probes, by.x = "cpg", by.y = "cpg", all.x = T, all.y = F)
final_cleaned_sign_results_meta_analyzed_results_trimester1_infection_model2 <- merge(cleaned_sign_results_meta_analyzed_results_trimester1_infection_model2, flagged_probes, by.x = "cpg", by.y = "cpg", all.x = T, all.y = F)

final_cleaned_sign_results_meta_analyzed_results_trimester2_infection_model1 <- merge(cleaned_sign_results_meta_analyzed_results_trimester2_infection_model1, flagged_probes, by.x = "cpg", by.y = "cpg", all.x = T, all.y = F)
final_cleaned_sign_results_meta_analyzed_results_trimester2_infection_model2 <- merge(cleaned_sign_results_meta_analyzed_results_trimester2_infection_model2, flagged_probes, by.x = "cpg", by.y = "cpg", all.x = T, all.y = F)

final_cleaned_sign_results_meta_analyzed_results_trimester3_infection_model1 <- merge(cleaned_sign_results_meta_analyzed_results_trimester3_infection_model1, flagged_probes, by.x = "cpg", by.y = "cpg", all.x = T, all.y = F)
final_cleaned_sign_results_meta_analyzed_results_trimester3_infection_model2 <- merge(cleaned_sign_results_meta_analyzed_results_trimester3_infection_model2, flagged_probes, by.x = "cpg", by.y = "cpg", all.x = T, all.y = F)

## save new datasets
setwd("path_to_cleaned_postprocessing_results")

library(writexl)

write_xlsx(final_sign_results_meta_analyzed_results_total_infection_model1, 'final_sign_results_meta_analyzed_results_total_infection_model1.xlsx')
write_xlsx(final_sign_results_meta_analyzed_results_total_infection_model2, 'final_sign_results_meta_analyzed_results_total_infection_model2.xlsx')

write_xlsx(final_cleaned_sign_results_meta_analyzed_results_trimester1_infection_model1, 'final_cleaned_sign_results_meta_analyzed_results_trimester1_infection_model1.xlsx')
write_xlsx(final_cleaned_sign_results_meta_analyzed_results_trimester1_infection_model2, 'final_cleaned_sign_results_meta_analyzed_results_trimester1_infection_model2.xlsx')

write_xlsx(final_cleaned_sign_results_meta_analyzed_results_trimester2_infection_model1, 'final_cleaned_sign_results_meta_analyzed_results_trimester2_infection_model1.xlsx')
write_xlsx(final_cleaned_sign_results_meta_analyzed_results_trimester2_infection_model2, 'final_cleaned_sign_results_meta_analyzed_results_trimester2_infection_model2.xlsx')

write_xlsx(final_cleaned_sign_results_meta_analyzed_results_trimester3_infection_model1, 'final_cleaned_sign_results_meta_analyzed_results_trimester3_infection_model1.xlsx')
write_xlsx(final_cleaned_sign_results_meta_analyzed_results_trimester3_infection_model2, 'final_cleaned_sign_results_meta_analyzed_results_trimester3_infection_model2.xlsx')

###----------------STEP 8: Manhattan + q-q plots----------------###
## run this part on local r studio 
 
# set wd 
setwd("path_to_annotated_meta_analysis")

# load data 
library(readr)

merged_meta_analyzed_results_total_infection_model2 <- read_csv("merged_meta_analyzed_results_total_infection_model2.csv")
merged_meta_analyzed_results_trimester1_infection_model2 <- read_csv("merged_meta_analyzed_results_trimester1_infection_model2.csv")
merged_meta_analyzed_results_trimester2_infection_model2 <- read_csv("merged_meta_analyzed_results_trimester2_infection_model2.csv")
merged_meta_analyzed_results_trimester3_infection_model2 <- read_csv("merged_meta_analyzed_results_trimester3_infection_model2.csv")

df1 <- merged_meta_analyzed_results_total_infection_model2
df2 <- merged_meta_analyzed_results_trimester1_infection_model2
df3 <- merged_meta_analyzed_results_trimester2_infection_model2
df4 <- merged_meta_analyzed_results_trimester3_infection_model2

## create manhattan plots
# https://danielroelfs.com/blog/how-i-create-manhattan-plots-using-ggplot/
library(ggplot2)
library(tidyverse)
library(ggtext)

# we will use negative log p value for manhattan plot
# advantages of this include: scale transformation, enhanced visual separation, more intuititive to itnerpret (e.g. height 2 represents a p val of 0.01, height 3 represents p val of 0.001 etc), easier highlights significance and reduces overplotting 

setwd("path_to_figures_folder")

generate_manhattan_plot <- function(df1) {
  # Calculate cumulative base pair position for each chromosome
  data_cum <- df1 %>%
    group_by(chr) %>%
    summarise(max_bp = max(pos)) %>%
    mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>%
    select(chr, bp_add)
  
  # Merge with the original dataframe and calculate cumulative base pair position
  new_df1 <- df1 %>%
    inner_join(data_cum, by = "chr") %>%
    mutate(bp_cum = pos + bp_add)
  
  # Calculate the center position for each chromosome
  axis_set <- new_df1 %>%
    group_by(chr) %>%
    summarize(center = mean(bp_cum))
  
  # Calculate y-axis limit for the plot
  ylim <- new_df1 %>%
    filter(meta_pval == min(meta_pval)) %>%
    mutate(ylim = abs(floor(log10(meta_pval))) + 2) %>%
    pull(ylim)
  
  # Set significance thresholds
  sig <- 0.05 / nrow(new_df1) #bonferroni correction
  sig2 <- 0.00005 #suggestive significance threshold 
  
  # Create Manhattan plot
  ggplot(new_df1, aes(
    x = bp_cum, y = -log10(meta_pval),
    color = as_factor(chr), size = -log10(meta_pval)
  )) +
    geom_hline(
      yintercept = -log10(sig), color = "red",
      linetype = "dashed"
    ) +
    geom_hline(
      yintercept = -log10(sig2), color = "purple",
      linetype = "dashed"
    ) +
    geom_point(alpha = 0.75) +
    scale_x_continuous(name = 'Chromosome',
                       breaks = axis_set$center,
                       labels = c("chr1", 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22')
    ) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, ylim)) +
    scale_color_manual(values = rep(
      c("#bc5090", "#ffa600"),
      unique(length(axis_set$chr))
    )) +
    scale_size_continuous(range = c(0.5, 3)) +
    labs(
      x = NULL,
      y = "-log<sub>10</sub>(p)"
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.title.y = element_markdown(),
      axis.text.x = element_text(angle = 45, size = 9, vjust = 0.5),
      panel.background = element_rect(fill = 'white', colour = 'black', size = 0.5)
    )
}

generate_manhattan_plot(df1)
generate_manhattan_plot(df2)
generate_manhattan_plot(df3)
generate_manhattan_plot(df4)

## create q-q plots 
#https://genome.sph.umich.edu/wiki/Code_Sample:_Generating_QQ_Plots_in_R
library(lattice)

qqunif.plot<-function(pvalues, 
                      should.thin=T, thin.obs.places=2, thin.exp.places=2, 
                      xlab=expression(paste("Expected (",-log[10], " p-value)")),
                      ylab=expression(paste("Observed (",-log[10], " p-value)")), 
                      draw.conf=TRUE, conf.points=1000, conf.col="#ffa600", conf.alpha=.05,
                      already.transformed=FALSE, pch=20, aspect="iso", prepanel=prepanel.qqunif,
                      par.settings=list(superpose.symbol=list(pch=pch)), ...) {
  
  
  #error checking
  if (length(pvalues)==0) stop("pvalue vector is empty, can't draw plot")
  if(!(class(pvalues)=="numeric" || 
       (class(pvalues)=="list" && all(sapply(pvalues, class)=="numeric"))))
    stop("pvalue vector is not numeric, can't draw plot")
  if (any(is.na(unlist(pvalues)))) stop("pvalue vector contains NA values, can't draw plot")
  if (already.transformed==FALSE) {
    if (any(unlist(pvalues)==0)) stop("pvalue vector contains zeros, can't draw plot")
  } else {
    if (any(unlist(pvalues)<0)) stop("-log10 pvalue vector contains negative values, can't draw plot")
  }
  
  
  grp<-NULL
  n<-1
  exp.x<-c()
  if(is.list(pvalues)) {
    nn<-sapply(pvalues, length)
    rs<-cumsum(nn)
    re<-rs-nn+1
    n<-min(nn)
    if (!is.null(names(pvalues))) {
      grp=factor(rep(names(pvalues), nn), levels=names(pvalues))
      names(pvalues)<-NULL
    } else {
      grp=factor(rep(1:length(pvalues), nn))
    }
    pvo<-pvalues
    pvalues<-numeric(sum(nn))
    exp.x<-numeric(sum(nn))
    for(i in 1:length(pvo)) {
      if (!already.transformed) {
        pvalues[rs[i]:re[i]] <- -log10(pvo[[i]])
        exp.x[rs[i]:re[i]] <- -log10((rank(pvo[[i]], ties.method="first")-.5)/nn[i])
      } else {
        pvalues[rs[i]:re[i]] <- pvo[[i]]
        exp.x[rs[i]:re[i]] <- -log10((nn[i]+1-rank(pvo[[i]], ties.method="first")-.5)/(nn[i]+1))
      }
    }
  } else {
    n <- length(pvalues)+1
    if (!already.transformed) {
      exp.x <- -log10((rank(pvalues, ties.method="first")-.5)/n)
      pvalues <- -log10(pvalues)
    } else {
      exp.x <- -log10((n-rank(pvalues, ties.method="first")-.5)/n)
    }
  }
  
  
  #this is a helper function to draw the confidence interval
  panel.qqconf<-function(n, conf.points=1000, conf.col="gray", conf.alpha=.05, ...) {
    require(grid)
    conf.points = min(conf.points, n-1);
    mpts<-matrix(nrow=conf.points*2, ncol=2)
    for(i in seq(from=1, to=conf.points)) {
      mpts[i,1]<- -log10((i-.5)/n)
      mpts[i,2]<- -log10(qbeta(1-conf.alpha/2, i, n-i))
      mpts[conf.points*2+1-i,1]<- -log10((i-.5)/n)
      mpts[conf.points*2+1-i,2]<- -log10(qbeta(conf.alpha/2, i, n-i))
    }
    grid.polygon(x=mpts[,1],y=mpts[,2], gp=gpar(fill=conf.col, lty=0), default.units="native")
  }
  
  #reduce number of points to plot
  if (should.thin==T) {
    if (!is.null(grp)) {
      thin <- unique(data.frame(pvalues = round(pvalues, thin.obs.places),
                                exp.x = round(exp.x, thin.exp.places),
                                grp=grp))
      grp = thin$grp
    } else {
      thin <- unique(data.frame(pvalues = round(pvalues, thin.obs.places),
                                exp.x = round(exp.x, thin.exp.places)))
    }
    pvalues <- thin$pvalues
    exp.x <- thin$exp.x
  }
  gc()
  
  prepanel.qqunif= function(x,y,...) {
    A = list()
    A$xlim = range(x, y)*1.02
    A$xlim[1]=0
    A$ylim = A$xlim
    return(A)
  }
  
  #draw the plot
  xyplot(pvalues~exp.x, groups=grp, xlab=xlab, ylab=ylab, aspect=aspect,
         prepanel=prepanel, scales=list(axs="i"), pch=pch,
         panel = function(x, y, ...) {
           if (draw.conf) {
             panel.qqconf(n, conf.points=conf.points, 
                          conf.col=conf.col, conf.alpha=conf.alpha)
           };
           panel.xyplot(x,y, ..., col = '#58508d');
           panel.abline(0,1, col = 'black');
         }, par.settings=par.settings, ...
  )
}


qqunif.plot(df1$meta_pval)
qqunif.plot(df2$meta_pval)
qqunif.plot(df3$meta_pval)
qqunif.plot(df4$meta_pval)

#\ 

###----------------STEP 9: Post-hoc analysis----------------###
#In this posthoc analysis we will create a sumscore in genR, similar to ALSPAC, only for influenza and urinary tract infections, and associate this with the suggestive hits from the EWAs

#First, load GenR EWAS sumstats (nominal significant results; cross-reactive cpgs are removed)
setwd('path_to_input_files_MPS')

library(readxl)
res_totalinf <- read_excel('final_sign_results_meta_analyzed_results_total_infection_model1.xlsx')

# remove unnecessary cols
res_totalinf <- dplyr::select(res_totalinf, -c('X', 'flag'))

# create vectors for suggestive sig cpgs
sigCpGs_totalinf <- res_totalinf$cpg[res_totalinf$meta_pval<5e-5]
length(sigCpGs_totalinf) #33

#load infection pheno data and create new infection score 
setwd('path_to_store_results_posthoc_analysis')

library(foreign)
pheno <- read.spss('Second_hits_DF.sav', to.data.frame = T)
pheno$alspac_infections <- pheno$flu_tri1 + pheno$flu_tri2 + pheno$flu_tri3 + pheno$UWI_tri1 + pheno$UWI_tri2 + pheno$UWI_tri3

setwd('path_to_data_folder')
df_final_450k <- readRDS('df_final_450k.rds')
df_final_epic <- readRDS('df_final_epic.rds')

pheno_all_450k <- merge(df_final_450k, pheno, by = 'IDC', all.x = T)
pheno_all_epic <- merge(df_final_epic, pheno, by = 'IDC', all.x = T)

#load genR methylation data
overlap_df_450k <- readRDS('overlap_df_450k.rds') 
overlap_df_epic <- readRDS('overlap_df_epic.rds') 

library(tibble)
overlap_df_450k <- tibble::rownames_to_column(overlap_df_450k, "Sample_ID")

# Merge methylation dataframe with pheno 
df_450k <- merge(pheno_all_450k, overlap_df_450k, by = 'Sample_ID', all.x = T)
df_epic <- cbind(pheno_all_epic, overlap_df_epic)

# Run regression between alspac like infection score in genR and top hits in genR (450k)
results_df_totalinf <- data.frame() 

for (x in sigCpGs_totalinf){
  # specify lm
  m1 <- paste0(x,"~ scale(alspac_infections) + GENDER.x + SMOKE_ALL.x + AGE_M_v2.x + PARITY + EDUCM_3groups + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC")
  m2 <- paste0(x,"~ scale(alspac_infections) + GENDER.x + SMOKE_ALL.x + AGE_M_v2.x + PARITY + EDUCM_3groups + GESTBIR.x + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC")
  # extract coefficients
  coef_totalinf_m1 <- summary(lm(as.formula(m1), data = df_450k))$coefficients[2,]
  coef_totalinf_m2 <- summary(lm(as.formula(m2), data = df_450k))$coefficients[2,]
  # save output 
  results_df_totalinf[x,c(1:4)] <- coef_totalinf_m1
  results_df_totalinf[x,c(5:8)] <- coef_totalinf_m2
}

colnames(results_df_totalinf) <- c("bval_m1", "seval_m1", "tval_m1","pval_m1","bval_m2", "seval_m2", "tval_m2","pval_m2")

results_df_totalinf$outcome <- sigCpGs_totalinf

results_df_totalinf$sign_m1 = "" 
results_df_totalinf$sign_m1[results_df_totalinf$pval_m1 < 0.05] = "*"
results_df_totalinf$sign_m1[results_df_totalinf$pval_m1 < 0.01] = "**"
results_df_totalinf$sign_m1[results_df_totalinf$pval_m1 < 0.001] = "***"
results_df_totalinf$sign_m2 = "" 
results_df_totalinf$sign_m2[results_df_totalinf$pval_m2 < 0.05] = "*"
results_df_totalinf$sign_m2[results_df_totalinf$pval_m2 < 0.01] = "**"
results_df_totalinf$sign_m2[results_df_totalinf$pval_m2 < 0.001] = "***"

library(writexl)
write_xlsx(results_df_totalinf, path = "alspac_like_infections_genr_top_hits_genr_450k.xlsx")

# Run regression between alspac like infection score in genR and top hits in genR (epic)
results_df_totalinf2 <- data.frame() 

for (x in sigCpGs_totalinf){
  # specify lm
  m1 <- paste0(x,"~ scale(alspac_infections) + GENDER.x + SMOKE_ALL.x + AGE_M_v2.x + PARITY + EDUCM_3groups + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC")
  m2 <- paste0(x,"~ scale(alspac_infections) + GENDER.x + SMOKE_ALL.x + AGE_M_v2.x + PARITY + EDUCM_3groups + GESTBIR.x + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC")
  # extract coefficients
  coef_totalinf_m1 <- summary(lm(as.formula(m1), data = df_epic))$coefficients[2,]
  coef_totalinf_m2 <- summary(lm(as.formula(m2), data = df_epic))$coefficients[2,]
  # save output 
  results_df_totalinf2[x,c(1:4)] <- coef_totalinf_m1
  results_df_totalinf2[x,c(5:8)] <- coef_totalinf_m2
}

colnames(results_df_totalinf2) <- c("bval_m1", "seval_m1", "tval_m1","pval_m1","bval_m2", "seval_m2", "tval_m2","pval_m2")

results_df_totalinf2$outcome <- sigCpGs_totalinf

results_df_totalinf2$sign_m1 = "" 
results_df_totalinf2$sign_m1[results_df_totalinf2$pval_m1 < 0.05] = "*"
results_df_totalinf2$sign_m1[results_df_totalinf2$pval_m1 < 0.01] = "**"
results_df_totalinf2$sign_m1[results_df_totalinf2$pval_m1 < 0.001] = "***"
results_df_totalinf2$sign_m2 = "" 
results_df_totalinf2$sign_m2[results_df_totalinf2$pval_m2 < 0.05] = "*"
results_df_totalinf2$sign_m2[results_df_totalinf2$pval_m2 < 0.01] = "**"
results_df_totalinf2$sign_m2[results_df_totalinf2$pval_m2 < 0.001] = "***"

write_xlsx(results_df_totalinf2, path = "alspac_like_infections_genr_top_hits_genr_epic.xlsx")

# Meta-analyze output 450k & epic
df1 <- read_xlsx('alspac_like_infections_genr_top_hits_genr_epic.xlsx')
df2<- read_xlsx('alspac_like_infections_genr_top_hits_genr_450k.xlsx')

combined_results <- rbind(df1, df2)
combined_results <- as.data.frame(combined_results)

library(metafor)
meta_analyze <- function(data) {
  rma(yi = data$bval_m1, sei = data$seval_m1, method = "FE")
}

library(pbmcapply)
results <- pbmclapply(combined_results$outcome, function(outcome) {
  subset_data <- combined_results[combined_results$outcome == outcome,] 
  result <- meta_analyze(subset_data) 
  c(outcome, result$beta, result$se, result$pval) 
}, mc.cores = 1) 


results <- as.data.frame(matrix(unlist(results), ncol = 4, byrow = TRUE))
colnames(results) <- c("outcome","meta_beta", "meta_se", "meta_pval")
results <- subset(results, !duplicated(outcome))

results$sign = "" 
results$sign[results$meta_pval < 0.05] = "*"
results$sign[results$meta_pval < 0.01] = "**"
results$sign[results$meta_pval < 0.001] = "***"

write_xlsx(results, 'MA_alspac_like_infections_genr_top_hits_genr.xlsx')

#\ END OF SCRIPT 
