#########################################################################
#####################PROJECT: PRENATAL INFECTION EWAS####################
#########################################################################

#Author: Anna Suleri

### Parts in this script
# Part 0: Loading relevant libraries, source functions file and setting working directory
# Part 1: Loading and infection and covariate data 
# Part 2: Inclusion and exclusion criteria
# Part 3: Multiple imputation of covariates 
# Part 4: Baseline table 

#O f Note, in this script we prepare the data for the EWAS. The next part of the script is the EWAS script and that part is run on the server. 

###-------------------------------------------------###
########################PART 0#########################
###-------------------------------------------------###
## clears the environment
rm(list = ls()) 

# specify library path 

## Set seed 
set.seed(2023)

## Load relevant libraries 
libraries <- c('foreign', 'haven', 'tidyverse', 'dplyr', 'ggplot2', 'mice', 'miceadds', 'write_xl')

invisible(lapply(libraries, require, character.only = T))

## Set wd 
setwd('set_path_to_scripts')
wd <- getwd()

## Source functions file
source('0.Functions_script.R')

#\

###-------------------------------------------------###
########################PART 1#########################
###-------------------------------------------------###
setwd('set_path_to_data')

## Loading all files of interest 
# exposure and covariates 
df1 <- read.spss('Second_hits_DF.sav', to.data.frame = T)
ddf1 <- dplyr::select(df1, c('IDC', 'IDM', 'MOTHER', 'GENDER', 'SMOKE_ALL', "WEIGHT", "GESTBIR", "AGE_M_v2",'sumscore_inf_tot', 'sumscore_inf_tri1', 'sumscore_inf_tri2', 'sumscore_inf_tri3')) #selecting vars we need 
df2 <- read.spss('PARTUS_17072015.sav', to.data.frame = T)
ddf2 <- dplyr::select(df2, c('IDM', 'PARITY')) #selecting vars we need 
df3 <- read.csv('Covariates_behavior_long-behav_paper.csv')
ddf3 <- dplyr::select(df3, c('IDC', 'EDUCM_3groups'))
ddf3$EDUCM_3groups <- as.factor(ddf3$EDUCM_3groups) #from character to factor 
df4 <- read.spss("ETHNICITYCHILDsurinameandafricansubdivided_19072017.sav", to.data.frame = T)
df4$ethnicity <- as.factor(ifelse(df4$ETHNINFv5 == 'Dutch' | df4$ETHNINFv5 == 'American,western' | df4$ETHNINFv5 == 'European', 'caucasian', 'not caucasian')) #create binary caucasian variable

# 450k specific 
df5 <- read.spss('Selection_GENR_450kmeth_release3_birth_20180125.sav', to.data.frame = T)
ddf5 <- dplyr::select(df5, -c('Array_number', 'Position')) #remove columns we don't need
ddf5$Sample_ID <- gsub(" ", "", ddf5$Sample_ID, fixed = T) #remove spaces 
df6 <- read.csv('GENR_450kmeth_release3_birth_Salas.csv') #450k celltype 
df6 <- rename(df6, 'Sample_ID' = Rownames) #to align with later files 
df6$Sample_ID <- gsub(" ", "", df6$Sample_ID, fixed = T) #remove spaces 

# EPIC specific 
df7 <- read.spss('Selection_GENR_MethylEPIC_release1_birth_20230627.sav', to.data.frame = T) #recode to different numbers so that batch plates for 450k do not have same number plate as batch plates of epic assay
ddf7 <- subset(df7, EwasChildMethylEpic == 1) #only include children with suitable cpgs for ewas 
load('GENR_EPICv1METH_birth_CellTypes_combined.Rdata') #load epic cell type proportions
df8 <- GENR_EPICv1METH_birth_CellTypes_combined.data

## Merging all files into a final dataframe for 450k and epic seperately 
merge1 <- merge(ddf1, ddf2, by = 'IDM', all.x = T)
merge2 <- merge(merge1, ddf3, by = 'IDC', all.x = T)
merge3 <- merge(merge2, df4, by = 'IDC', all.x = T)

merge_450k_1 <- merge(merge3, ddf5, by = 'IDC', all.x = T)
merge_450k_2 <- merge(merge_450k_1, df6, by = 'Sample_ID', all.x = T) #full 450k file

merge_epic_1 <- merge(merge3, ddf7, by = 'IDC', all.x = T)
merge_epic_2 <- merge(merge_epic_1, df8, by = 'SampleID', all.x = T) #full epic file 
#\

###-------------------------------------------------###
########################PART 2#########################
###-------------------------------------------------###

# inclusion criteria: good qc epigenetic data 
# exclusion criteria: both twins and one sibling for each pair 
# we start with n = 9901 children 

### 450k sample----
# Inclusion criteria 
include_450k_df1 <- subset(merge_450k_2, complete.cases(Sample_ID)) #n = 1396
include_450k_df2 <- subset(include_450k_df1, ethnicity == 'caucasian') #n =1382

# Exclusion criteria
include_450k_df3 <- include_450k_df2[sample(nrow(include_450k_df2)),]
include_450k_df3$na_count <- apply(include_450k_df3, 1, function(x) sum(is.na(x))) 
include_450k_df4 <- include_450k_df3[order(include_450k_df3$na_count),] 
include_450k_df5 <- include_450k_df4[!duplicated(include_450k_df4$MOTHER, fromLast = T),] # n = 1367
summary(duplicated(include_450k_df5$IDM)) # double check that all twins are removed

### EPIC sample----
include_epic_df1 <- subset(merge_epic_2, complete.cases(SampleID)) #n = 1115 
include_epic_df2 <- subset(include_epic_df1, ethnicity == 'caucasian') #n = 985

# Exclusion criteria
include_epic_df3 <- include_epic_df2[sample(nrow(include_epic_df2)),]
include_epic_df3$na_count <- apply(include_epic_df3, 1, function(x) sum(is.na(x))) 
include_epic_df4 <- include_epic_df3[order(include_epic_df3$na_count),] 
include_epic_df5 <- include_epic_df4[!duplicated(include_epic_df4$MOTHER, fromLast = T),] # n = 971
summary(duplicated(include_epic_df5$IDM)) # double check that all twins are removed

#\ 

### Creating baseline table 
library(xlsx)
setwd('set_path_to_results')
baselinevars <- c('GENDER', 'GESTBIR', 'WEIGHT', 'AGE_M_v2','PARITY', 'EDUCM_3groups' ,'SMOKE_ALL', 'sumscore_inf_tot', 'sumscore_inf_tri1', 'sumscore_inf_tri2', 'sumscore_inf_tri3')

baseline_table_function(baselinevars, include_450k_df5)
baseline_table_function(baselinevars, include_epic_df5)

### calculate infection type frequencies for plot later 
setwd("set_path_to_daata")

# load dataframe that has all individual infection types 
infections_df <- read.spss('Second_hits_DF.sav', to.data.frame = T)

# only select participants in that df that are in epigenetics paper
infections_df2 <- as.data.frame(infections_df[infections_df$IDC %in% include_450k_df5$IDC,])

infections_df3 <- as.data.frame(infections_df[infections_df$IDC %in% include_epic_df5$IDC,])

# calculate percentage of 'yes' for each infection type per trimester
infection_columns <- c('clean_upper_resp_inf_tri1', 'clean_upper_resp_inf_tri2', 'clean_upper_resp_inf_tri3', 'clean_GI_inf_tri1', 'clean_GI_inf_tri2', 'clean_GI_inf_tri3', 'flu_tri1', 'flu_tri2', 'flu_tri3', 'clean_uwi_tri1', 'clean_uwi_tri2', 'clean_uwi_tri3', 'dermatitis_tri1', 'dermatitis_tri2', 'dermatitis_tri3', 'clean_lower_resp_inf_tri1', 'clean_lower_resp_inf_tri2', 'clean_lower_resp_inf_tri3', 'herpeszoster_tri1', 'herpeszoster_tri2', 'herpeszoster_tri3', 'jaundice_tri1', 'jaundice_tri2', 'jaundice_tri3', 'STD_tri1', 'STD_tri2', 'STD_tri3', 'fever_tri1', 'fever_tri2', 'fever_tri3')

infections_df2[, infection_columns] <- lapply(infections_df2[, infection_columns], factor) #from numeric to factor, so we have a yes/no variable 
infections_df3[, infection_columns] <- lapply(infections_df3[, infection_columns], factor) 

for(i in infection_columns) {
  yes_count <- sum(infections_df2[[i]] == 1, na.rm = T)
  total_responses <- sum(!is.na(infections_df2[[i]]))
  percentage_yes <- (yes_count / total_responses) * 100
  message(i)
  print(percentage_yes)
}

for(i in infection_columns) {
  yes_count <- sum(infections_df3[[i]] == 1, na.rm = T)
  total_responses <- sum(!is.na(infections_df3[[i]]))
  percentage_yes <- (yes_count / total_responses) * 100
  message(i)
  print(percentage_yes)
}

#\ 

###-------------------------------------------------###
########################PART 3#########################
###-------------------------------------------------###

# We will impute infection score for each trimester and covariates; then we will select last imputed dataset and use that for the downstream analyses 

### 450k sample ----
# first check missingness 
miss_values_function(include_450k_df5)

# check str of dataframe if all vars are coded correct and recode if necessary and remove unnecessary vars
include_450k_df6 <- dplyr::select(include_450k_df5, -c('ETHNINFv5', 'na_count'))
str(include_450k_df6)

# apply mice function 
imputation_function(data = include_450k_df6, exclude_imp_vars = c(1:4, 16:24), exclude_predictors = c(1:4, 16:24), method = "default")

imp_450k_df <- readRDS('imputedData_include_450k_df6.rds') #load imputed df 

# check quality imputation 
plot(imp_450k_df)
imp_450k_df$loggedEvents

# select last imputed df
df_final_450k <- complete(imp_450k_df, 30)

# scale infection scores
sd_vars <- c('sumscore_inf_tot', 'sumscore_inf_tri1', 'sumscore_inf_tri2', 'sumscore_inf_tri3')

for (x in sd_vars){ 
  t <- df_final_450k[x]
  colname <- paste0(colnames(t), "_standardized")
  df_final_450k[,colname] <- as.numeric(scale(t))
}

#c heck if there are indeed no missingness after the imputation (if there are, mice did not succeed/converge)
miss_values_function(df_final_450k)

# save final df to continue EWAS with 
saveRDS(df_final_450k, 'df_final_450k.rds')

#\ 

### EPIC sample ----
# first check missingness 
miss_values_function(include_epic_df5)

# check str of dataframe if all vars are coded correct and recode if necessary and remove unnecessary vars
include_epic_df6 <- dplyr::select(include_epic_df5, -c('ETHNINFv5', 'na_count', 'EwasChildMethylEpic', 'na_count'))
str(include_epic_df6)
include_epic_df6$Sample_Plate <- as.numeric(include_epic_df6$Sample_Plate) 

# apply mice function 
imputation_function(data = include_epic_df6, exclude_imp_vars = c(1:4, 16:24), exclude_predictors = c(1:4, 16:24), method = "default")

imp_epic_df <- readRDS('imputedData_include_epic_df6.rds') #load imputed df

# check quality imputation 
plot(imp_epic_df)
imp_epic_df$loggedEvents

# select last imputed df
df_final_epic <- complete(imp_epic_df, 30)

# scale infection scores
for (x in sd_vars){ 
  t <- df_final_epic[x]
  colname <- paste0(colnames(t), "_standardized")
  df_final_epic[,colname] <- as.numeric(scale(t))
}

# check if there are indeed no missingness after the imputation (if there are, mice did not succeed/converge)
miss_values_function(df_final_epic)

# save final df to continue EWAS with 
saveRDS(df_final_epic, 'df_final_epic.rds')

#\ END OF THIS SCRIPT.
