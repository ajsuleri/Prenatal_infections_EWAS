#########################################################################
#####################PROJECT: PRENATAL INFECTION EWAS####################
#########################################################################

# Author: Anna Suleri & Mannan Luo 

#In this part of the script we will create a methylation profile score for infections based on the EWAS output. We will then validate the score with self-reported infections and some other child phenotypes within this discovery sample.

# script_p1 (all steps in the discovery sample): pre-selecting CpGs form EWAS sumstats; re-weighting effect size using ENR in training set; calculating MPS and checking prediction performance in testing set.

#Some point to note. 
#1)	Preselect CpGs from EWAS sumstats based on a p-value threshold (if you choose 0.05, it might be a heavy burden for computation, but you can try).Then reweight effect size using ENR in training set. Then calculate MPS. Then check prediction performance in test set. 
#2)	Create a covariate for array types, so that you can create arrays-equal-distribution data for training and test
#3)	Have Sample_ID (i.e., child's ID for DNAm) in your phenotype data, so that you can merge it with DNAm data. 

# Check reference links
#' glmnetUtils package https://cran.r-project.org/web/packages/glmnetUtils/vignettes/intro.html
#' glmnet package https://glmnet.stanford.edu/articles/glmnet.html
#' caret package  http://topepo.github.io/caret/index.html
#' https://www.science.smith.edu/~jcrouser/SDS293/labs/lab10-r.html
#' https://stats.stackexchange.com/questions/299653/caret-glmnet-vs-cv-glmnet
#' Interpreting results https://medium.com/analytics-vidhya/multicollinearity-ridge-lasso-elastic-net-regression-using-r-6582cbabf7f3
#' Helpful tutorial: https://www.analyticsvidhya.com/blog/2017/06/a-comprehensive-guide-for-linear-ridge-and-lasso-regression/

# Run this script on the GENA server!

###----------------STEP 1: Preselecting CpGs from EWAS----------------###

## We do this based on a specific p-value threshold; we start with nominal sig; so p < 0.05 

## Clean environment 
rm(list=ls())

# Load packages
library(data.table)
library(dplyr)
library(caret) 
library(glmnetUtils) # for cva.glmnet() function
library(readxl)

## Load EWAS sumstats (nominal significant results; cross-reactive cpgs are removed)
# set wd
setwd('set_path_to_input_files')

# load data
res_totalinf <- read_excel('final_sign_results_meta_analyzed_results_total_infection_model1.xlsx')
res_trim1 <- read_excel('final_cleaned_sign_results_meta_analyzed_results_trimester1_infection_model1.xlsx')
res_trim2 <- read_excel('final_cleaned_sign_results_meta_analyzed_results_trimester2_infection_model1.xlsx')
res_trim3 <- read_excel('final_cleaned_sign_results_meta_analyzed_results_trimester3_infection_model1.xlsx')

# remove unnecessary cols
res_totalinf <- dplyr::select(res_totalinf, -c('X', 'flag'))
res_trim1 <- dplyr::select(res_trim1, -c('flag'))
res_trim2 <- dplyr::select(res_trim2, -c('flag'))
res_trim3 <- dplyr::select(res_trim3, -c('flag'))

# create vectors for sig cpgs
sigCpGs_totalinf <- res_totalinf$cpg[res_totalinf$meta_pval<0.05]
length(sigCpGs_totalinf) #20.000
sigCpGs_trimester1 <- res_trim1$cpg[res_trim1$meta_pval<0.05]
length(sigCpGs_trimester1) #18.542
sigCpGs_trimester2 <- res_trim2$cpg[res_trim2$meta_pval<0.05]
length(sigCpGs_trimester2) #18.506
sigCpGs_trimester3 <- res_trim3$cpg[res_trim3$meta_pval<0.05]
length(sigCpGs_trimester3)#19.666

#\ 

###----------------STEP 2: Data preparation----------------###

## Step a) load phenotype data (i.e., your exposure/outcome + Sample_ID variable from methylation data + all covariates in your EWAS model )
pheno.dat.450k <- readRDS('df_final_450k.rds')
pheno.dat.epic <- readRDS('df_final_epic.rds')

# create assay vector for each df
pheno.dat.450k[, "assay"] <- rep('450k', length.out = nrow(pheno.dat.450k)) %>% factor()
pheno.dat.epic[, "assay"] <- rep('epic', length.out = nrow(pheno.dat.epic)) %>% factor()

# align sample ID names in both dataframes
pheno.dat.epic <- rename(pheno.dat.epic, 'Sample_ID' = SampleID)

# merge phenotype data
pheno.dat <- rbind(pheno.dat.450k, pheno.dat.epic)

# imputed data but double check that there is no missing because we want to use complete data 
pheno <-  pheno.dat[complete.cases(pheno.dat), ]

## Step b) load methylation data

# load methylation files for 450k and epic
load("450k_methylation_data") #beta matrix is called x
load("epic_methylation_data") 
epic_df <- GENR_EPICv1METH_Norm_Betas_birth_3IQRwinsorized.data #rename beta matrix of epic 

# Check dimensions to see if rownames and colnames are id and cpgs (so identitical for both dataframes, if not, make it identitical)
dim(x) #rows = 458563, columns = 1396
dim(epic_df) #rows = 808183, columns = 1115

# Create df of cpgs that are available in both 450k assay and epic assay 
overlapping_betas <- as.data.frame(x[rownames(x) %in% rownames(epic_df),]) #393360 cpgs 
overlapping_betas2 <- as.data.frame(epic_df[rownames(epic_df) %in% rownames(x),])

# combine methylation data 
methylation.dat <- cbind(overlapping_betas, overlapping_betas2)
dim(methylation.dat) #393360 rows and 2511 columns 

# check data in environment
ls()

## Step c) subset DNAm data based the preselected CpGs from EWAS
# first select cpgs from the methylation data that reached nominal significance 
meth.sub.totalinf<- methylation.dat[rownames(methylation.dat) %in% sigCpGs_totalinf,]

meth.sub.trimester1<- methylation.dat[rownames(methylation.dat) %in% sigCpGs_trimester1,]

meth.sub.trimester2<- methylation.dat[rownames(methylation.dat) %in% sigCpGs_trimester2,]

meth.sub.trimester3<- methylation.dat[rownames(methylation.dat) %in% sigCpGs_trimester3,]

# Check dimensions; rows are cpg sites and cols are sample id
dim(meth.sub.totalinf) 

dim(meth.sub.trimester1) 

dim(meth.sub.trimester2) 

dim(meth.sub.trimester3) 

# Get an example what the data looks like and if it is correct 
meth.sub.totalinf[1:5,1:5]

meth.sub.trimester1[1:5,1:5]

meth.sub.trimester2[1:5,1:5]

meth.sub.trimester3[1:5,1:5]

## Step d) match DNAm data and pheno data based on Sample_ID
# Make sure colnames in your meth data refers to Sample_ID, not CpG name
# Include identical individuals in meth and pheno files
meth <- meth.sub.totalinf[,na.omit(match(pheno$Sample_ID,colnames(meth.sub.totalinf)))]   
pheno <- pheno[match(colnames(meth),pheno$Sample_ID),]

meth2 <- meth.sub.trimester1[,na.omit(match(pheno$Sample_ID,colnames(meth.sub.trimester1)))]   
pheno2 <- pheno[match(colnames(meth2),pheno$Sample_ID),]

meth3<- meth.sub.trimester2[,na.omit(match(pheno$Sample_ID,colnames(meth.sub.trimester2)))]   
pheno3 <- pheno[match(colnames(meth3),pheno$Sample_ID),]

meth4 <- meth.sub.trimester3[,na.omit(match(pheno$Sample_ID,colnames(meth.sub.trimester3)))]   
pheno4 <- pheno[match(colnames(meth4),pheno$Sample_ID),]

# Double check if both data match well
identical(colnames(meth), pheno$Sample_ID) #needs to be TRUE
all(pheno$Sample_ID==colnames(meth))

identical(colnames(meth2), pheno2$Sample_ID)
all(pheno2$Sample_ID==colnames(meth2))

identical(colnames(meth3), pheno3$Sample_ID)
all(pheno3$Sample_ID==colnames(meth3))

identical(colnames(meth4), pheno4$Sample_ID)
all(pheno4$Sample_ID==colnames(meth4))

## Step e) Transpose methylation data so that rows are sample ID(person) and columns are probes(cpg sites)
# Transpose
meth.t <- t(meth)

meth.t2 <- t(meth2)

meth.t3 <- t(meth3)

meth.t4 <- t(meth4)

# check if done correctly 
ncol(meth.t)
sum(pheno$Sample_ID == rownames(meth.t))/nrow(meth.t) #should be 1
all.equal(as.character(pheno$Sample_ID),rownames(meth.t))  # must be TRUE

ncol(meth.t2)
sum(pheno$Sample_ID == rownames(meth.t2))/nrow(meth.t2) #should be 1
all.equal(as.character(pheno2$Sample_ID),rownames(meth.t2))  # must be TRUE

ncol(meth.t3)
sum(pheno$Sample_ID == rownames(meth.t3))/nrow(meth.t3) #should be 1
all.equal(as.character(pheno3$Sample_ID),rownames(meth.t3))  # must be TRUE

ncol(meth.t4)
sum(pheno$Sample_ID == rownames(meth.t4))/nrow(meth.t4) #should be 1
all.equal(as.character(pheno4$Sample_ID),rownames(meth.t4))  # must be TRUE

## Step f) Z-score standize all CpGs before run glmnet model
library(tibble)
library(psych)

# standardize cpg sites 
meth.t <- as.data.frame(meth.t)
meth.z <- scale(meth.t)
meth.z[1:5,1:5] # check if it went well 

meth.t2 <- as.data.frame(meth.t2)
meth.z_tri1 <- scale(meth.t2)
meth.z_tri1[1:5,1:5] 

meth.t3 <- as.data.frame(meth.t3)
meth.z_tri2 <- scale(meth.t3)
meth.z_tri2[1:5,1:5] 

meth.t4 <- as.data.frame(meth.t4)
meth.z_tri3 <- scale(meth.t4)
meth.z_tri3[1:5,1:5] 

# transform to dataframe 
meth.z2 <- as.data.frame(meth.z)
meth.z3 <-  tibble::rownames_to_column(meth.z2, "Sample_ID")

meth.z_tri1 <- as.data.frame(meth.z_tri1)
meth.z_tri1 <-  tibble::rownames_to_column(meth.z_tri1, "Sample_ID")

meth.z_tri2 <- as.data.frame(meth.z_tri2)
meth.z_tri2 <-  tibble::rownames_to_column(meth.z_tri2, "Sample_ID")

meth.z_tri3 <- as.data.frame(meth.z_tri3)
meth.z_tri3 <-  tibble::rownames_to_column(meth.z_tri3, "Sample_ID")

## Step g) Merge to two data for sample splitting next
identical(pheno$Sample_ID,meth.z3$Sample_ID) #should be true
df <-merge(pheno,meth.z3,by = "Sample_ID")
str(df)

identical(pheno2$Sample_ID,meth.z_tri1$Sample_ID) #should be true
df_tri1 <-merge(pheno2,meth.z_tri1,by = "Sample_ID")

identical(pheno3$Sample_ID,meth.z_tri2$Sample_ID) #should be true
df_tri2 <-merge(pheno3,meth.z_tri2,by = "Sample_ID")

identical(pheno4$Sample_ID,meth.z_tri3$Sample_ID) #should be true
df_tri3 <-merge(pheno4,meth.z_tri3,by = "Sample_ID")

#\ 

###----------------STEP 3: Splitting training and test subsets in GenR----------------###

## We do not want to split the sample at random, but we want to have an equal distribution of samples for assay across both sets
# @Anna, you need to creat a binary variable "array" (450K=1,Epic =2),set as.factor. then you will have similar distribution for both arrays in train and test.

## Step a) split the data into training and test set

# First set seed; important for reproducibility 
set.seed(2023) 

# CreateDataPartition will give equal distribution of samples based on arrays when splitting your data in an 80/20 train/test set. 
trainIndex <- createDataPartition(df$assay,
                                  p = 0.8, 
                                  list = FALSE)  

set.seed(2023) 
trainIndex_tri1 <- createDataPartition(df_tri1$assay,
                                  p = 0.8, 
                                  list = FALSE)  

set.seed(2023) 
trainIndex_tri2 <- createDataPartition(df_tri2$assay,
                                  p = 0.8, 
                                  list = FALSE)  

set.seed(2023) 
trainIndex_tri3 <- createDataPartition(df_tri3$assay,
                                  p = 0.8, 
                                  list = FALSE)  

# Create train set for model estimation #(n=1871)
datTrain <- df[trainIndex, ]   
str(datTrain)
dim(datTrain) 

datTrain_tri1 <- df_tri1[trainIndex_tri1, ]   
str(datTrain_tri1)
dim(datTrain_tri1)

datTrain_tri2 <- df_tri2[trainIndex_tri2, ]   
str(datTrain_tri2)
dim(datTrain_tri2)

datTrain_tri3 <- df_tri3[trainIndex_tri3, ]   
str(datTrain_tri3)
dim(datTrain_tri3)

# Create test set for model validation (n=467)
datTest <- df[-trainIndex, ] 
str(datTest)
dim(datTest)

datTest_tri1 <- df_tri1[-trainIndex_tri1, ] 
str(datTest_tri1)
dim(datTest_tri1)

datTest_tri2 <- df_tri2[-trainIndex_tri2, ] 
str(datTest_tri2)
dim(datTest_tri2)

datTest_tri3 <- df_tri3[-trainIndex_tri3, ] 
str(datTest_tri3)
dim(datTest_tri3)

# save datasets 
setwd('path_to_train_sets')

saveRDS(datTrain,"train_totalinf.rds")
saveRDS(datTest,"test_totalinf.rds")

saveRDS(datTrain_tri1,"train_tri1.rds")
saveRDS(datTest_tri1,"test_tri1.rds")

saveRDS(datTrain_tri2,"train_tri2.rds")
saveRDS(datTest_tri2,"test_tri2.rds")

saveRDS(datTrain_tri3,"train_tri3.rds")
saveRDS(datTest_tri3,"test_tri3.rds")

# Removing the sample id, only keep predictors and outcome in data
row.names(datTrain) <- datTrain$Sample_ID
datTrain[1] <- NULL
datTrain[1:5,1:5]
dim(datTrain)

row.names(datTest) <- datTest$Sample_ID
datTest[1] <- NULL
datTest[1:5,1:5]
dim(datTest)

row.names(datTrain_tri1) <- datTrain_tri1$Sample_ID
datTrain_tri1[1] <- NULL
datTrain_tri1[1:5,1:5]
dim(datTrain_tri1)

row.names(datTest_tri1) <- datTest_tri1$Sample_ID
datTest_tri1[1] <- NULL
datTest_tri1[1:5,1:5]
dim(datTest_tri1)

row.names(datTrain_tri2) <- datTrain_tri2$Sample_ID
datTrain_tri2[1] <- NULL
datTrain_tri2[1:5,1:5]
dim(datTrain_tri2)

row.names(datTest_tri2) <- datTest_tri2$Sample_ID
datTest_tri2[1] <- NULL
datTest_tri2[1:5,1:5]
dim(datTest_tri2)

row.names(datTrain_tri3) <- datTrain_tri3$Sample_ID
datTrain_tri3[1] <- NULL
datTrain_tri3[1:5,1:5]
dim(datTrain_tri3)

row.names(datTest_tri3) <- datTest_tri3$Sample_ID
datTest_tri3[1] <- NULL
datTest_tri3[1:5,1:5]
dim(datTest_tri3)

#\ 

###----------------STEP 4: Applying cva.glmnet in training set----------------###
# Elastic net cross-validation for both alpha and lambda: 
# 10-Fold cross-validation to determine optimal combination of alpha and lambda
# https://search.r-project.org/CRAN/refmans/glmnetUtils/html/cva.glmnet.html

## Step a) check dataset; we only want to include variables that we will use for glmnet model and we want to remove any unnecessary variables 
names(datTrain)
str(datTrain)
datTrain2 <- dplyr::select(datTrain, -c('IDC', 'IDM', 'MOTHER', 'sumscore_inf_tri1', 'sumscore_inf_tri2', 'sumscore_inf_tri3', 'sumscore_inf_tot', 'sumscore_inf_tri1_standardized', 'sumscore_inf_tri2_standardized', 'sumscore_inf_tri3_standardized', 'GENDER', 'SMOKE_ALL', 'WEIGHT', 'GESTBIR', 'PARITY', 'EDUCM_3groups', 'ethnicity', 'Sample_Plate', 'Bcell', 'CD4T', 'CD8T', 'Gran', 'Mono', 'NK', 'nRBC', 'assay', 'AGE_M_v2'))

names(datTrain_tri1)
str(datTrain_tri1)
datTrain_tri1b <- dplyr::select(datTrain_tri1, -c('IDC', 'IDM', 'MOTHER', 'sumscore_inf_tri1', 'sumscore_inf_tri2', 'sumscore_inf_tri3', 'sumscore_inf_tot', 'sumscore_inf_tot_standardized', 'sumscore_inf_tri2_standardized', 'sumscore_inf_tri3_standardized', 'GENDER', 'SMOKE_ALL', 'WEIGHT', 'GESTBIR', 'PARITY', 'EDUCM_3groups', 'ethnicity', 'Sample_Plate', 'Bcell', 'CD4T', 'CD8T', 'Gran', 'Mono', 'NK', 'nRBC', 'assay','AGE_M_v2'))

names(datTrain_tri2)
str(datTrain_tri2)
datTrain_tri2b <- dplyr::select(datTrain_tri2, -c('IDC', 'IDM', 'MOTHER', 'sumscore_inf_tri1', 'sumscore_inf_tri2', 'sumscore_inf_tri3', 'sumscore_inf_tot', 'sumscore_inf_tot_standardized', 'sumscore_inf_tri1_standardized', 'sumscore_inf_tri3_standardized', 'GENDER', 'SMOKE_ALL', 'WEIGHT', 'GESTBIR', 'PARITY', 'EDUCM_3groups', 'ethnicity', 'Sample_Plate', 'Bcell', 'CD4T', 'CD8T', 'Gran', 'Mono', 'NK', 'nRBC', 'assay','AGE_M_v2'))

names(datTrain_tri3)
str(datTrain_tri3)
datTrain_tri3b <- dplyr::select(datTrain_tri3, -c('IDC', 'IDM', 'MOTHER', 'sumscore_inf_tri1', 'sumscore_inf_tri2', 'sumscore_inf_tri3', 'sumscore_inf_tot', 'sumscore_inf_tot_standardized', 'sumscore_inf_tri2_standardized', 'sumscore_inf_tri1_standardized', 'GENDER', 'SMOKE_ALL', 'WEIGHT', 'GESTBIR', 'PARITY', 'EDUCM_3groups', 'ethnicity', 'Sample_Plate', 'Bcell', 'CD4T', 'CD8T', 'Gran', 'Mono', 'NK', 'nRBC', 'assay','AGE_M_v2'))

## Step b) run glmnet with 10 cross-validation
set.seed(2023)
train.cvafit <- cva.glmnet(sumscore_inf_tot_standardized ~ ., data = datTrain2, family = "gaussian")

set.seed(2023)
train.cvafit_tri1 <- cva.glmnet(sumscore_inf_tri1_standardized ~ ., data = datTrain_tri1b, family = "gaussian")

set.seed(2023)
train.cvafit_tri2 <- cva.glmnet(sumscore_inf_tri2_standardized ~ ., data = datTrain_tri2b, family = "gaussian")

set.seed(2023)
train.cvafit_tri3 <- cva.glmnet(sumscore_inf_tri3_standardized ~ ., data = datTrain_tri3b, family = "gaussian")

# Save the plot object as a PNG file
setwd('path_to_output_files')

png("lambda_plot_totalinf.png")
plot <- plot(train.cvafit)
print(plot)     # Print the plot to the device
dev.off()       # Close the device and save the PNG file

png("lambda_plot_trimester1.png")
plot_tri1 <- plot(train.cvafit_tri1)
print(plot_tri1)     
dev.off()

png("lambda_plot_trimester2.png")
plot_tri2 <- plot(train.cvafit_tri2)
print(plot_tri2)     
dev.off() 

png("lambda_plot_trimester3.png")
plot_tri3 <- plot(train.cvafit_tri3)
print(plot_tri3)     
dev.off()   

## Step c) take a close look at object returned by cva.glmnet
print(train.cvafit)
class(train.cvafit)
str(train.cvafit)
train.cvafit$modlist[[1]] 
train.cvafit$modlist[[1]]$glmnet.fit
train.cvafit$alpha

print(train.cvafit_tri1)
class(train.cvafit_tri1)
str(train.cvafit_tri1)
train.cvafit_tri1$modlist[[1]] 
train.cvafit_tri1$modlist[[1]]$glmnet.fit
train.cvafit_tri1$alpha

print(train.cvafit_tri2)
class(train.cvafit_tri2)
str(train.cvafit_tri2)
train.cvafit_tri2$modlist[[1]] 
train.cvafit_tri2$modlist[[1]]$glmnet.fit
train.cvafit_tri2$alpha

print(train.cvafit_tri3)
class(train.cvafit_tri3)
str(train.cvafit_tri3)
train.cvafit_tri3$modlist[[1]] 
train.cvafit_tri3$modlist[[1]]$glmnet.fit
train.cvafit_tri3$alpha

# check optimal alpha
get_alpha <- function(fit) {
  alpha <- fit$alpha
  error <- sapply(fit$modlist, function(mod) {min(mod$cvm)})
  alpha[which.min(error)]
}

#cvm. The mean cross-validated error - a vector of length length(lambda).
# so this function helps to pick up the alpha with lowest mean CV error
get_alpha(train.cvafit) # In this case that is 0.027 
get_alpha(train.cvafit_tri1) #In this case that is 1
get_alpha(train.cvafit_tri2) #in this case that is 0.001
get_alpha(train.cvafit_tri3) #In this case that is 0.008

# Get all parameters
get_model_params <- function(fit) {
  alpha <- fit$alpha
  lambdaMin <- sapply(fit$modlist, `[[`, "lambda.min")
  lambdaSE <- sapply(fit$modlist, `[[`, "lambda.1se")
  error <- sapply(fit$modlist, function(mod) {min(mod$cvm)})
  best <- which.min(error)
  data.frame(alpha = alpha[best], lambdaMin = lambdaMin[best],
             lambdaSE = lambdaSE[best], eror = error[best])
}

# alpha lambdaMin  lambdaSE  eror
#lambda.min gives minimum mean cross-validated error, 
# while lambda.1se gives the most regularized model such that the cross-validated error is within one standard error of the minimum

get_model_params(train.cvafit) # in this case, alpha = 0.027, lambdaMin = 3.224541, lambdaSE = 5.134381, error = 1.00824
get_model_params(train.cvafit_tri1) #alpha = 1, lambdaMin = 0.141458, lambdaSE = 0.1481936, error = 0.9730963
get_model_params(train.cvafit_tri2) #alpha = 0.001, lambdaMin = 133.737, lambdaSE = 133.737, error = 0.9903397
get_model_params(train.cvafit_tri3) #alpha = 0.008, lambdaMin = 7.251467, lambdaSE = 20.17764, error = 1.011117

# Extract non-zero coefficients at optimal alpha and lambda that you get from alpha & lambdaMin from above step (so these have the best cross validated performance)

coef <-coef(train.cvafit,alpha = 0.027, s = 3.224541)
coef_tri1 <-coef(train.cvafit_tri1,alpha = 1, s = 0.141458)
coef_tri2 <-coef(train.cvafit_tri2,alpha = 0.001, s = 133.737)
coef_tri3 <-coef(train.cvafit_tri3,alpha = 0.008, s = 7.251467)

# Removing predictors with zero coef.
coef[coef == 0] <- NA
coef <-as.matrix(coef)
coef <-as.data.frame(coef)
str(coef)
dim(coef) 
coef <-  tibble::rownames_to_column(coef, "probeid")
coef$beta <- coef$s1
coef$s1 <- NULL
coef <-  coef[complete.cases(coef), ] 
str(coef)

coef_tri1[coef_tri1 == 0] <- NA
coef_tri1 <-as.matrix(coef_tri1)
coef_tri1 <-as.data.frame(coef_tri1)
str(coef_tri1)
dim(coef_tri1) 
coef_tri1 <-  tibble::rownames_to_column(coef_tri1, "probeid")
coef_tri1$beta <- coef_tri1$s1
coef_tri1$s1 <- NULL
coef_tri1 <-  coef_tri1[complete.cases(coef_tri1), ] 
str(coef_tri1)

coef_tri2[coef_tri2 == 0] <- NA
coef_tri2 <-as.matrix(coef_tri2)
coef_tri2 <-as.data.frame(coef_tri2)
str(coef_tri2)
dim(coef_tri2) 
coef_tri2 <-  tibble::rownames_to_column(coef_tri2, "probeid")
coef_tri2$beta <- coef_tri2$s1
coef_tri2$s1 <- NULL
coef_tri2 <-  coef_tri2[complete.cases(coef_tri2), ] 
str(coef_tri2)

coef_tri3[coef_tri3 == 0] <- NA
coef_tri3 <-as.matrix(coef_tri3)
coef_tri3 <-as.data.frame(coef_tri3)
str(coef_tri3)
dim(coef_tri3) 
coef_tri3 <-  tibble::rownames_to_column(coef_tri3, "probeid")
coef_tri3$beta <- coef_tri3$s1
coef_tri3$s1 <- NULL
coef_tri3 <-  coef_tri3[complete.cases(coef_tri3), ] 
str(coef_tri3)

# Removing intercept, then save coef as coef.cpg (we end with 133 cpg sites)
coef.cpg_totalinf <- coef[coef$probeid !='(Intercept)',]
coef.cpg_tri1 <- coef_tri1[coef_tri1$probeid !='(Intercept)',]
coef.cpg_tri2 <- coef_tri2[coef_tri2$probeid !='(Intercept)',]
coef.cpg_tri3 <- coef_tri3[coef_tri3$probeid !='(Intercept)',]

write.csv(coef.cpg_totalinf, 'ENR_weights_total_infection_GenR.csv', row.names=F, quote=F)
write.csv(coef.cpg_tri1, 'ENR_weights_trimester1_GenR.csv', row.names=F, quote=F)
write.csv(coef.cpg_tri2, 'ENR_weights_trimester2_GenR.csv', row.names=F, quote=F)
write.csv(coef.cpg_tri3, 'ENR_weights_trimester3_GenR.csv', row.names=F, quote=F)

library(writexl)
write_xlsx(coef.cpg_totalinf, 'ENR_weights_total_infection_GenR.xlsx')
write_xlsx(coef.cpg_tri1, 'ENR_weights_trimester1_GenR.xlsx')
write_xlsx(coef.cpg_tri2, 'ENR_weights_trimester2_GenR.xlsx')
write_xlsx(coef.cpg_tri3, 'ENR_weights_trimester3_GenR.xlsx')

# Examine the shrinking of coefficients if you like
best_apha.fit  <- glmnet(sumscore_inf_tot_standardized ~ ., data = datTrain2, family = "gaussian")
print(best_apha.fit)
png("best_alpha_fit_plot_totalinf.png")
plot2 <- plot(best_apha.fit)
print(plot2)     
dev.off()       

best_apha.fit_tri1  <- glmnet(sumscore_inf_tri1_standardized ~ ., data = datTrain_tri1b, family = "gaussian")
print(best_apha.fit_tri1)
png("best_alpha_fit_plot_trimester1.png")
plot_alpha_tri1 <- plot(best_apha.fit_tri1)
print(plot_alpha_tri1)     
dev.off()       

best_apha.fit_tri2  <- glmnet(sumscore_inf_tri2_standardized ~ ., data = datTrain_tri2b, family = "gaussian")
print(best_apha.fit_tri2)
png("best_alpha_fit_plot_trimester2.png")
plot_alpha_tri2 <- plot(best_apha.fit_tri2)
print(plot_alpha_tri2)     
dev.off()

best_apha.fit_tri3  <- glmnet(sumscore_inf_tri3_standardized ~ ., data = datTrain_tri3b, family = "gaussian")
print(best_apha.fit_tri3)
png("best_alpha_fit_plot_trimester3.png")
plot_alpha_tri3 <- plot(best_apha.fit_tri3)
print(plot_alpha_tri3)     
dev.off()    

#\ 

###----------------STEP 5: Constructing methylation profile score----------------###

rm(list=ls())

## Step a) load re-weighted effect size of cpGs
setwd('path_to_output_files')

library(readxl)

ENR_w_totalinf <- read_excel('ENR_weights_total_infection_GenR.xlsx') # yields 133 cpg sites 
ENR_W_tri1 <- read_excel('ENR_weights_trimester1_GenR.xlsx') # yields 1 cpg site 
ENR_W_tri2 <- read_excel('ENR_weights_trimester2_GenR.xlsx') # yields 0 cpg sites 
ENR_W_tri3 <- read_excel('ENR_weights_trimester3_GenR.xlsx') # yields 996 cpg sites 

ENR_w_totalinf<-ENR_w_totalinf[order(ENR_w_totalinf$probeid),]
ENR_W_tri3<-ENR_W_tri3[order(ENR_W_tri3$probeid),]

# Load & check the testing set
setwd('path_to_test_sets')

datTest_totalinf <- readRDS('test_totalinf.rds')
datTest_tri3 <- readRDS('test_tri3.rds')

str(datTest_totalinf)
str(datTest_tri3)

## Step b) transpose the meth data in datTest

# Only select CpGs, excluded phenotypes 
meth.test.totalinf <- dplyr::select(datTest_totalinf, -c('Sample_ID', 'IDC', 'IDM', 'MOTHER', 'sumscore_inf_tri1', 'sumscore_inf_tri2', 'sumscore_inf_tri3', 'sumscore_inf_tot', 'sumscore_inf_tot_standardized', 'sumscore_inf_tri2_standardized', 'sumscore_inf_tri1_standardized', 'sumscore_inf_tri3_standardized','GENDER', 'SMOKE_ALL', 'WEIGHT', 'GESTBIR', 'PARITY', 'EDUCM_3groups', 'ethnicity', 'Sample_Plate', 'Bcell', 'CD4T', 'CD8T', 'Gran', 'Mono', 'NK', 'nRBC', 'assay','AGE_M_v2')) 

meth.test.tri3 <- dplyr::select(datTest_tri3, -c('Sample_ID', 'IDC', 'IDM', 'MOTHER', 'sumscore_inf_tri1', 'sumscore_inf_tri2', 'sumscore_inf_tri3', 'sumscore_inf_tot', 'sumscore_inf_tot_standardized', 'sumscore_inf_tri2_standardized', 'sumscore_inf_tri1_standardized', 'sumscore_inf_tri3_standardized','GENDER', 'SMOKE_ALL', 'WEIGHT', 'GESTBIR', 'PARITY', 'EDUCM_3groups', 'ethnicity', 'Sample_Plate', 'Bcell', 'CD4T', 'CD8T', 'Gran', 'Mono', 'NK', 'nRBC', 'assay','AGE_M_v2')) 

# Transpose dataset 
meth.test.totalinf.t  <-t(meth.test.totalinf)
meth.test.tri3.t <- t(meth.test.tri3)

# Check dimensions of dataset 
meth.test.totalinf.t[1:5,1:5]
dim(meth.test.totalinf.t)

meth.test.tri3.t[1:5,1:5]
dim(meth.test.tri3.t)

## Step c) match DNAm data with ENR-based weight data, based on probeID(i.e.CpG name)
# only include identical CpGs between two files
methsub <- meth.test.totalinf.t[na.omit(match(ENR_w_totalinf$probeid,rownames(meth.test.totalinf.t))),]
dim(methsub)
ENR_w_totalinf2 <- ENR_w_totalinf[match(ENR_w_totalinf$probeid,rownames(methsub)),]
dim(ENR_w_totalinf2)

methsub_tri3 <- meth.test.tri3.t[na.omit(match(ENR_W_tri3$probeid,rownames(meth.test.tri3.t))),]
dim(methsub_tri3)
ENR_W_tri3 <- ENR_W_tri3[match(ENR_W_tri3$probeid,rownames(methsub_tri3)),]
dim(ENR_W_tri3)

# Check that all rows are the same - need to be in the same order
identical(rownames(methsub), ENR_w_totalinf2$probeid) #needs to be TRUE
methsub[1:5,1:5]
ENR_w_totalinf2[1:5,]

identical(rownames(methsub_tri3), ENR_W_tri3$probeid) #needs to be TRUE
methsub_tri3[1:5,1:5]
ENR_W_tri3[1:5,]

## Step d) Matrix multiplication (meth*weight)
# This multiples each CpG by its weight(ie.effect size from ENR) 
wmatrix = (methsub * ENR_w_totalinf2$beta) 
dim(wmatrix)

wmatrix3 = (methsub_tri3 * ENR_W_tri3$beta) 
dim(wmatrix3)

# To double check 
wmatrix[1:5,1:5]
methsub[1:5,1:5]
ENR_w_totalinf2[1:5,]
colsum = colSums(wmatrix, na.rm=TRUE)
mps = as.data.frame(colsum)
head(mps)

wmatrix3[1:5,1:5]
methsub_tri3[1:5,1:5]
ENR_W_tri3[1:5,]
colsum3 = colSums(wmatrix3, na.rm=TRUE)
mps_tri3 = as.data.frame(colsum3)
head(mps_tri3)

# To convert row names into first column
mps <- tibble::rownames_to_column(mps, "Sample_ID")
colnames(mps)[2] <- "MPS_total_infection"
head(mps)

mps_tri3 <- tibble::rownames_to_column(mps_tri3, "Sample_ID")
colnames(mps_tri3)[2] <- "MPS_trimester3_infection"
head(mps_tri3)

## Step e) Save mps files 
setwd('path_to_results_MPS')

library(writexl)

write_xlsx(mps, 'MPS_total_infection_GenR_testset.xlsx')
write_xlsx(mps_tri3, 'MPS_trimester3_infection_GenR_testset.xlsx')

#\

###----------------STEP 6: Assessing model performance in test set----------------###

# Run this part on local r studio

## Step a) combine MPS file with phenotype file 
# Load phenotype and MPS data files
setwd("path_to_test_files")
datTest_totalinf <- readRDS('test_totalinf.rds')
datTest_tri3 <- readRDS('test_tri3.rds')

setwd("path_to_results_MPS")
library(readxl)
mps <- read_excel('MPS_total_infection_GenR_testset.xlsx')
mps_tri3 <- read_excel('MPS_trimester3_infection_GenR_testset.xlsx')

# Subset the phenotype data; exclude CpGs
phenodat.totalinf <- dplyr::select(datTest_totalinf, c('IDC', 'IDM', 'MOTHER', 'sumscore_inf_tri1', 'sumscore_inf_tri2', 'sumscore_inf_tri3', 'sumscore_inf_tot', 'sumscore_inf_tot_standardized', 'sumscore_inf_tri2_standardized', 'sumscore_inf_tri1_standardized', 'sumscore_inf_tri3_standardized','GENDER', 'SMOKE_ALL', 'WEIGHT', 'GESTBIR', 'PARITY', 'EDUCM_3groups', 'ethnicity', 'Sample_Plate', 'Bcell', 'CD4T', 'CD8T', 'Gran', 'Mono', 'NK', 'nRBC', 'assay','AGE_M_v2')) 
pheno.totalinf <- tibble::rownames_to_column(phenodat.totalinf, "Sample_ID")
head(pheno.totalinf)

phenodat.trimester3 <- dplyr::select(datTest_tri3, c('IDC', 'IDM', 'MOTHER', 'sumscore_inf_tri1', 'sumscore_inf_tri2', 'sumscore_inf_tri3', 'sumscore_inf_tot', 'sumscore_inf_tot_standardized', 'sumscore_inf_tri2_standardized', 'sumscore_inf_tri1_standardized', 'sumscore_inf_tri3_standardized','GENDER', 'SMOKE_ALL', 'WEIGHT', 'GESTBIR', 'PARITY', 'EDUCM_3groups', 'ethnicity', 'Sample_Plate', 'Bcell', 'CD4T', 'CD8T', 'Gran', 'Mono', 'NK', 'nRBC', 'assay','AGE_M_v2')) 
pheno.trimester3 <- tibble::rownames_to_column(phenodat.trimester3, "Sample_ID")
head(pheno.trimester3)

# Merge testing phenotype data and MPS file
head(mps)
df_totalinfections <- merge(pheno.totalinf,mps, by = "Sample_ID")

head(mps_tri3)
df_trimester3 <- merge(pheno.trimester3,mps_tri3, by = "Sample_ID")

## Step b) associate prenatal infections with MPS of infections 
# b1. Linear regression 
setwd("path_to_MPS_validation")
lm_totalinf <- summary(lm(scale(MPS_total_infection) ~  sumscore_inf_tot_standardized+ GENDER + SMOKE_ALL + AGE_M_v2 + PARITY + EDUCM_3groups + GESTBIR + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC, data =df_totalinfections))$coefficients[2,] 

lm_tri3 <- summary(lm(scale(MPS_trimester3_infection) ~  sumscore_inf_tri3_standardized+ GENDER + SMOKE_ALL + AGE_M_v2 + PARITY + EDUCM_3groups + GESTBIR + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC, data =df_trimester3))$coefficients[2,] 

library(tibble) # create df with results to later save to excel 
lm_total <- rbind(lm_totalinf, lm_tri3)
lm_total <- as.data.frame(lm_total)
lm_total <-  tibble::rownames_to_column(lm_total, "Exposure")

# b2. Calculate incremental r2 --> this indicates prediction peformance of MPS 
# mod 1=basic modle, including all covariates in your EWAS
mod1 <- lm(sumscore_inf_tot_standardized ~  GENDER + SMOKE_ALL + AGE_M_v2 + PARITY + EDUCM_3groups + GESTBIR + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC, df_totalinfections)
r.mod1 = summary(mod1)$r.squared   

mod2 <- lm(sumscore_inf_tot_standardized ~ scale(MPS_total_infection) +  GENDER + SMOKE_ALL + AGE_M_v2 + PARITY + EDUCM_3groups + GESTBIR + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC, df_totalinfections)
r.mod2 = summary(mod2)$r.squared
r2change =  r.mod2 - r.mod1 #0.0001581543

mod3 <- lm(sumscore_inf_tri3_standardized ~  GENDER + SMOKE_ALL + AGE_M_v2 + PARITY + EDUCM_3groups + GESTBIR + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC, df_trimester3)
r.mod3 = summary(mod3)$r.squared   

mod4 <- lm(sumscore_inf_tri3_standardized ~ scale(MPS_trimester3_infection) +  GENDER + SMOKE_ALL + AGE_M_v2 + PARITY + EDUCM_3groups + GESTBIR + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC, df_trimester3)
r.mod4 = summary(mod4)$r.squared
r2change2 =  r.mod4 - r.mod3 #0.001341319

lm_total$r2_change <- c(0.0001581543,0.001341319) # add results to lm validation df 

# b3. ROC curves & calculating AUC to determine discriminative ability 
library(pROC)

resid<- lm(MPS_total_infection~GENDER + SMOKE_ALL + AGE_M_v2 + PARITY + EDUCM_3groups + GESTBIR + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC, data =df_totalinfections) 
resid2<- lm(MPS_trimester3_infection ~ GENDER + SMOKE_ALL + AGE_M_v2 + PARITY + EDUCM_3groups + GESTBIR + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC, data =df_trimester3) 

df_totalinfections$MPS_total_infection_resid <- residuals(resid)
df_trimester3$MPS_trimester3_infection_resid <- residuals(resid2)

mean(df_totalinfections$sumscore_inf_tot_standardized) # 0.013
sd(df_totalinfections$sumscore_inf_tot_standardized) # 0.986
mean(df_trimester3$sumscore_inf_tri3_standardized) # 0.025
sd(df_trimester3$sumscore_inf_tri3_standardized) # 0.990

df_totalinfections$infections_categorical <- as.factor(ifelse(df_totalinfections$sumscore_inf_tot_standardized > 0.999, 1, 0)) #Categorize infections (based on +1 SD of mean)
df_trimester3$infections_categorical <- as.factor(ifelse(df_trimester3$sumscore_inf_tri3_standardized > 1.016, 1, 0)) #Categorize infections (based on +1 SD of mean)

roc <- roc(infections_categorical ~ MPS_total_infection_resid, data = df_totalinfections) #auc = 0.51
roc <- roc(infections_categorical ~ MPS_trimester3_infection_resid, data = df_trimester3) #auc = 0.49

library(writexl)
lm_total$auc <- c(0.5101, 0.4871)
write_xlsx(lm_total, 'lm_infection_and_mps.xlsx') # save results of lm and r2 to excel 

png("roc_total_infections.png")
roc_plot<- plot(roc)
print(roc_plot)     
dev.off()

png("roc_trimester3_infections.png")
roc_plot2<- plot(roc)
print(roc_plot2)     
dev.off()

# b4. Associating MPS of infections to child outcomes
setwd("path_to_MPS_validation_data")

# First load and combine child pheno data
library(foreign)
library(dplyr)
second_hits_df <- read.spss('Second_hits_DF.sav', to.data.frame = T) 
second_hits_df <- dplyr::select(second_hits_df, c('IDC', 'IDM', 'sum_int_14', 'sum_ext_14', 'cbcl_sum_14'))
crp_cordblood <- read.spss('CHILDCORDBLOOD-CRP_23092013.sav', to.data.frame = T)
crp_cordblood <- dplyr::select(crp_cordblood, c('IDC', "CRP_birth"))
crp_child5 <- read.spss('CHILDCRP5_24092013.sav', to.data.frame = T)
crp_child5 <- dplyr::select(crp_child5, c('IDC', "CRPCHILD5"))

child.pheno <- merge(second_hits_df, crp_cordblood, by = 'IDC', all.x =  T)
child.pheno <- merge(child.pheno, crp_child5, by = 'IDC', all.x =  T)

# Combine child pheno with mps data 
child.mps.pheno <- merge(df_totalinfections, child.pheno, by = 'IDC', all.x = T)
child.mps_tri3.pheno <- merge(df_trimester3, child.pheno, by = 'IDC', all.x = T)

# Run regressions
child_pheno <- c('sum_int_14', 'sum_ext_14', 'cbcl_sum_14', 'CRP_birth', 'CRPCHILD5') # create vector with all outcomes to loop over 

results_df <- data.frame() # Create an empty data frame to store results

for(x in child_pheno) {
  # lm formula
  a <- paste0(x, "~ scale(MPS_total_infection) + GENDER + SMOKE_ALL + AGE_M_v2 + PARITY + EDUCM_3groups + GESTBIR + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC")
  
  b <- paste0(x, "~ scale(MPS_trimester3_infection) + GENDER + SMOKE_ALL + AGE_M_v2 + PARITY + EDUCM_3groups + GESTBIR + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC")
  
  #extract coefficients 
  coef_totalinf <- summary(lm(as.formula(a), data = child.mps.pheno))$coefficients[2,]
  coef_tri3inf <- summary(lm(as.formula(b), data = child.mps_tri3.pheno))$coefficients[2,]
  
  #Create a data frame with the results
  result_row <- data.frame(
    Phenotype = x,
    Total_Infection_Coefficient = coef_totalinf,
    Tri3_Infection_Coefficient = coef_tri3inf
  )
  
  # Append the results to the main data frame
  results_df <- rbind(results_df, result_row)
}

# Write the results to Excel
setwd("path_to_MPS_validation_results")

results_df2 <- t(results_df)
results_df2 <- as.data.frame(results_df2)

write_xlsx(results_df2, path = "child_pheno_valid_mps_results.xlsx")

#\ END OF SCRIPT. 
