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

# Of note, in the prior script we used a preselection threshold of nominal significance and the calibration/internal validation of the score was not good. Hence, we now make a score using suggestive significance p value as threshold and input those CpG sites into ENR. 

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
setwd('path_to_data')

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

# create vectors for suggestive sig cpgs
sigCpGs_totalinf <- res_totalinf$cpg[res_totalinf$meta_pval<5e-5]
length(sigCpGs_totalinf) #33
sigCpGs_trimester1 <- res_trim1$cpg[res_trim1$meta_pval<5e-5]
length(sigCpGs_trimester1) #20
sigCpGs_trimester2 <- res_trim2$cpg[res_trim2$meta_pval<5e-5]
length(sigCpGs_trimester2) #24
sigCpGs_trimester3 <- res_trim3$cpg[res_trim3$meta_pval<5e-5]
length(sigCpGs_trimester3)#26

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
setwd('path_to_results')

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
get_alpha(train.cvafit) # In this case that is 0.729
get_alpha(train.cvafit_tri1) #In this case that is 0.001
get_alpha(train.cvafit_tri2) #in this case that is 0.001
get_alpha(train.cvafit_tri3) #In this case that is 0.125

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

get_model_params(train.cvafit) # in this case, alpha = 0.729, lambdaMin = 0.0004544753, lambdaSE = 0.05734828, error = 0.980779
get_model_params(train.cvafit_tri1) #alpha = 0.001, lambdaMin = 0.01944439, lambdaSE = 18.99736, error = 0.99672
get_model_params(train.cvafit_tri2) #alpha = 0.001, lambdaMin = 0.08742814, lambdaSE = 48.87941, error = 0.9807181
get_model_params(train.cvafit_tri3) #alpha = 0.125, lambdaMin = 0.0009160025, lambdaSE = 0.6168499, error = 0.9880418

# Extract non-zero coefficients at optimal alpha and lambda that you get from alpha & lambdaMin from above step (so these have the best cross validated performance)

coef <-coef(train.cvafit,alpha = 0.729, s = 0.0004544753)
coef_tri1 <-coef(train.cvafit_tri1,alpha = 0.001, s = 0.01944439)
coef_tri2 <-coef(train.cvafit_tri2,alpha = 0.001, s = 0.08742814)
coef_tri3 <-coef(train.cvafit_tri3,alpha = 0.125, s = 0.0009160025)

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

write.csv(coef.cpg_totalinf, 'ENR_weights_total_infection_GenR_suggesivepval.csv', row.names=F, quote=F)
write.csv(coef.cpg_tri1, 'ENR_weights_trimester1_GenR_suggesivepval.csv', row.names=F, quote=F)
write.csv(coef.cpg_tri2, 'ENR_weights_trimester2_GenR_suggesivepval.csv', row.names=F, quote=F)
write.csv(coef.cpg_tri3, 'ENR_weights_trimester3_GenR_suggesivepval.csv', row.names=F, quote=F)

library(writexl)
write_xlsx(coef.cpg_totalinf, 'ENR_weights_total_infection_GenR_suggesivepval.xlsx')
write_xlsx(coef.cpg_tri1, 'ENR_weights_trimester1_GenR_suggesivepval.xlsx')
write_xlsx(coef.cpg_tri2, 'ENR_weights_trimester2_GenR_suggesivepval.xlsx')
write_xlsx(coef.cpg_tri3, 'ENR_weights_trimester3_GenR_suggesivepval.xlsx')

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

## Step a) load re-weighted effect size of cpGs
setwd('path_to_data')

library(readxl)

ENR_w_totalinf <- read_excel('ENR_weights_total_infection_GenR_suggesivepval.xlsx') # yields 30 cpg sites 
ENR_W_tri1 <- read_excel('ENR_weights_trimester1_GenR_suggesivepval.xlsx') # yields 120 cpg site 
ENR_W_tri2 <- read_excel('ENR_weights_trimester2_GenR_suggesivepval.xlsx') # yields 24 cpg sites 
ENR_W_tri3 <- read_excel('ENR_weights_trimester3_GenR_suggesivepval.xlsx') # yields 26 cpg sites 

ENR_w_totalinf<-ENR_w_totalinf[order(ENR_w_totalinf$probeid),]
ENR_W_tri1<-ENR_W_tri1[order(ENR_W_tri1$probeid),]
ENR_W_tri2<-ENR_W_tri2[order(ENR_W_tri2$probeid),]
ENR_W_tri3<-ENR_W_tri3[order(ENR_W_tri3$probeid),]

## Step b) transpose the meth data in datTest

# Only select CpGs, excluded phenotypes 
meth.test.totalinf <- dplyr::select(datTest, -c('IDC', 'IDM', 'MOTHER', 'sumscore_inf_tri1', 'sumscore_inf_tri2', 'sumscore_inf_tri3', 'sumscore_inf_tot', 'sumscore_inf_tot_standardized', 'sumscore_inf_tri2_standardized', 'sumscore_inf_tri1_standardized', 'sumscore_inf_tri3_standardized','GENDER', 'SMOKE_ALL', 'WEIGHT', 'GESTBIR', 'PARITY', 'EDUCM_3groups', 'ethnicity', 'Sample_Plate', 'Bcell', 'CD4T', 'CD8T', 'Gran', 'Mono', 'NK', 'nRBC', 'assay','AGE_M_v2')) 

meth.test.tri1 <- dplyr::select(datTest_tri1, -c('IDC', 'IDM', 'MOTHER', 'sumscore_inf_tri1', 'sumscore_inf_tri2', 'sumscore_inf_tri3', 'sumscore_inf_tot', 'sumscore_inf_tot_standardized', 'sumscore_inf_tri2_standardized', 'sumscore_inf_tri1_standardized', 'sumscore_inf_tri3_standardized','GENDER', 'SMOKE_ALL', 'WEIGHT', 'GESTBIR', 'PARITY', 'EDUCM_3groups', 'ethnicity', 'Sample_Plate', 'Bcell', 'CD4T', 'CD8T', 'Gran', 'Mono', 'NK', 'nRBC', 'assay','AGE_M_v2')) 

meth.test.tri2 <- dplyr::select(datTest_tri2, -c('IDC', 'IDM', 'MOTHER', 'sumscore_inf_tri1', 'sumscore_inf_tri2', 'sumscore_inf_tri3', 'sumscore_inf_tot', 'sumscore_inf_tot_standardized', 'sumscore_inf_tri2_standardized', 'sumscore_inf_tri1_standardized', 'sumscore_inf_tri3_standardized','GENDER', 'SMOKE_ALL', 'WEIGHT', 'GESTBIR', 'PARITY', 'EDUCM_3groups', 'ethnicity', 'Sample_Plate', 'Bcell', 'CD4T', 'CD8T', 'Gran', 'Mono', 'NK', 'nRBC', 'assay','AGE_M_v2')) 

meth.test.tri3 <- dplyr::select(datTest_tri3, -c('IDC', 'IDM', 'MOTHER', 'sumscore_inf_tri1', 'sumscore_inf_tri2', 'sumscore_inf_tri3', 'sumscore_inf_tot', 'sumscore_inf_tot_standardized', 'sumscore_inf_tri2_standardized', 'sumscore_inf_tri1_standardized', 'sumscore_inf_tri3_standardized','GENDER', 'SMOKE_ALL', 'WEIGHT', 'GESTBIR', 'PARITY', 'EDUCM_3groups', 'ethnicity', 'Sample_Plate', 'Bcell', 'CD4T', 'CD8T', 'Gran', 'Mono', 'NK', 'nRBC', 'assay','AGE_M_v2')) 

# Transpose dataset 
meth.test.totalinf.t  <-t(meth.test.totalinf)
meth.test.tri1.t <- t(meth.test.tri1)
meth.test.tri2.t <- t(meth.test.tri2)
meth.test.tri3.t <- t(meth.test.tri3)

## Step c) match DNAm data with ENR-based weight data, based on probeID(i.e.CpG name)
# only include identical CpGs between two files
methsub <- meth.test.totalinf.t[na.omit(match(ENR_w_totalinf$probeid,rownames(meth.test.totalinf.t))),]
dim(methsub)
ENR_w_totalinf2 <- ENR_w_totalinf[match(ENR_w_totalinf$probeid,rownames(methsub)),]
dim(ENR_w_totalinf2)

methsub_tri1 <- meth.test.tri1.t[na.omit(match(ENR_W_tri1$probeid,rownames(meth.test.tri1.t))),]
dim(methsub_tri1)
ENR_W_tri1 <- ENR_W_tri1[match(ENR_W_tri1$probeid,rownames(methsub_tri1)),]
dim(ENR_W_tri1)

methsub_tri2 <- meth.test.tri2.t[na.omit(match(ENR_W_tri2$probeid,rownames(meth.test.tri2.t))),]
dim(methsub_tri2)
ENR_W_tri2 <- ENR_W_tri2[match(ENR_W_tri2$probeid,rownames(methsub_tri2)),]
dim(ENR_W_tri2)

methsub_tri3 <- meth.test.tri3.t[na.omit(match(ENR_W_tri3$probeid,rownames(meth.test.tri3.t))),]
dim(methsub_tri3)
ENR_W_tri3 <- ENR_W_tri3[match(ENR_W_tri3$probeid,rownames(methsub_tri3)),]
dim(ENR_W_tri3)

# Check that all rows are the same - need to be in the same order
identical(rownames(methsub), ENR_w_totalinf2$probeid) #needs to be TRUE
identical(rownames(methsub_tri1), ENR_W_tri1$probeid) 
identical(rownames(methsub_tri2), ENR_W_tri2$probeid) 
identical(rownames(methsub_tri3), ENR_W_tri3$probeid) 

## Step d) Matrix multiplication (meth*weight)
# This multiples each CpG by its weight(ie.effect size from ENR) 
wmatrix = (methsub * ENR_w_totalinf2$beta) 
dim(wmatrix)

wmatrix1 = (methsub_tri1 * ENR_W_tri1$beta) 
dim(wmatrix1)

wmatrix2 = (methsub_tri2 * ENR_W_tri2$beta) 
dim(wmatrix2)

wmatrix3 = (methsub_tri3 * ENR_W_tri3$beta) 
dim(wmatrix3)

# To double check 
wmatrix[1:5,1:5]
methsub[1:5,1:5]
ENR_w_totalinf2[1:5,]
colsum = colSums(wmatrix, na.rm=TRUE)
mps = as.data.frame(colsum)
head(mps)

wmatrix1[1:5,1:5]
methsub_tri1[1:5,1:5]
ENR_W_tri1[1:5,]
colsum1 = colSums(wmatrix1, na.rm=TRUE)
mps_tri1 = as.data.frame(colsum1)
head(mps_tri1)

wmatrix2[1:5,1:5]
methsub_tri2[1:5,1:5]
ENR_W_tri2[1:5,]
colsum2 = colSums(wmatrix2, na.rm=TRUE)
mps_tri2 = as.data.frame(colsum2)
head(mps_tri2)

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

mps_tri1 <- tibble::rownames_to_column(mps_tri1, "Sample_ID")
colnames(mps_tri1)[2] <- "MPS_trimester1_infection"
head(mps_tri1)

mps_tri2 <- tibble::rownames_to_column(mps_tri2, "Sample_ID")
colnames(mps_tri2)[2] <- "MPS_trimester2_infection"
head(mps_tri2)

mps_tri3 <- tibble::rownames_to_column(mps_tri3, "Sample_ID")
colnames(mps_tri3)[2] <- "MPS_trimester3_infection"
head(mps_tri3)

## Step e) Save mps files 
library(writexl)

write_xlsx(mps, 'MPS_total_infection_GenR_testset_suggestivepval.xlsx')
write_xlsx(mps_tri1, 'MPS_trimester1_infection_GenR_testset_suggestivepval.xlsx')
write_xlsx(mps_tri2, 'MPS_trimester2_infection_GenR_testset_suggestivepval.xlsx')
write_xlsx(mps_tri3, 'MPS_trimester3_infection_GenR_testset_suggestivepval.xlsx')

#\

###----------------STEP 6: Assessing model performance in test set----------------###

# Run this part on local r studio

## step 0) Test which cpgs in suggestive hits mps are also in nominal sig mps 
# not doing this for trimester 2 because no cpg sites in ENR with nominal significance 
setwd('path_to_data')
library(readxl)
ENR_weights_totalinf_nom_sig <- read_excel('ENR_weights_total_infection_GenR.xlsx')
ENR_weights_trimester1_nom_sig <- read_excel('ENR_weights_trimester1_GenR.xlsx')
ENR_weights_trimester3_nom_sig <- read_excel('ENR_weights_trimester3_GenR.xlsx')

setwd('path_to_results')
ENR_weights_totalinf_sug_sig <- read_excel('ENR_weights_total_infection_GenR_suggesivepval.xlsx')
ENR_weights_trimester1_sug_sig <- read_excel('ENR_weights_trimester1_GenR_suggesivepval.xlsx')
ENR_weights_trimester3_sug_sig <- read_excel('ENR_weights_trimester3_GenR_suggesivepval.xlsx')

overlapping_cpgs1 <- ENR_weights_totalinf_nom_sig[ENR_weights_totalinf_nom_sig$probeid %in% ENR_weights_totalinf_sug_sig$probeid,]
overlapping_cpgs2 <- ENR_weights_trimester1_nom_sig[ENR_weights_trimester1_nom_sig$probeid %in% ENR_weights_trimester1_sug_sig$probeid,]
overlapping_cpgs3 <- ENR_weights_trimester3_nom_sig[ENR_weights_trimester3_nom_sig$probeid %in% ENR_weights_trimester3_sug_sig$probeid,]

## Step a) combine MPS file with phenotype file 
# Load phenotype and MPS data files

datTest_totalinf <- readRDS('test_totalinf.rds')
datTest_tri1 <- readRDS('test_tri1.rds')
datTest_tri2 <- readRDS('test_tri2.rds')
datTest_tri3 <- readRDS('test_tri3.rds')

mps <- read_excel('MPS_total_infection_GenR_testset_suggestivepval.xlsx')
mps_tri1 <- read_excel('MPS_trimester1_infection_GenR_testset_suggestivepval.xlsx')
mps_tri2 <- read_excel('MPS_trimester2_infection_GenR_testset_suggestivepval.xlsx')
mps_tri3 <- read_excel('MPS_trimester3_infection_GenR_testset_suggestivepval.xlsx')

# create descriptive plots on mps
library(plotly)

plot_ly(data = mps, x = ~MPS_total_infection, type = "histogram", marker = list(color = "purple")) %>%
  layout(title = "A. Histogram of MPS total infection",
         xaxis = list(title = "MPS total infection"),
         yaxis = list(title = "Frequency"))

plot_ly(data = mps_tri1, x = ~MPS_trimester1_infection, type = "histogram", marker = list(color = "purple")) %>%
  layout(title = "B. Histogram of MPS trimester 1",
         xaxis = list(title = "MPS trimester 1"),
         yaxis = list(title = "Frequency"))

plot_ly(data = mps_tri2, x = ~MPS_trimester2_infection, type = "histogram", marker = list(color = "purple")) %>%
  layout(title = "C. Histogram of MPS trimester 2",
         xaxis = list(title = "MPS trimester 2"),
         yaxis = list(title = "Frequency"))

plot_ly(data = mps_tri3, x = ~MPS_trimester3_infection, type = "histogram", marker = list(color = "purple")) %>%
  layout(title = "D. Histogram of MPS trimester 3",
         xaxis = list(title = "MPS trimester 3"),
         yaxis = list(title = "Frequency"))

# Subset the phenotype data; exclude CpGs
phenodat.totalinf <- dplyr::select(datTest_totalinf, c('IDC', 'IDM', 'MOTHER', 'sumscore_inf_tri1', 'sumscore_inf_tri2', 'sumscore_inf_tri3', 'sumscore_inf_tot', 'sumscore_inf_tot_standardized', 'sumscore_inf_tri2_standardized', 'sumscore_inf_tri1_standardized', 'sumscore_inf_tri3_standardized','GENDER', 'SMOKE_ALL', 'WEIGHT', 'GESTBIR', 'PARITY', 'EDUCM_3groups', 'ethnicity', 'Sample_Plate', 'Bcell', 'CD4T', 'CD8T', 'Gran', 'Mono', 'NK', 'nRBC', 'assay','AGE_M_v2')) 

phenodat.trimester1 <- dplyr::select(datTest_tri1, c('IDC', 'IDM', 'MOTHER', 'sumscore_inf_tri1', 'sumscore_inf_tri2', 'sumscore_inf_tri3', 'sumscore_inf_tot', 'sumscore_inf_tot_standardized', 'sumscore_inf_tri2_standardized', 'sumscore_inf_tri1_standardized', 'sumscore_inf_tri3_standardized','GENDER', 'SMOKE_ALL', 'WEIGHT', 'GESTBIR', 'PARITY', 'EDUCM_3groups', 'ethnicity', 'Sample_Plate', 'Bcell', 'CD4T', 'CD8T', 'Gran', 'Mono', 'NK', 'nRBC', 'assay','AGE_M_v2')) 

phenodat.trimester2 <- dplyr::select(datTest_tri2, c('IDC', 'IDM', 'MOTHER', 'sumscore_inf_tri1', 'sumscore_inf_tri2', 'sumscore_inf_tri3', 'sumscore_inf_tot', 'sumscore_inf_tot_standardized', 'sumscore_inf_tri2_standardized', 'sumscore_inf_tri1_standardized', 'sumscore_inf_tri3_standardized','GENDER', 'SMOKE_ALL', 'WEIGHT', 'GESTBIR', 'PARITY', 'EDUCM_3groups', 'ethnicity', 'Sample_Plate', 'Bcell', 'CD4T', 'CD8T', 'Gran', 'Mono', 'NK', 'nRBC', 'assay','AGE_M_v2')) 

phenodat.trimester3 <- dplyr::select(datTest_tri3, c('IDC', 'IDM', 'MOTHER', 'sumscore_inf_tri1', 'sumscore_inf_tri2', 'sumscore_inf_tri3', 'sumscore_inf_tot', 'sumscore_inf_tot_standardized', 'sumscore_inf_tri2_standardized', 'sumscore_inf_tri1_standardized', 'sumscore_inf_tri3_standardized','GENDER', 'SMOKE_ALL', 'WEIGHT', 'GESTBIR', 'PARITY', 'EDUCM_3groups', 'ethnicity', 'Sample_Plate', 'Bcell', 'CD4T', 'CD8T', 'Gran', 'Mono', 'NK', 'nRBC', 'assay','AGE_M_v2')) 

# Merge testing phenotype data and MPS file
df_totalinfections <- cbind(mps, phenodat.totalinf)
df_trimester1 <- cbind(mps_tri1, phenodat.trimester1)
df_trimester2 <- cbind(mps_tri2, phenodat.trimester2)
df_trimester3 <- cbind(mps_tri3, phenodat.trimester3)

## Step b) associate prenatal infections with MPS of infections 
# Create plots of linear regresion for total infections 
scatter <- ggplot(df_totalinfections, aes(x = sumscore_inf_tot_standardized, y = scale(MPS_total_infection))) +
  geom_point(color = "#bc5090") +
  labs(title = "Regression Plot", x = "Prenatal infection sum score", y = "Methylation profile score of infections")

reg_line <- geom_smooth(data = df_totalinfections, method = "lm", se = FALSE, color = "#58508d")

plot <- ggplotly(scatter + reg_line)
plot <- plot %>% layout(plot_bgcolor = "lightblue")
plot 

# b1. Linear regression 
lm_totalinf <- summary(lm(scale(MPS_total_infection) ~  sumscore_inf_tot_standardized+ GENDER + SMOKE_ALL + AGE_M_v2 + PARITY + EDUCM_3groups + GESTBIR + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC, data =df_totalinfections))$coefficients[2,] 

lm_tri1 <- summary(lm(scale(MPS_trimester1_infection) ~  sumscore_inf_tri1_standardized+ GENDER + SMOKE_ALL + AGE_M_v2 + PARITY + EDUCM_3groups + GESTBIR + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC, data =df_trimester1))$coefficients[2,] 

lm_tri2 <- summary(lm(scale(MPS_trimester2_infection) ~  sumscore_inf_tri2_standardized+ GENDER + SMOKE_ALL + AGE_M_v2 + PARITY + EDUCM_3groups + GESTBIR + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC, data =df_trimester2))$coefficients[2,] 

lm_tri3 <- summary(lm(scale(MPS_trimester3_infection) ~  sumscore_inf_tri3_standardized+ GENDER + SMOKE_ALL + AGE_M_v2 + PARITY + EDUCM_3groups + GESTBIR + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC, data =df_trimester3))$coefficients[2,] 

library(tibble) # create df with results to later save to excel 
lm_total <- rbind(lm_totalinf, lm_tri1, lm_tri2, lm_tri3)
lm_total <- as.data.frame(lm_total)
lm_total <-  tibble::rownames_to_column(lm_total, "Exposure")

# b2. Calculate incremental r2 --> this indicates prediction peformance of MPS 
# mod 1=basic model, including all covariates in your EWAS
mod1 <- lm(sumscore_inf_tot_standardized ~  GENDER + SMOKE_ALL + AGE_M_v2 + PARITY + EDUCM_3groups + GESTBIR + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC, df_totalinfections)
r.mod1 = summary(mod1)$r.squared   

mod2 <- lm(sumscore_inf_tot_standardized ~ scale(MPS_total_infection) +  GENDER + SMOKE_ALL + AGE_M_v2 + PARITY + EDUCM_3groups + GESTBIR + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC, df_totalinfections)
r.mod2 = summary(mod2)$r.squared
r2change =  r.mod2 - r.mod1 #0.04852256

mod3 <- lm(sumscore_inf_tri1_standardized ~  GENDER + SMOKE_ALL + AGE_M_v2 + PARITY + EDUCM_3groups + GESTBIR + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC, df_trimester1)
r.mod3 = summary(mod3)$r.squared   

mod4 <- lm(sumscore_inf_tri1_standardized ~ scale(MPS_trimester1_infection) +  GENDER + SMOKE_ALL + AGE_M_v2 + PARITY + EDUCM_3groups + GESTBIR + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC, df_trimester1)
r.mod4 = summary(mod4)$r.squared
r2change2 =  r.mod4 - r.mod3 #0.02322675

mod5 <- lm(sumscore_inf_tri2_standardized ~  GENDER + SMOKE_ALL + AGE_M_v2 + PARITY + EDUCM_3groups + GESTBIR + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC, df_trimester2)
r.mod5 = summary(mod5)$r.squared   

mod6 <- lm(sumscore_inf_tri2_standardized ~ scale(MPS_trimester2_infection) +  GENDER + SMOKE_ALL + AGE_M_v2 + PARITY + EDUCM_3groups + GESTBIR + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC, df_trimester2)
r.mod6 = summary(mod6)$r.squared
r2change3 =  r.mod6 - r.mod5 #0.01999763

mod7 <- lm(sumscore_inf_tri3_standardized ~  GENDER + SMOKE_ALL + AGE_M_v2 + PARITY + EDUCM_3groups + GESTBIR + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC, df_trimester3)
r.mod7 = summary(mod7)$r.squared   

mod8 <- lm(sumscore_inf_tri3_standardized ~ scale(MPS_trimester3_infection) +  GENDER + SMOKE_ALL + AGE_M_v2 + PARITY + EDUCM_3groups + GESTBIR + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC, df_trimester3)
r.mod8 = summary(mod8)$r.squared
r2change4 =  r.mod8 - r.mod7 #0.01785841

lm_total$r2_change <- c(r2change, r2change2, r2change3, r2change4) # add results to lm validation df 

# b3. ROC curves & calculating AUC to determine discriminative ability 
library(pROC)

resid<- lm(scale(MPS_total_infection)~GENDER + SMOKE_ALL + AGE_M_v2 + PARITY + EDUCM_3groups + GESTBIR + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC, data =df_totalinfections) 
resid2<- lm(scale(MPS_trimester1_infection) ~ GENDER + SMOKE_ALL + AGE_M_v2 + PARITY + EDUCM_3groups + GESTBIR + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC, data =df_trimester1) 
resid3<- lm(scale(MPS_trimester2_infection) ~ GENDER + SMOKE_ALL + AGE_M_v2 + PARITY + EDUCM_3groups + GESTBIR + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC, data =df_trimester2) 
resid4<- lm(scale(MPS_trimester3_infection) ~ GENDER + SMOKE_ALL + AGE_M_v2 + PARITY + EDUCM_3groups + GESTBIR + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC, data =df_trimester3) 

df_totalinfections$MPS_total_infection_resid <- residuals(resid)
df_trimester1$MPS_trimester1_infection_resid <- residuals(resid2)
df_trimester2$MPS_trimester2_infection_resid <- residuals(resid3)
df_trimester3$MPS_trimester3_infection_resid <- residuals(resid4)

mean(df_totalinfections$sumscore_inf_tot_standardized) + sd(df_totalinfections$sumscore_inf_tot_standardized) # 0.9995902
mean(df_trimester1$sumscore_inf_tri1_standardized) + sd(df_trimester1$sumscore_inf_tri1_standardized) # 1.061525
mean(df_trimester2$sumscore_inf_tri2_standardized) + sd(df_trimester2$sumscore_inf_tri2_standardized) # 0.9886368 
mean(df_trimester3$sumscore_inf_tri3_standardized) + sd(df_trimester3$sumscore_inf_tri3_standardized) # 1.015416

df_totalinfections$infections_categorical <- as.factor(ifelse(df_totalinfections$sumscore_inf_tot_standardized > 0.999, 1, 0)) #Categorize infections (based on +1 SD of mean)
df_trimester1$infections_categorical <- as.factor(ifelse(df_trimester1$sumscore_inf_tri1_standardized > 1.016, 1, 0)) 
df_trimester2$infections_categorical <- as.factor(ifelse(df_trimester2$sumscore_inf_tri2_standardized > 0.9886368, 1, 0)) 
df_trimester3$infections_categorical <- as.factor(ifelse(df_trimester3$sumscore_inf_tri3_standardized > 1.015416, 1, 0)) 

roc1 <- roc(infections_categorical ~ MPS_total_infection_resid, data = df_totalinfections) #auc = 0.63
roc2 <- roc(infections_categorical ~ MPS_trimester1_infection_resid, data = df_trimester1) #auc = 0.64
roc3 <- roc(infections_categorical ~ MPS_trimester2_infection_resid, data = df_trimester2) #auc = 0.61
roc4 <- roc(infections_categorical ~ MPS_trimester3_infection_resid, data = df_trimester3) #auc = 0.47

library(writexl)
lm_total$auc <- c(roc1$auc, roc2$auc, roc3$auc, roc4$auc)
write_xlsx(lm_total, 'lm_infection_and_mps.xlsx') # save results of lm and r2 to excel 

# create a plot ROC for total infections 
resid <- glm(infections_categorical ~ MPS_total_infection_resid, data = df_totalinfections, family = "binomial")

roc1 <- roc(df_totalinfections$infections_categorical, resid$fitted.values)

auc_value <- auc(roc1)

roc_plot <- plot_ly(data = data.frame(specificity = roc1$specificities, sensitivity = roc1$sensitivities),
                    x = ~specificity, y = ~sensitivity, type = 'scatter', mode = 'lines',
                    line = list(color = '#bc5090', width = 2), name = 'ROC Curve')

roc_plot <- roc_plot %>%
  add_trace(type = 'scatter', x = c(0.7), y = c(0.2), mode = 'text',
            text = paste('AUC =', round(auc_value, 3)), showlegend = FALSE)

roc_plot <- layout(roc_plot, title = "ROC Curve",
                   xaxis = list(title = "1 - Specificity"),
                   yaxis = list(title = "Sensitivity"),
                   plot_bgcolor = 'lightblue')

print(roc_plot)

# b4. Associating MPS of infections to child outcomes
setwd('path_to_data')

# First load and combine child pheno data
library(foreign)
library(dplyr)
second_hits_df <- read.spss('Second_hits_DF.sav', to.data.frame = T) 
second_hits_df <- dplyr::select(second_hits_df, c('IDC', 'IDM', 'sum_int_14', 'sum_ext_14', 'cbcl_sum_14'))

child_obesity_f13 <- read.spss('CHILDGROWTH13_10122020.sav', to.data.frame = T)
child_obesity_f13 <- dplyr::select(child_obesity_f13, c('IDC', 'sdsbmiforage13childT'))
child_obesity_f13$bmi_f13 <- as.numeric(child_obesity_f13$sdsbmiforage13childT)

child_resp2 <- read.spss('GR1093_D1-D10_15052020.sav', to.data.frame = T)
child_resp_f13 <- dplyr::select(child_resp2, c('IDC', 'D001001R9301_c'))
child_resp_f13 <- rename(child_resp_f13, 'asthma_f13' = D001001R9301_c)
child_resp_f13$asthma_f13_cat <- ifelse(child_resp_f13$asthma_f13 == 'Yes', 1,0)
child_resp_f13$asthma_f13 <- as.numeric(child_resp_f13$asthma_f13)

child.pheno <- merge(second_hits_df, child_resp_f13, by = 'IDC', all.x =  T)
child.pheno <- merge(child.pheno, child_obesity_f13, by = 'IDC', all.x = T)

# Combine child pheno with mps data 
child.mps.pheno <- merge(df_totalinfections, child.pheno, by = 'IDC', all.x = T)
child.mps_tri1.pheno <- merge(df_trimester1, child.pheno, by = 'IDC', all.x = T)
child.mps_tri2.pheno <- merge(df_trimester2, child.pheno, by = 'IDC', all.x = T)
child.mps_tri3.pheno <- merge(df_trimester3, child.pheno, by = 'IDC', all.x = T)

# Run regressions (categorical)
summary(glm(asthma_f13_cat ~ scale(MPS_total_infection) +  GENDER + SMOKE_ALL + AGE_M_v2 + PARITY + EDUCM_3groups + GESTBIR + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC, data = child.mps.pheno,family = 'binomial'), conf.int = T)
summary(glm(asthma_f13_cat ~ scale(MPS_trimester1_infection) +  GENDER + SMOKE_ALL + AGE_M_v2 + PARITY + EDUCM_3groups + GESTBIR + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC, data = child.mps_tri1.pheno,family = 'binomial'), conf.int = T)
summary(glm(asthma_f13_cat ~ scale(MPS_trimester2_infection) +  GENDER + SMOKE_ALL + AGE_M_v2 + PARITY + EDUCM_3groups + GESTBIR + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC, data = child.mps_tri2.pheno,family = 'binomial'), conf.int = T)
summary(glm(asthma_f13_cat ~ scale(MPS_trimester3_infection) +  GENDER + SMOKE_ALL + AGE_M_v2 + PARITY + EDUCM_3groups + GESTBIR + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC, data = child.mps_tri3.pheno,family = 'binomial'), conf.int = T)

# Run regressions (continuous)
child_pheno <- c('sum_int_14', 'sum_ext_14', 'cbcl_sum_14', 'bmi_f13') # create vector with all outcomes to loop over 

results_df <- data.frame() # Create an empty data frame to store results

for(x in child_pheno) {
  # lm formula
  a <- paste0(x, "~ scale(MPS_total_infection) + GENDER + SMOKE_ALL + AGE_M_v2 + PARITY + EDUCM_3groups + GESTBIR + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC")
  
  b <- paste0(x, "~ scale(MPS_trimester1_infection) + GENDER + SMOKE_ALL + AGE_M_v2 + PARITY + EDUCM_3groups + GESTBIR + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC")
  
  c <- paste0(x, "~ scale(MPS_trimester2_infection) + GENDER + SMOKE_ALL + AGE_M_v2 + PARITY + EDUCM_3groups + GESTBIR + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC")
  
  d <- paste0(x, "~ scale(MPS_trimester3_infection) + GENDER + SMOKE_ALL + AGE_M_v2 + PARITY + EDUCM_3groups + GESTBIR + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC")
  
  #extract coefficients 
  coef_totalinf <- summary(lm(as.formula(a), data = child.mps.pheno))$coefficients[2,]
  coef_tri1inf <- summary(lm(as.formula(b), data = child.mps_tri1.pheno))$coefficients[2,]
  coef_tri2inf <- summary(lm(as.formula(c), data = child.mps_tri2.pheno))$coefficients[2,]
  coef_tri3inf <- summary(lm(as.formula(d), data = child.mps_tri3.pheno))$coefficients[2,]
  
  #Create a data frame with the results
  result_row <- data.frame(
    Phenotype = x,
    Total_Infection_Coefficient = coef_totalinf,
    Tri1_Infection_Coefficient = coef_tri1inf,
    Tri2_Infection_Coefficient = coef_tri2inf,
    Tri3_Infection_Coefficient = coef_tri3inf
  )
  
  # Append the results to the main data frame
  results_df <- rbind(results_df, result_row)
}

# Write the results to Excel
setwd('path_to_results')

results_df2 <- t(results_df)
results_df2 <- as.data.frame(results_df2)

write_xlsx(results_df2, path = "child_pheno_valid_mps_results_updated.xlsx")

## Step c) create correlation plot in genR testset 
total_df <- merge(child.mps.pheno, mps_tri1, by = 'Sample_ID', all.x = T)
total_df2 <-  merge(total_df, mps_tri2, by = 'Sample_ID', all.x = T)
total_df3 <-  merge(total_df2, mps_tri3, by = 'Sample_ID', all.x = T)

corr_vars <- c("GENDER", "GESTBIR", "EDUCM_3groups", "AGE_M_v2", "cbcl_sum_14", "sum_int_14", "sum_ext_14", "bmi_f13",  "SMOKE_ALL", "sumscore_inf_tri1", "sumscore_inf_tri2", "sumscore_inf_tri3", "MPS_trimester1_infection", "MPS_trimester2_infection", "MPS_trimester3_infection", "MPS_total_infection", "sumscore_inf_tot", 'asthma_f13') 

selected_df <- total_df3[corr_vars]
selected_df <- rename(selected_df, "13. Child sex" = GENDER, "9. Gestational age at birth" = GESTBIR, "11. Maternal education" = EDUCM_3groups, "10. Maternal age" = AGE_M_v2, "14. CBCL total behavioral problems (14)" = cbcl_sum_14, "15. CBCL internalizing problems (14)" = sum_int_14, "16. CBCL externalizing problems (14)" = sum_ext_14, "17. BMI (age 14)" = bmi_f13, "12. Maternal smoking" = SMOKE_ALL, "2. Infection sum score (trimester 1)" = sumscore_inf_tri1, "3. Infection sum score (trimester 2)" = sumscore_inf_tri2, "4. Infection sum score (trimester 3)" = sumscore_inf_tri3, "1. Total infection sum score" = sumscore_inf_tot, "5. MPS (total infections)" = MPS_total_infection, "6. MPS (trimester 1)" = MPS_trimester1_infection, "7. MPS (trimester 2)" = MPS_trimester2_infection, "8. MPS (trimester 3)" = MPS_trimester3_infection, '18. Asthma (age 14)' = asthma_f13) 

selected_df$`13. Child sex` <- as.numeric(selected_df$`13. Child sex`)
selected_df$`11. Maternal education` <- as.numeric(selected_df$`11. Maternal education`)
selected_df$`17. BMI (age 14)` <- as.numeric(selected_df$`17. BMI (age 14)`)
selected_df$`12. Maternal smoking` <- as.numeric(selected_df$`12. Maternal smoking`)

selected_df <- subset(selected_df, complete.cases(selected_df))

library(corrplot)

# Specify the desired order of variables
desired_order <- c(
  "1. Total infection sum score",
  "2. Infection sum score (trimester 1)",
  "3. Infection sum score (trimester 2)",
  "4. Infection sum score (trimester 3)",
  "5. MPS (total infections)",
  "6. MPS (trimester 1)",
  "7. MPS (trimester 2)",
  "8. MPS (trimester 3)",
  "9. Gestational age at birth",
  "10. Maternal age",
  "11. Maternal education",
  "12. Maternal smoking",
  "13. Child sex",
  "14. CBCL total behavioral problems (14)",
  "15. CBCL internalizing problems (14)",
  "16. CBCL externalizing problems (14)",
  "17. BMI (age 14)",
  "18. Asthma (age 14)"
)

# Reorder the correlation matrix based on the desired order
correlation_ordered <- correlation[desired_order, desired_order]

# Plot the correlation matrix with the specified order
corrplot::corrplot(
  correlation_ordered, 
  method = 'color',
  addCoef.col = "black", 
  number.cex = 0.6,
  type = 'lower', 
  diag = FALSE,  
  tl.col = 'black', 
  tl.cex = 0.7, 
  sig.level = 0.05,
  col = colorRampPalette(c("midnightblue", "white", "darkred"))(100),
  colnames = colnames(correlation_ordered)
)

cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

p.mat <- cor.mtest(correlation)

### Post hoc analysis: create alspac like infection score in genR and associate that to MPS in GenR
# create scores 
setwd('path_to_data')
second_hits_df <- read.spss('Second_hits_DF.sav', to.data.frame = T) 
second_hits_df$alspac_infections_tot <- second_hits_df$flu_tri1 + second_hits_df$flu_tri2 + second_hits_df$flu_tri3 + second_hits_df$UWI_tri1 + second_hits_df$UWI_tri2 + second_hits_df$UWI_tri3
second_hits_df$alspac_infections_tri1 <- second_hits_df$flu_tri1 + second_hits_df$UWI_tri1
second_hits_df$alspac_infections_tri2 <-  second_hits_df$flu_tri2 + second_hits_df$UWI_tri2
second_hits_df$alspac_infections_tri3 <-  second_hits_df$flu_tri3 + second_hits_df$UWI_tri3
second_hits_df <- dplyr::select(second_hits_df, c('IDC', 'alspac_infections_tot', "alspac_infections_tri1", "alspac_infections_tri2", "alspac_infections_tri3"))

df_totalinfections2 <- merge(df_totalinfections,second_hits_df, by = 'IDC', all.x=T)
df_totalinfections3 <- merge(df_totalinfections2, mps_tri1, by = 'Sample_ID', all.x = T)
df_totalinfections4 <- merge(df_totalinfections3, mps_tri2, by = 'Sample_ID', all.x = T)
df_totalinfections5 <- merge(df_totalinfections4, mps_tri3, by = 'Sample_ID', all.x = T)

# lm with infections 
a_lm_totalinf <- summary(lm(scale(MPS_total_infection) ~  scale(alspac_infections_tot)+ GENDER + SMOKE_ALL + AGE_M_v2 + PARITY + EDUCM_3groups + GESTBIR + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC, data =df_totalinfections5))$coefficients[2,] 

a_lm_tri1 <- summary(lm(scale(MPS_trimester1_infection) ~  scale(alspac_infections_tri1)+ GENDER + SMOKE_ALL + AGE_M_v2 + PARITY + EDUCM_3groups + GESTBIR + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC, data =df_totalinfections5))$coefficients[2,] 

a_lm_tri2 <- summary(lm(scale(MPS_trimester2_infection) ~  scale(alspac_infections_tri2)+ GENDER + SMOKE_ALL + AGE_M_v2 + PARITY + EDUCM_3groups + GESTBIR + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC, data =df_totalinfections5))$coefficients[2,] 

a_lm_tri3 <- summary(lm(scale(MPS_trimester3_infection) ~  scale(alspac_infections_tri3)+ GENDER + SMOKE_ALL + AGE_M_v2 + PARITY + EDUCM_3groups + GESTBIR + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC, data =df_totalinfections5))$coefficients[2,] 

library(tibble) # create df with results to later save to excel 
a_lm_total <- rbind(a_lm_totalinf, a_lm_tri1, a_lm_tri2, a_lm_tri3)
a_lm_total <- as.data.frame(a_lm_total)
a_lm_total <-  tibble::rownames_to_column(a_lm_total, "Exposure")

# calculate r2
mod1 <- lm(scale(alspac_infections_tot) ~  GENDER + SMOKE_ALL + AGE_M_v2 + PARITY + EDUCM_3groups + GESTBIR + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC, df_totalinfections5)
r.mod1 = summary(mod1)$r.squared   

mod2 <- lm(scale(alspac_infections_tot) ~ scale(MPS_total_infection) +  GENDER + SMOKE_ALL + AGE_M_v2 + PARITY + EDUCM_3groups + GESTBIR + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC, df_totalinfections5)
r.mod2 = summary(mod2)$r.squared
r2change =  r.mod2 - r.mod1 

mod3 <- lm(scale(alspac_infections_tri1) ~  GENDER + SMOKE_ALL + AGE_M_v2 + PARITY + EDUCM_3groups + GESTBIR + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC, df_totalinfections5)
r.mod3 = summary(mod3)$r.squared   

mod4 <- lm(scale(alspac_infections_tri1) ~ scale(MPS_trimester1_infection) +  GENDER + SMOKE_ALL + AGE_M_v2 + PARITY + EDUCM_3groups + GESTBIR + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC, df_totalinfections5)
r.mod4 = summary(mod4)$r.squared
r2change2 =  r.mod4 - r.mod3

mod5 <- lm(scale(alspac_infections_tri2) ~  GENDER + SMOKE_ALL + AGE_M_v2 + PARITY + EDUCM_3groups + GESTBIR + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC, df_totalinfections5)
r.mod5 = summary(mod5)$r.squared   

mod6 <- lm(scale(alspac_infections_tri2) ~ scale(MPS_trimester2_infection) +  GENDER + SMOKE_ALL + AGE_M_v2 + PARITY + EDUCM_3groups + GESTBIR + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC, df_totalinfections5)
r.mod6 = summary(mod6)$r.squared
r2change3 =  r.mod6 - r.mod5 

mod7 <- lm(scale(alspac_infections_tri3) ~  GENDER + SMOKE_ALL + AGE_M_v2 + PARITY + EDUCM_3groups + GESTBIR + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC, df_totalinfections5)
r.mod7 = summary(mod7)$r.squared   

mod8 <- lm(scale(alspac_infections_tri3) ~ scale(MPS_trimester3_infection) +  GENDER + SMOKE_ALL + AGE_M_v2 + PARITY + EDUCM_3groups + GESTBIR + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC, df_totalinfections5)
r.mod8 = summary(mod8)$r.squared
r2change4 =  r.mod8 - r.mod7 

a_lm_total$r2_change <- c(r2change, r2change2, r2change3, r2change4)

# save results 
setwd('path_to_results')

write_xlsx(a_lm_total, path = "lm_alspac-like-infection_and_mps.xlsx")

#\ END OF SCRIPT. 
