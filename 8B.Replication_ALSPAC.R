#########################################################################
#####################PROJECT: PRENATAL INFECTION EWAS####################
#########################################################################

#Authors: Anna Suleri & Rosa Mulder

#This script is to:
#'a) select variables of interest in dataframe, apply inclusion and exclusion criteria, apply multiple imputation
#'b) Create methylation profile score of infections based on weights from GenR-train set
#'c) Associate MPS-infection with asthma/bmi/behavior as child phenotype (obtain beta, SE and p-value)

#############---------Step 1: Prepare pheno dataframe---------#############
### Set wd & load libraries
set.seed(2023)

libraries <- c('foreign', 'tidyverse', 'devtools', 'ALSPAC.helpR', 'mice', 'readxl', 'tibble', 'writexl', 'data.table', 'haven', 'dplyr', 'MASS', 'Matrix', 'DescTools', "meffil", 'pROC', 'plotly')

install_github("SereDef/ALSPAC.helpR")

invisible(lapply(libraries, require, character.only = T)) 

setwd('set_path_to_data') 

### Read in data
phenotype.data <- read.spss('path_to_data', to.data.frame = T) 

## Select variables of interest 
# id variable: cidB4195
pheno <- dplyr::select(phenotype.data, c('cidB4195','kz021', 'e171', 'e173', 'e175', 'e176', 'b670', 'e178', 'bestgest', 'd994', 'qlet', 'c642', 'c630', 'c631', "c632" ,'c641', 'kb611', 'FJMR022a', 'fh6876', 'fh6878', 'fh6870', 'fh6874', 'fh6873','tb1070')) 

#check asthma
summary(pheno$tb1070)
pheno$asthma <- NA
pheno$asthma[pheno$tb1070=="Yes, asthma" | pheno$tb1070=="Yes, asthma & eczema"] <- "yes"
pheno$asthma[pheno$tb1070=="Other text answer" | pheno$tb1070=="Yes, eczema" | pheno$tb1070=="No"] <- "no"

pheno$asthma <- as.factor(pheno$asthma)

table(pheno$asthma, useNA="always")
table(as.numeric(pheno$asthma), useNA="always")

saveRDS(pheno, 'pheno_240320.rds')

## Re categorize data and check structure of all variables 
## Create maternal smoking variable
pheno$msmoke <- NA
pheno$msmoke[pheno$e171 == "Y" | pheno$e173 == "Y" |  pheno$e175 == "Y" | pheno$e176 == "Y"] <- "3_continued"             
pheno$msmoke[pheno$b670 != 0  & pheno$e178 == "Not at all"] <- "2_stopped"
pheno$msmoke[pheno$e171 == "N" & pheno$e173 == "N" & pheno$e175 == "N" & pheno$e176 == "N" & pheno$b670 == 0 & pheno$e178 == "Not at all"] <- "1_no smoking"
pheno$msmoke <- as.factor(pheno$msmoke)

# check
table(pheno$msmoke, useNA="always")

str(pheno$msmoke)

pheno[1:10,c("e171","e173","e175","e176","b670","e178","msmoke")]

# drop columns we don't need anymore 
pheno <- dplyr::select(pheno, -c('e171', 'e173', 'e175', 'e176', 'b670', 'e178')) 

## Create continuous maternal education variable
pheno$EDUCM <- NA
pheno$EDUCM[pheno$c642 == "Yes"]  <- 0 #No
pheno$EDUCM[pheno$c630 == "Yes"] <- 1 #CSE
pheno$EDUCM[pheno$c631 == "Yes" | pheno$c632 == "Yes"] <- 2 #O-level
pheno$EDUCM[pheno$c641 == "Yes"]  <- 3 #University

summary(pheno$EDUCM)
str(pheno$EDUCM)

# drop columns we don't need anymore 
pheno <- dplyr::select(pheno, -c('c642', 'c630', "c632" ,'c631', 'c641')) 

## Rename variables
pheno <- pheno %>% rename('GENDER' = kz021, 'm_age' = d994, 'twin' = kb611, 'bmi_f17' = FJMR022a, 'depression_f15' = fh6876, 'anxiety_f15' = fh6878, 'adhd_f15' = fh6870, 'cd_f15' = fh6874, 'odd_f15' = fh6873) 

## Check structure of vars and adapt if needed  
str(pheno)

pheno$bestgest <- as.numeric(as.character(pheno$bestgest))
pheno$m_age <- as.numeric(as.character(pheno$m_age))
pheno$bmi_f17 <- as.numeric(as.character(pheno$bmi_f17))

pheno$depression_f15 <- as.numeric(pheno$depression_f15)
pheno$anxiety_f15 <- as.numeric(pheno$anxiety_f15)
pheno$adhd_f15 <- as.numeric(pheno$adhd_f15)
pheno$cd_f15 <- as.numeric(pheno$cd_f15)
pheno$odd_f15 <- as.numeric(pheno$odd_f15)

#check
summary(pheno$bestgest)
summary(pheno$m_age)
summary(pheno$bmi_f17)
summary(pheno$depression_f15)
summary(pheno$anxiety_f15)
summary(pheno$adhd_f15)
summary(pheno$cd_f15)
summary(pheno$odd_f15)

#save
saveRDS(pheno, 'final_pheno_240320.rds')

### Merge link file 
link.data <- read.spss("path_to_data", to.data.frame=T) 
phenotype_link.data <- merge(pheno, link.data, by = c("cidB4195", "qlet"))

dim(phenotype_link.data)

### Merge methylation meta data
#samples <- read.csv("~/Alspac_Tempo/methylation/data/samplesheet/samplesheet.csv") 
load("/home/r.mulder/GENR3/Projects/AFAR_ALSPAC/dnam_450_g0m_g1/samplesheet/data.Robj")   # samplesheet
samples <- samplesheet
samplesheet$qlet <- samplesheet$QLET

dim(samplesheet)
table(samplesheet$duplicate.rm, useNA="always")
samples <- samplesheet[is.na(samplesheet$duplicate.rm),]
dim(samples)

# Filter for birth
table(samples$time_point, useNA="always")
samples_birth.data <- samples[samples$time_point == "cord", ] #913 participants (913 rows and 33 columns)
dim(samples_birth.data)
phenotype_link_sample.data <- merge(phenotype_link.data, samples_birth.data, by = c("dnam_450_g0m_g1", 'qlet'))  
dim(phenotype_link_sample.data)

# Merge cell count data
cell_counts_combined.data <- readRDS("/home/r.mulder/GENR3/Projects/ALSPAC_WBCs/dnam_450_g0m_g1_salas_combinedcordblood_240307.rds") #

dim(cell_counts_combined.data)

#phenotype_link_sample_counts.data <- merge(phenotype_link_sample.data, cell_counts_combined.data, by.x = "Sample_Name", by.y = "IID") 
phenotype_link_sample_counts.data <- merge(phenotype_link_sample.data, cell_counts_combined.data, by.x = "Sample_Name", by.y = "row.names")

dim(phenotype_link_sample_counts.data)

### Apply exclusion criteria
summary(phenotype_link_sample_counts.data$twin)

table(duplicated(phenotype_link_sample_counts.data$cidB4195), useNA="always")

# Remove twins
df1 <- subset(phenotype_link_sample_counts.data, twin == 'No') # -41 participants in tempo dataset (n=872) [-5 twins, -35 NA, 1 consent withdrawn by mother]

#check
summary(df1$twin)
dim(df1)
nrow(df1)-nrow(phenotype_link_sample_counts.data)

# Remove 1 sibling
table(duplicated(df1$cidB4195), useNA="always")

df1$id <- make_idc(mom.id = 'cidB4195', parity='qlet', data=df1)
df2    <- rm_siblings(method = 'missing', data = df1, idc='id') 

dim(df2)                                                   
nrow(df2)-nrow(df1)

# check final participants 
dim(df2) 

saveRDS(df2, 'df_after_exclusion_240320.rds')

#############---------Step 2: Prepare methylation dataframe---------#############
### Load and merge in methylation data 

## Load birth methylation data ALSPAC
#dnam_alspac_orig.data <- meffil.gds.methylation("~/Alspac_Tempo/methylation/data/betas/450.gds")
load("/home/r.mulder/GENR3/Projects/AFAR_ALSPAC/dnam_450_g0m_g1/betas/data.Robj")

dnam_alspac_orig.data <- betas

rm(betas)

dim(dnam_alspac_orig.data)

## Keep only birth samples
#samples <- read.csv("~/Alspac_Tempo/methylation/data/samplesheet/samplesheet.csv") 
load("/home/r.mulder/GENR3/Projects/AFAR_ALSPAC/dnam_450_g0m_g1/samplesheet/data.Robj")   # samplesheet
samples <- samplesheet
samples$qlet <- samples$QLET

dim(samples)

## Filter for birth
table(samples$time_point, useNA="always")

samples_birth.data <- samples[samples$time_point == "cord", ]

dim(samples_birth.data)

## Filter the 450K methylation data for birth
dnam_alspac_cord.data <- dnam_alspac_orig.data[ ,samples_birth.data$Sample_Name]

dim(dnam_alspac_cord.data) #rows = cpg's and columns = id 

# ## Transpose (now cpgs will be columns and ids will be rows)    
# dnam_alspac_cord <- as.data.frame(t(dnam_alspac_cord.data))
# dim(dnam_alspac_cord)
# 
# ## Save dataset
# saveRDS(dnam_alspac_cord, 'dnam_alspac_cord_240320.rds')

#\




#############---------Step 3: Compute MPSes---------#############
### Create MPS with weights from GenR-train set
## Step 1: Construct MPS based on ENR-weights from GenR. The base data = weights/beta from elastic net regression in GENR & target Data: Methylation at birth in ALSPAC. 

# 1)load GENR ENR-weight file 
#setwd("path_to_data_MPS")
S1 <- read_xlsx("ENR_weights_total_infection_GenR_suggesivepval.xlsx")  
S2 <- read_xlsx("ENR_weights_trimester1_GenR_suggesivepval.xlsx")    
S3 <- read_xlsx("ENR_weights_trimester2_GenR_suggesivepval.xlsx")    
S4 <- read_xlsx("ENR_weights_trimester3_GenR_suggesivepval.xlsx")    

dim(S1)
dim(S2)
dim(S3)
dim(S4)

S1 <- as.data.frame(S1)
S2 <- as.data.frame(S2)
S3 <- as.data.frame(S3)
S4 <- as.data.frame(S4)

# 2)Prepare the target data: the ALSPAC methlyation at birth 
# load  methylation data 
meth <- as.data.frame(dnam_alspac_cord.data)

# 3) match target data (DNAm in ALSPAC) with ENR-weight file, based on probeID.
# only include overlapped CpGs between target and base data
cpg_totalinf <- S1$probeid
cpg_tri1_inf <- S2$probeid
cpg_tri2_inf <- S3$probeid
cpg_tri3_inf <- S4$probeid

methsub1 <- meth[na.omit(match(cpg_totalinf,rownames(meth))),]
methsub2 <- meth[na.omit(match(cpg_tri1_inf,rownames(meth))),]
methsub3 <- meth[na.omit(match(cpg_tri2_inf,rownames(meth))),]
methsub4 <- meth[na.omit(match(cpg_tri3_inf,rownames(meth))),]

# Important! Make sure all rows keep the same (in the same order)
# needs to be TRUE
identical(rownames(methsub1), S1$probeid) 
identical(rownames(methsub2), S2$probeid)
identical(rownames(methsub3), S3$probeid)
identical(rownames(methsub4), S4$probeid)

# 4) MATRIX MULTIPLICATION (meth*weight)
#this multiples each CpG by the ENR weight
wmatrix1 = (methsub1 * S1$beta) 
dim(wmatrix1)
wmatrix2 = (methsub2 * S2$beta) 
dim(wmatrix2)
wmatrix3 = (methsub3 * S3$beta) 
dim(wmatrix3)
wmatrix4 = (methsub4 * S4$beta) 
dim(wmatrix4)

# 5) Create MPS _the sum scores (with all selected cpg) per person (summarizing weighted meth x each PP)
# MPS based on ENR-weight1
wmatrix1[1:5,1:5]
dim(wmatrix1)
colsum1 = colSums(wmatrix1, na.rm=TRUE)
mps_S1 = as.data.frame(colsum1)
head(mps_S1)
summary(mps_S1$colsum1)

# MPS based on ENR-weight2
wmatrix2[1:5,1:5]
dim(wmatrix2)
colsum2 = colSums(wmatrix2, na.rm=TRUE)
mps_S2 = as.data.frame(colsum2)
head(mps_S2)
summary(mps_S2$colsum2)

# MPS based on ENR-weight3
wmatrix3[1:5,1:5]
dim(wmatrix3)
colsum3 = colSums(wmatrix3, na.rm=TRUE)
mps_S3 = as.data.frame(colsum3)
head(mps_S3)
summary(mps_S3$colsum3)

# MPS based on ENR-weight4
wmatrix4[1:5,1:5]
dim(wmatrix4)
colsum4 = colSums(wmatrix4, na.rm=TRUE)
mps_S4 = as.data.frame(colsum4)
head(mps_S4)
summary(mps_S4$colsum4)

# to convert row names into first column
mps_S1 <- tibble::rownames_to_column(mps_S1, "Sample_ID")
colnames(mps_S1)[2] <- "MPS_total_infections"
head(mps_S1)

mps_S2 <- tibble::rownames_to_column(mps_S2, "Sample_ID")
colnames(mps_S2)[2] <- "MPS_tri1_infections"
head(mps_S2)

mps_S3 <- tibble::rownames_to_column(mps_S3, "Sample_ID")
colnames(mps_S3)[2] <- "MPS_tri2_infections"
head(mps_S3)

mps_S4 <- tibble::rownames_to_column(mps_S4, "Sample_ID")
colnames(mps_S4)[2] <- "MPS_tri3_infections"
head(mps_S4)

# 6) Save all MPS into one file
setwd('path_to_results')

MPS.infect <- Reduce(merge, list(mps_S1, mps_S2, mps_S3, mps_S4))

saveRDS(MPS.infect, "MPS_infections_ALSPAC_240320.rds")

#\

### Merge pheno file with mps file
mps_pheno_df <- merge(df2, MPS.infect, by.x = 'Sample_Name', by.y = 'Sample_ID', all.x = T)
dim(mps_pheno_df)

#remove unnecessary variables
mps_pheno_df <- dplyr::select(mps_pheno_df, -c("dnam_450_g0m_g1", "twin", "tb1070", "gi_1000g_g0m_g1", "ge_ht12_g1", "QLET",
                                               "Slide", "sentrix_row", "sentrix_col", "time_code", "time_point", "Sex", "BCD_plate", "sample_type", "additive",
                                               "age", "duplicate.rm", "genotypeQCkids", "genotypeQCmums", "id")) 

dim(mps_pheno_df)

##check structure
str(mps_pheno_df)

saveRDS(mps_pheno_df, "mps_pheno_df_240320.rds")

#\

#############---------Step 4: Descriptives non-imputed set---------#############

baselinevars <- c('asthma','depression_f15', 'anxiety_f15', 'adhd_f15', 'cd_f15', 'odd_f15',  'bmi_f17')

#function for continuous variables
summary_continuous <- function(x){
  standev <- sd(x, na.rm = T)
  meanvar <- mean(x, na.rm = T)
  print(paste0(round(meanvar, 1), '(', round(standev, 1), ')'))
}

#function for categorical variables
summary_categorical <- function(x){
  tab1 <- prop.table(table(x, useNA = 'always'))
  tab2 <- table(x, useNA = "always")
  print(paste(round(tab1 * 100, 1), '%', names(tab1), collapse = ','))
  print(paste(tab2, names(tab2)))
}


for(i in baselinevars){
  
  #x = i vars that are columns in the dataframe df
  x <- mps_pheno_df[, i]  #replace df_descriptives with your df name
  
  #show column name as heading per output
  message(i)
  
  #if else to apply correct function for vars type
  if (class(x) == 'numeric') {
    summary_continuous(x)
  }else{
    summary_categorical(x)
  }
}

#\

#############---------Step 5: Multiple imputation---------#############

## First check missingness 
miss_values_function <- function(mps_pheno_df){
  missvalues <- cbind("# NA" = sort(colSums(is.na(mps_pheno_df))),
                      "% NA" = round(sort(colMeans(is.na(mps_pheno_df))) * 100, 2))
  print(missvalues)
}

miss_values_function(mps_pheno_df)

## Impute missing variables 
imputation_function <- function(data, exclude_imp_vars, exclude_predictors, method = "default") {
  
  # Running setup imputation run
  if (method == "rf") { 
    imp0 <- mice(data, maxit = 0, method = "rf")
  } else {
    imp0 <- mice(data, maxit = 0, defaultMethod = c("norm", "logreg", "polyreg", "polr"))
  }
  
  # Imputation method matrix
  meth <- imp0$method
  meth[exclude_imp_vars] <- ""
  
  # Predictor matrix
  pred <- imp0$predictorMatrix
  pred[exclude_predictors] <- 0
  
  # Visit sequence
  visSeq <- imp0$visitSequence
  
  # Performing the imputation
  imp.test <- mice(data, method = meth, predictorMatrix = pred, visitSequence = visSeq, maxit = 50, m = 30, printFlag = TRUE, seed = 2023)  
  
  #get the dataframe name 
  dataName <- deparse(substitute(data))
  
  #assigning a dynamic name to the imp.test object
  imputedDataName <- paste0("imputedData_", dataName)
  assign(imputedDataName, imp.test)
  
  # Saving the imputed dataset as .RDS file
  saveRDS(imp.test, file = paste0(imputedDataName, "_240320", ".rds"))
  
  # Output 
  return(list(imp0 = imp0, imp.test = imp.test))
}

exclude_vars     <- c("Sample_Name", "qlet", "cidB4195")
exclude_vars_nrs <- which(colnames(mps_pheno_df) %in% exclude_vars)
exclude_vars_nrs

# impute
imputation_function(data = mps_pheno_df, exclude_imp_vars = exclude_vars_nrs, exclude_predictors = exclude_vars_nrs, method = "pmm")  

#load it
imp_df <- readRDS("imputedData_mps_pheno_df_240320.rds")

# Select last imputed df (on this we do our analyses) 
mps_pheno_imp30_df <- complete(imp_df, 30)  

#check 
any(is.na(mps_pheno_imp30_df))

#check structure
str(mps_pheno_imp30_df)

saveRDS(mps_pheno_imp30_df, "mps_pheno_imp30_df_240320.rds")

#\


#############---------Step 6: Descriptives imputed set---------#############

baselinevars <- c('asthma','depression_f15', 'anxiety_f15', 'adhd_f15', 'cd_f15', 'odd_f15',  'bmi_f17')

#function for continuous variables
summary_continuous <- function(x){
  standev <- sd(x, na.rm = T)
  meanvar <- mean(x, na.rm = T)
  print(paste0(round(meanvar, 1), '(', round(standev, 1), ')'))
}

#function for categorical variables
summary_categorical <- function(x){
  tab1 <- prop.table(table(x, useNA = 'always'))
  tab2 <- table(x, useNA = "always")
  print(paste(round(tab1 * 100, 1), '%', names(tab1), collapse = ','))
  print(paste(tab2, names(tab2)))
}


for(i in baselinevars){
  
  #x = i vars that are columns in the dataframe df
  x <- mps_pheno_imp30_df[, i]  #replace df_descriptives with your df name
  
  #show column name as heading per output
  message(i)
  
  #if else to apply correct function for vars type
  if (class(x) == 'numeric') {
    summary_continuous(x)
  }else{
    summary_categorical(x)
  }
}

#\


#############---------Step 7: Compute SVAs---------#############
covariates <- mps_pheno_imp30_df[,c("m_age","GENDER","msmoke","EDUCM","Bcell","CD4T","CD8T","Gran","Mono","NK","nRBC","MPS_tri1_infections","MPS_tri2_infections","MPS_tri3_infections")]
variable   <- mps_pheno_imp30_df[,"MPS_total_infections"]

#check
any(is.na(covariates))
any(is.na(variable))
identical(mps_pheno_imp30_df$Sample_Name, colnames(meth))        

#adjust meth
meth2 <- meth[,match(mps_pheno_imp30_df$Sample_Name, colnames(meth))]
identical(mps_pheno_imp30_df$Sample_Name, colnames(meth2))              

#set no. of SVs, most variable probes
n.sv = 20
most.variable = 10000 

#select autosomal probes
beta.sva <- as.matrix(meth2)
autosomal.sites <- meffil.get.autosomal.sites("450k")                    
length(autosomal.sites) 

autosomal.sites <- intersect(autosomal.sites, rownames(beta.sva))
length(autosomal.sites) 

beta.sva <- beta.sva[autosomal.sites,]
dim(beta.sva)

#check
any(is.na(beta.sva))

#select most variable probes
var.idx  <- order(rowVars(beta.sva, na.rm=T), decreasing=T)[1:most.variable] 
beta.sva <- beta.sva[var.idx,,drop=F]

dim(beta.sva)

#define model
cov.frame <- model.frame(~., data.frame(covariates, stringsAsFactors=F), na.action=na.pass)
mod0 <- model.matrix(~., cov.frame)

mod <- cbind(mod0, variable)

#sva
set.seed(2023)
sva.ret <- sva(beta.sva, mod=mod, mod0=mod0, n.sv=n.sv)

mps_pheno_imp30_svs_df  <- data.frame(mps_pheno_imp30_df, sva.ret$sv, stringsAsFactors=F)  

dim(mps_pheno_imp30_svs_df)

str(mps_pheno_imp30_svs_df)

#save
saveRDS(mps_pheno_imp30_svs_df, file = "mps_pheno_imp30_svs_df_240320.rds")

#\

#############---------Step 8: Run replication analyses---------#############

### Associate MPS-infection with child phenotypes
#setwd('path_to_where_you_want_to_store_results')

# First: logistic regression for asthma 
child_pheno <- c('asthma') #@Rosa replace with name of your asthma variable

asthma_results_df_Rosa <- data.frame() 

for(x in child_pheno) {
  # lm formula
  a <- paste0(x, " ~ scale(MPS_total_infections) + m_age + GENDER + msmoke + EDUCM + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + 
              X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20")      
  
  b <- paste0(x, " ~ scale(MPS_tri1_infections) + m_age + GENDER + msmoke + EDUCM + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + 
              X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20")
  
  c <- paste0(x, " ~ scale(MPS_tri2_infections) + m_age + GENDER + msmoke + EDUCM + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + 
              X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20")
  
  d <- paste0(x, "~ scale(MPS_tri3_infections) + m_age + GENDER + msmoke + EDUCM + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + 
              X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20")
  
  # Extract coefficients 
  coef_totalinf <- summary(glm(as.formula(a), data = mps_pheno_imp30_svs_df, family = 'binomial'))$coefficients[2,]
  coef_tri1inf <- summary(glm(as.formula(b), data = mps_pheno_imp30_svs_df, family = 'binomial'))$coefficients[2,]
  coef_tri2inf <- summary(glm(as.formula(c), data = mps_pheno_imp30_svs_df, family = 'binomial'))$coefficients[2,]
  coef_tri3inf <- summary(glm(as.formula(d), data = mps_pheno_imp30_svs_df, family = 'binomial'))$coefficients[2,]
  
  #Create a data frame with the results
  result_row <- data.frame(
    Phenotype = x,
    Total_Infection_Coefficient = coef_totalinf,
    Tri1_Infection_Coefficient = coef_tri1inf,
    Tri2_Infection_Coefficient = coef_tri2inf,
    Tri3_Infection_Coefficient = coef_tri3inf
  )
  
  # Append the results to the main data frame
  asthma_results_df_Rosa <- rbind(asthma_results_df_Rosa, result_row)
}

asthma_results_df_Rosa2 <- t(asthma_results_df_Rosa)
asthma_results_df_Rosa2 <- as.data.frame(asthma_results_df_Rosa2)

write_xlsx(asthma_results_df_Rosa2, path = "replication_mps_child_phenotypes-asthma_Rosa_240320.xlsx")

#\

# Second: linear regression for BMI and behavior
child_pheno2 <- c('depression_f15', 'anxiety_f15', 'adhd_f15', 'cd_f15', 'odd_f15',  'bmi_f17') 

rosa_results_bmi_behavior_df <- data.frame() 

for(x in child_pheno2) {
  # lm formula
  a <- paste0(x, "~ scale(MPS_total_infections) + m_age + GENDER + msmoke + EDUCM + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + 
              X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20")
  
  b <- paste0(x, "~ scale(MPS_tri1_infections) + m_age + GENDER + msmoke + EDUCM + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + 
              X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20")
  
  c <- paste0(x, "~ scale(MPS_tri2_infections) + m_age + GENDER + msmoke + EDUCM + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + 
              X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20")
  
  d <- paste0(x, "~ scale(MPS_tri3_infections) + m_age + GENDER + msmoke + EDUCM + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + 
              X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20")
  
  #extract coefficients 
  coef_totalinf <- summary(lm(as.formula(a), data = mps_pheno_imp30_svs_df))$coefficients[2,]
  coef_tri1inf <- summary(lm(as.formula(b), data = mps_pheno_imp30_svs_df))$coefficients[2,]
  coef_tri2inf <- summary(lm(as.formula(c), data = mps_pheno_imp30_svs_df))$coefficients[2,]
  coef_tri3inf <- summary(lm(as.formula(d), data = mps_pheno_imp30_svs_df))$coefficients[2,]
  
  #Create a data frame with the results
  result_row <- data.frame(
    Phenotype = x,
    Total_Infection_Coefficient = coef_totalinf,
    Tri1_Infection_Coefficient = coef_tri1inf,
    Tri2_Infection_Coefficient = coef_tri2inf,
    Tri3_Infection_Coefficient = coef_tri3inf
  )
  
  # Append the results to the main data frame
  rosa_results_bmi_behavior_df <- rbind(rosa_results_bmi_behavior_df, result_row)
}

results_df2 <- t(rosa_results_bmi_behavior_df)
results_df2 <- as.data.frame(results_df2)

write_xlsx(results_df2, path = "replication_mps_child_phenotypes-bmi-behavior_Rosa_240320.xlsx")

# Scaled lm's
child_pheno2 <- c('depression_f15', 'anxiety_f15', 'adhd_f15', 'cd_f15', 'odd_f15',  'bmi_f17') 

rosa_results_bmi_behavior_scaled_df <- data.frame() 

for(x in child_pheno2) {
  # lm formula
  a <- paste0("scale(",x,")", "~ scale(MPS_total_infections) + m_age + GENDER + msmoke + EDUCM + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + 
              X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20")
  
  b <- paste0("scale(",x,")", "~ scale(MPS_tri1_infections) + m_age + GENDER + msmoke + EDUCM + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + 
              X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20")
  
  c <- paste0("scale(",x,")", "~ scale(MPS_tri2_infections) + m_age + GENDER + msmoke + EDUCM + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + 
              X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20")
  
  d <- paste0("scale(",x,")", "~ scale(MPS_tri3_infections) + m_age + GENDER + msmoke + EDUCM + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + 
              X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20")
  
  #extract coefficients 
  coef_totalinf <- summary(lm(as.formula(a), data = mps_pheno_imp30_svs_df))$coefficients[2,]
  coef_tri1inf <- summary(lm(as.formula(b), data = mps_pheno_imp30_svs_df))$coefficients[2,]
  coef_tri2inf <- summary(lm(as.formula(c), data = mps_pheno_imp30_svs_df))$coefficients[2,]
  coef_tri3inf <- summary(lm(as.formula(d), data = mps_pheno_imp30_svs_df))$coefficients[2,]
  
  #Create a data frame with the results
  result_row <- data.frame(
    Phenotype = x,
    Total_Infection_Coefficient = coef_totalinf,
    Tri1_Infection_Coefficient = coef_tri1inf,
    Tri2_Infection_Coefficient = coef_tri2inf,
    Tri3_Infection_Coefficient = coef_tri3inf
  )
  
  # Append the results to the main data frame
  rosa_results_bmi_behavior_scaled_df <- rbind(rosa_results_bmi_behavior_scaled_df, result_row)
}

results_df2_scaled <- t(rosa_results_bmi_behavior_scaled_df)
results_df2_scaled <- as.data.frame(results_df2_scaled)

write_xlsx(results_df2_scaled, path = "replication_mps_child_phenotypes-bmi-behavior-SCALED_Rosa_240320.xlsx")

#\ END OF THIS SCRIPT.
