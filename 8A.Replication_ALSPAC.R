#########################################################################
#####################PROJECT: PRENATAL INFECTION EWAS####################
#########################################################################

#Author: Anna Suleri

#In this script we will replicate the analyses from Generation R in the ALSPAC cohort. This script includes the following parts:
#' Step 1: Select variables of interest in dataframe and prepare dataframe for our specific project (also create infection sum score), apply inclusion and exclusion criteria, create baseline table, apply multiple imputation
#' Step 2: Prepare methylation data and create gestational clocks 
#' Step 3: Run replication analyses:
#' 3.1. Replicate top hits (suggestive significance) identified in GenR by associating them with infections
#' 3.2. Dmrff with ewas 450k results genR and methylation data alspac? (NA - no sign findings in GenR)
#' 3.3. Create MPS with weights from GenR-train set
#' 3.4. Associate MPS-infection with infections (lm, incremental R2)
#' 3.5. Associate MPS-infection with child phenotypes 
#' 3.6 Associate infection with gestational age clocks 
#' 
#' Run this on GENA server 

#############---------Step 1: Prepare pheno dataframe---------#############
### Set wd & load libraries
set.seed(2023)

libraries <- c('foreign', 'tidyverse', 'devtools', 'ALSPAC.helpR', 'mice', 'readxl', 'tibble', 'writexl', 'data.table', 'haven', 'dplyr', 'MASS', 'Matrix', 'DescTools', "meffil", 'pROC', 'plotly')

invisible(lapply(libraries, require, character.only = T)) 

setwd('set_path_to_data')

### Read in data
phenotype.data <- read.spss('path_to_alspac_data', to.data.frame = T)

### Build phenotype dataframe  
## First select all variables we need for all the analyses 
pheno <- dplyr::select(phenotype.data, c('cidB3361','b053', 'b055', 'c059', 'c057', 'e106', 'e105', 'kz021', 'e171', 'e173', 'e175', 'e176', 'b670', 'e178', 'bestgest', 'd994', 'qlet', 'c642', 'c630', 'c631', "c632" ,'c641', 'kb611', 'fh6876', 'fh6878', 'fh6870', 'fh6874', 'fh6873', 'cf062', 'cf064', 'cf066', 'cf069', 'f7ms026a', 'f9ms026a', 'fdms026a', 'fems026a', 'FJMR022a', 'FKMS1040'))

saveRDS(pheno, 'pheno.rds')

## Re categorize data and check structure of all variables 
## Create maternal smoking variable
pheno$msmoke <- NA
pheno$msmoke[pheno$e171 == "Y" | pheno$e173 == "Y" |  pheno$e175 == "Y" | pheno$e176 == "Y"] <- "continued"
pheno$msmoke[pheno$b670 != 0  & pheno$e178 == "Not at all"] <- "stopped"
pheno$msmoke[pheno$e171 == "N" & pheno$e173 == "N" & pheno$e175 == "N" & pheno$e176 == "N" & pheno$b670 == 0 & pheno$e178 == "Not at all"] <- "no smoking"
pheno$msmoke <- as.factor(pheno$msmoke)
#summary(pheno$msmoke)
#str(pheno$msmoke)

pheno <- dplyr::select(pheno, -c('e171', 'e173', 'e175', 'e176', 'b670', 'e178')) #drop columns we don't need anymore 

## Create continuous maternal education variable
pheno$EDUCM <- NA
pheno$EDUCM[pheno$c642 == "Yes"]  <- 0 #No
pheno$EDUCM[pheno$c630 == "Yes"] <- 1 #CSE
pheno$EDUCM[pheno$c631 == "Yes" | pheno$c632 == "Yes"] <- 2 #O-level
pheno$EDUCM[pheno$c641 == "Yes"]  <- 3 #University
#summary(pheno$EDUCM)
#str(pheno$EDUCM)

pheno <- dplyr::select(pheno, -c('c642', 'c630', "c632" ,'c631', 'c641')) #drop columns we don't need anymore 

## Create dimensions for infection sum score (but not yet sum up, that we will do after imputation)
pheno$influenza_tri1 <- as.factor(ifelse(pheno$b055 == 'Y in 1-3 MTHS', 1, ifelse(pheno$b055 == 'Y 4 MTHS to now' | pheno$b055 == 'not at all', 0, NA)))
pheno$influenza_tri2a <- as.factor(ifelse(pheno$b055 == 'Y 4 MTHS to now', 1, ifelse(pheno$b055 == 'Y in 1-3 MTHS' |  pheno$b055 == 'not at all', 0, NA)))
pheno$influenza_tri2b <- as.factor(ifelse(pheno$c059 == 'Y', 1, ifelse(pheno$c059 == 'N', 0, NA)))
pheno$influenza_tri2 <- as.factor(ifelse(pheno$influenza_tri2a == 1 | pheno$influenza_tri2b == 1, 1, ifelse(pheno$influenza_tri2a == 0 & pheno$influenza_tri2b == 0, 0, NA)))
pheno$influenza_tri3 <- as.factor(ifelse(pheno$e106 == 'Y', 1, ifelse(pheno$e106 == 'N', 0, NA))) 

pheno <- dplyr::select(pheno, -c('b055', 'c059', 'e106', 'influenza_tri2a', 'influenza_tri2b')) #drop columns we don't need anymore 

pheno$uti_tri1 <- as.factor(ifelse(pheno$b053 == 'Y in 1-3 MTHS', 1, ifelse(pheno$b053 == 'Y 4 MTHS to now' | pheno$b053 == 'not at all', 0, NA)))
pheno$uti_tri2a <- as.factor(ifelse(pheno$b053 == 'Y 4 MTHS to now', 1, ifelse(pheno$b053 == 'Y in 1-3 MTHS' |  pheno$b053 == 'not at all', 0, NA)))
pheno$uti_tri2b <- as.factor(ifelse(pheno$c057 == 'Y', 1, ifelse(pheno$c057 == 'N', 0, NA)))
pheno$uti_tri2 <- as.factor(ifelse(pheno$uti_tri2a == 1 | pheno$uti_tri2b == 1, 1, ifelse(pheno$uti_tri2a == 0 & pheno$uti_tri2b == 0, 0, NA)))

pheno$uti_tri3 <- as.factor(ifelse(pheno$e105 == 'Y', 1, ifelse(pheno$e105 == 'N', 0, NA))) 

pheno <- dplyr::select(pheno, -c('b053', 'c057', 'e105', 'uti_tri2a', 'uti_tri2b')) #drop columns we don't need anymore 

## Rename variables
pheno <- pheno %>% rename('GENDER' = kz021, 'm_age' = d994, 'depression_f15' = fh6876, 'anxiety_f15' = fh6878, 'adhd_f15' = fh6870, 'cd_f15' = fh6874, 'odd_f15' = fh6873, 'bmi_f15' = FJMR022a,'twin' = kb611)

## Check structure of vars and adapt if needed  
str(pheno)

pheno$bestgest <- as.numeric(as.character(pheno$bestgest))
pheno$m_age <- as.numeric(as.character(pheno$m_age))
pheno$bmi_f15 <- as.numeric(as.character(pheno$bmi_f17))
pheno$depression_f15 <- as.numeric(pheno$depression_f15)
pheno$anxiety_f15 <- as.numeric(pheno$anxiety_f15)
pheno$adhd_f15 <- as.numeric(pheno$adhd_f15)
pheno$cd_f15 <- as.numeric(pheno$cd_f15)
pheno$odd_f15 <- as.numeric(pheno$odd_f15)

saveRDS(pheno, 'final_pheno.rds')

### Merge link file 
link.data <- read.spss("path_to_data", to.data.frame=T)
phenotype_link.data <- merge(pheno, link.data, by = c("cidB3361", "qlet"))

### Merge methylation meta data
samples <- read.csv("~/Alspac_Tempo/methylation/data/samplesheet/samplesheet.csv")

# Filter for birth
samples_birth.data <- samples[samples$time_point == "cord", ]

phenotype_link_sample.data <- merge(phenotype_link.data, samples_birth.data, by = c("dnam_epic450_g0_g1", 'qlet'))

# Merge cell count data
cell_counts_combined.data <- read.csv("/data/Alspac_Tempo/methylation/data/derived/cellcounts/combined-cord-blood.txt", sep = "")
phenotype_link_sample_counts.data <- merge(phenotype_link_sample.data, cell_counts_combined.data, by.x = "Sample_Name", by.y = "IID")
# dim(phenotype_link_sample_counts.data)

### Apply exclusion criteria
# Remove twins 
df1 <- subset(phenotype_link_sample_counts.data, twin == 'No') 

# Remove 1 sibling
df2 <- df1[!duplicated(df1$cidB3361),]

# check final participants 
dim(df2) 

saveRDS(df2, 'df_after_exclusion.rds')

### Create baseline table 
baselinevars <- c('m_age', 'EDUCM', 'msmoke', 'bestgest', 'GENDER')

summary_continuous <- function(x)
{
  standev <- sd(x, na.rm = T)
  meanvar <- mean(x, na.rm = T)
  print(paste0(round(meanvar, 1), '(', round(standev, 1), ')'))
}

summary_categorical <- function(x)
{
  tab1 <- prop.table(table(x, useNA = 'always'))
  tab2 <- table(x, useNA = "always")
  print(paste(round(tab1 * 100, 1), '%', names(tab1), collapse = ','))
  print(paste(tab2, names(tab2)))
}

for(i in baselinevars){
  x <- df2[, i]
  message(i)
  if (class(x) == 'numeric') {
    summary_continuous(x)
  } else { summary_categorical(x) }
}

table(df2$EDUCM)

### Multiple imputation of exposure and covariates
## First check missingness 
miss_values_function <- function(df){
  missvalues <- cbind("# NA" = sort(colSums(is.na(df))),
                      "% NA" = round(sort(colMeans(is.na(df))) * 100, 2))
  print(missvalues)
}

miss_values_function(df2) 

## Impute missing exposure and covars
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
  imp.test <- mice(data, method = meth, predictorMatrix = pred, visitSequence = visSeq, maxit = 30, m = 30, printFlag = TRUE, seed = 2023)  
  
  #get the dataframe name 
  dataName <- deparse(substitute(data))
  
  #assigning a dynamic name to the imp.test object
  imputedDataName <- paste0("imputedData_", dataName)
  assign(imputedDataName, imp.test)
  
  # Saving the imputed dataset as .RDS file
  saveRDS(imp.test, file = paste0(imputedDataName, ".rds"))
  
  # Output 
  return(list(imp0 = imp0, imp.test = imp.test))
}

imputation_function(data = df2, exclude_imp_vars = c(1:4, 8:23, 32:70), exclude_predictors = c(1:2, 4, 14:17,19:23 ,32:70), method = "rf")

# Check quality imputation 
imp_df <- readRDS('imputedData_df2.rds')

plot(imp_df)

# Select last imputed df
df_final <- complete(imp_df, 30)

# Create total infection sum scores
df_final$infection_trimester1 <- as.numeric(ifelse(df_final$influenza_tri1 == 1 | df_final$uti_tri1 == 1, 1, 0))
df_final$infection_trimester2 <- as.numeric(ifelse(df_final$influenza_tri2 == 1 | df_final$uti_tri2 == 1, 1, 0))  
df_final$infection_trimester3 <- as.numeric(ifelse(df_final$influenza_tri3 == 1 | df_final$uti_tri3 == 1, 1, 0)) 

df_final$total_infections <- df_final$infection_trimester1 + df_final$infection_trimester2 + df_final$infection_trimester3

# Save final df 
saveRDS(df_final, 'df_final_imp.rds')

#\ 

#############---------Step 2: Prepare methylation dataframe---------#############
### Load and merge in methylation data 

## Load birth methylation data ALSPAC
dnam_alspac_orig.data <- meffil.gds.methylation("~/Alspac_Tempo/methylation/data/betas/450.gds")

## Keep only birth samples
samples <- read.csv("~/Alspac_Tempo/methylation/data/samplesheet/samplesheet.csv")

## Filter for birth
samples_birth.data <- samples[samples$time_point == "cord", ]

## Filter the 450K methylation data for birth
dnam_alspac_cord.data <- dnam_alspac_orig.data[ ,samples_birth.data$Sample_Name]

dim(dnam_alspac_cord.data) #rows = cpg's and columns = id 

## Transpose (now cpgs will be columns and ids will be rows)
dnam_alspac_cord <- as.data.frame(t(dnam_alspac_cord.data))
dim(dnam_alspac_cord)

## Save dataset
saveRDS(dnam_alspac_cord, 'dnam_alspac_cord.rds')

#\ 

#############---------Step 3: Run replication analyses---------#############
### 3.1. Replicate top hits (suggestive significance) identified in GenR by associating them with infections

## First, load GenR EWAS sumstats (nominal significant results; cross-reactive cpgs are removed)
# set wd
setwd('set_path_to_data_MPS_genR')

# load data
library(readxl)
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

## Load data set alspac 
setwd('set_path_to_data_alspac')

# First load methylation data
dnam_alspac <- readRDS('dnam_alspac_cord.rds')
dnam_alspac2 <-  tibble::rownames_to_column(dnam_alspac, "Sample_Name")

# only select cpg sites of interest
cpgs <- c(sigCpGs_totalinf, sigCpGs_trimester1, sigCpGs_trimester2, sigCpGs_trimester3)

dnam_alspac3 <- dnam_alspac2 %>%
  select(Sample_Name, all_of(cpgs))

# winsorize cpg
dnam_alspac3_winsorized <- dnam_alspac3 %>%
  mutate(across(-1, ~ Winsorize(., probs = c(0.05, 0.95), na.rm = TRUE)))

# Load imputed final pheno data
pheno <- readRDS('df_final_imp.rds')

# Merge methylation dataframe with pheno 
alspac_df <- merge(pheno, dnam_alspac3_winsorized, by = 'Sample_Name', all.x = T)

## Run replication analysis
# Check str dataframe
str(alspac_df)
alspac_df$qlet <- as.factor(alspac_df$qlet)

# Replication total infections 
# Of note, not adjusting for parity because no variation for all df (i.e., everyone is qlet=A)
setwd('set_path_to_results')

results_df_totalinf <- data.frame() # Create an empty data frame to store results

for (x in sigCpGs_totalinf){
  # specify lm
  m1 <- paste0(x,"~ scale(total_infections) + m_age + GENDER + msmoke + EDUCM + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + plate")
  m2 <- paste0(x,"~ scale(total_infections) + m_age + GENDER + msmoke + EDUCM + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + plate + bestgest")
  # extract coefficients
  coef_totalinf_m1 <- summary(lm(as.formula(m1), data = alspac_df))$coefficients[2,]
  coef_totalinf_m2 <- summary(lm(as.formula(m2), data = alspac_df))$coefficients[2,]
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

write_xlsx(results_df_totalinf, path = "replication_results_top_hits_ewas_totalinf.xlsx")

# Replication trimester 1 infections
results_df_trimester1 <- data.frame()

for (x in sigCpGs_trimester1){
  # specify lm
  m1 <- paste0(x,"~ scale(infection_trimester1) + m_age + GENDER + msmoke + EDUCM + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + plate")
  m2 <- paste0(x,"~ scale(infection_trimester1) + m_age + GENDER + msmoke + EDUCM + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + plate + bestgest")
  # extract coefficients
  coef_totalinf_m1 <- summary(lm(as.formula(m1), data = alspac_df))$coefficients[2,]
  coef_totalinf_m2 <- summary(lm(as.formula(m2), data = alspac_df))$coefficients[2,]
  # save output 
  results_df_trimester1[x,c(1:4)] <- coef_totalinf_m1
  results_df_trimester1[x,c(5:8)] <- coef_totalinf_m2
}

colnames(results_df_trimester1) <- c("bval_m1", "seval_m1", "tval_m1","pval_m1","bval_m2", "seval_m2", "tval_m2","pval_m2")

results_df_trimester1$outcome <- sigCpGs_trimester1

results_df_trimester1$sign_m1 = "" 
results_df_trimester1$sign_m1[results_df_trimester1$pval_m1 < 0.05] = "*"
results_df_trimester1$sign_m1[results_df_trimester1$pval_m1 < 0.01] = "**"
results_df_trimester1$sign_m1[results_df_trimester1$pval_m1 < 0.001] = "***"
results_df_trimester1$sign_m2 = "" 
results_df_trimester1$sign_m2[results_df_trimester1$pval_m2 < 0.05] = "*"
results_df_trimester1$sign_m2[results_df_trimester1$pval_m2 < 0.01] = "**"
results_df_trimester1$sign_m2[results_df_trimester1$pval_m2 < 0.001] = "***"

write_xlsx(results_df_trimester1, path = "replication_results_top_hits_ewas_tri1.xlsx")

# Replication trimester 2 infections
results_df_trimester2 <- data.frame()

for (x in sigCpGs_trimester2){
  # specify lm
  m1 <- paste0(x,"~ scale(infection_trimester2) + m_age + GENDER + msmoke + EDUCM + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + plate")
  m2 <- paste0(x,"~ scale(infection_trimester2) + m_age + GENDER + msmoke + EDUCM + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + plate + bestgest")
  # extract coefficients
  coef_totalinf_m1 <- summary(lm(as.formula(m1), data = alspac_df))$coefficients[2,]
  coef_totalinf_m2 <- summary(lm(as.formula(m2), data = alspac_df))$coefficients[2,]
  # save output 
  results_df_trimester2[x,c(1:4)] <- coef_totalinf_m1
  results_df_trimester2[x,c(5:8)] <- coef_totalinf_m2
}

colnames(results_df_trimester2) <- c("bval_m1", "seval_m1", "tval_m1","pval_m1","bval_m2", "seval_m2", "tval_m2","pval_m2")

results_df_trimester2$outcome <- sigCpGs_trimester2

results_df_trimester2$sign_m1 = "" 
results_df_trimester2$sign_m1[results_df_trimester2$pval_m1 < 0.05] = "*"
results_df_trimester2$sign_m1[results_df_trimester2$pval_m1 < 0.01] = "**"
results_df_trimester2$sign_m1[results_df_trimester2$pval_m1 < 0.001] = "***"
results_df_trimester2$sign_m2 = "" 
results_df_trimester2$sign_m2[results_df_trimester2$pval_m2 < 0.05] = "*"
results_df_trimester2$sign_m2[results_df_trimester2$pval_m2 < 0.01] = "**"
results_df_trimester2$sign_m2[results_df_trimester2$pval_m2 < 0.001] = "***"

write_xlsx(results_df_trimester2, path = "replication_results_top_hits_ewas_tri2.xlsx")

# Replication trimester 3 infections
results_df_trimester3 <- data.frame()

for (x in sigCpGs_trimester3){
  # specify lm
  m1 <- paste0(x,"~ scale(infection_trimester3) + m_age + GENDER + msmoke + EDUCM + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + plate")
  m2 <- paste0(x,"~ scale(infection_trimester3) + m_age + GENDER + msmoke + EDUCM + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + plate + bestgest")
  # extract coefficients
  coef_totalinf_m1 <- summary(lm(as.formula(m1), data = alspac_df))$coefficients[2,]
  coef_totalinf_m2 <- summary(lm(as.formula(m2), data = alspac_df))$coefficients[2,]
  # save output 
  results_df_trimester3[x,c(1:4)] <- coef_totalinf_m1
  results_df_trimester3[x,c(5:8)] <- coef_totalinf_m2
}

colnames(results_df_trimester3) <- c("bval_m1", "seval_m1", "tval_m1","pval_m1","bval_m2", "seval_m2", "tval_m2","pval_m2")

results_df_trimester3$outcome <- sigCpGs_trimester3

results_df_trimester3$sign_m1 = "" 
results_df_trimester3$sign_m1[results_df_trimester3$pval_m1 < 0.05] = "*"
results_df_trimester3$sign_m1[results_df_trimester3$pval_m1 < 0.01] = "**"
results_df_trimester3$sign_m1[results_df_trimester3$pval_m1 < 0.001] = "***"
results_df_trimester3$sign_m2 = "" 
results_df_trimester3$sign_m2[results_df_trimester3$pval_m2 < 0.05] = "*"
results_df_trimester3$sign_m2[results_df_trimester3$pval_m2 < 0.01] = "**"
results_df_trimester3$sign_m2[results_df_trimester3$pval_m2 < 0.001] = "***"

write_xlsx(results_df_trimester3, path = "replication_results_top_hits_ewas_tri3.xlsx")

#\ 

### 3.3. Create MPS with weights from GenR-train set
## Step 1: Construct MPS based on ENR-weights from GenR. The base data = weights/beta from elastic net regression in GENR & target Data: Methylation at birth in ALSPAC. 

# 1)load GENR ENR-weight file 
setwd('set_path_to_data')
S1 <- read_xlsx("ENR_weights_total_infection_GenR_suggesivepval.xlsx")
S2 <- read.csv("ENR_weights_trimester1_GenR_suggesivepval.csv")
S3 <- read.csv("ENR_weights_trimester2_GenR_suggesivepval.csv")
S4 <- read.csv("ENR_weights_trimester3_GenR_suggesivepval.csv")

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

# MPS based on ENR-weight2
wmatrix2[1:5,1:5]
dim(wmatrix2)
colsum2 = colSums(wmatrix2, na.rm=TRUE)
mps_S2 = as.data.frame(colsum2)
head(mps_S2)

# MPS based on ENR-weight3
wmatrix3[1:5,1:5]
dim(wmatrix3)
colsum3 = colSums(wmatrix3, na.rm=TRUE)
mps_S3 = as.data.frame(colsum3)
head(mps_S3)

# MPS based on ENR-weight4
wmatrix4[1:5,1:5]
dim(wmatrix4)
colsum4 = colSums(wmatrix4, na.rm=TRUE)
mps_S4 = as.data.frame(colsum4)
head(mps_S4)

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
setwd('set_path_to_results')

MPS.infect <- Reduce(merge, list(mps_S1, mps_S2, mps_S3, mps_S4))

saveRDS(MPS.infect, "MPS_infections_ALSPAC.rds")

#\ 

### 3.4. Associate MPS-infection with infections (lm, incremental R2)

## Step 1: merge pheno file with mps file
mps_pheno_df <- merge(pheno, MPS.infect, by.x = 'Sample_Name', by.y = 'Sample_ID', all.x = T)
dim(mps_pheno_df)
saveRDS(mps_pheno_df, "mps_pheno_df.rds")

#mps_pheno_df <- readRDS('mps_pheno_df.rds')

## Create descriptive figures 
# 1. Create histograms of MPS
library(plotly)

plot_ly(data = mps_pheno_df, x = ~MPS_total_infections, type = "histogram", marker = list(color = "purple")) %>%
  layout(title = "A. Histogram of MPS total infection",
         xaxis = list(title = "MPS total infection"),
         yaxis = list(title = "Frequency"))

plot_ly(data = mps_pheno_df, x = ~MPS_tri1_infections, type = "histogram", marker = list(color = "purple")) %>%
  layout(title = "B. Histogram of MPS trimester 1",
         xaxis = list(title = "MPS trimester 1"),
         yaxis = list(title = "Frequency"))

plot_ly(data = mps_pheno_df, x = ~MPS_tri2_infections, type = "histogram", marker = list(color = "purple")) %>%
  layout(title = "C. Histogram of MPS trimester 2",
         xaxis = list(title = "MPS trimester 2"),
         yaxis = list(title = "Frequency"))

plot_ly(data = mps_pheno_df, x = ~MPS_tri3_infections, type = "histogram", marker = list(color = "purple")) %>%
  layout(title = "D. Histogram of MPS trimester 3",
         xaxis = list(title = "MPS trimester 3"),
         yaxis = list(title = "Frequency"))


# 2. Run regressions between total infections and child outcomes
mps_pheno_df$bmi_f17 <- as.numeric(mps_pheno_df$bmi_f17)

child_pheno <- c('bmi_f17') 

results_df <- data.frame() 

#Initialize an empty dataframe to store results
results_df <- data.frame(Outcome = character(), stringsAsFactors = FALSE)

# Loop through each child phenotype
for (x in child_pheno) {
  # lm formula
  formula <- as.formula(paste0("scale(", x, ")", " ~ scale(total_infections) + m_age + GENDER + msmoke + EDUCM + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + plate"))
  
  # Fit linear model
  lm_model <- lm(formula, data = mps_pheno_df)
  
  # Extract coefficients 
  coef_totalinf <- summary(lm_model)$coefficients[2,]
  
  # Create a data frame with the results
  result_row <- data.frame(
    Outcome = x,
    t(coef_totalinf),  # Transpose to have coefficients as columns
    stringsAsFactors = FALSE
  )
  
  # Append the results to the main data frame
  results_df <- bind_rows(results_df, result_row)
}

# Print the results dataframe
print(results_df)

write_xlsx(results_df, path = "results_lm_infections_child_outcomes.xlsx")

## Step 2: linear regression with infections

# Create plots of linear regresion for total infections 
scatter <- ggplot(mps_pheno_df, aes(x = total_infections, y = scale(MPS_total_infections))) +
  geom_point(color = "#bc5090") +
  labs(title = "Regression Plot", x = "Prenatal infection sum score", y = "Methylation profile score of infections")

reg_line <- geom_smooth(data = mps_pheno_df, method = "lm", se = FALSE, color = "#58508d")

plot <- ggplotly(scatter + reg_line)
plot <- plot %>% layout(plot_bgcolor = "lightblue")
plot 

# 2.1 Run linear regression
lm_totalinf <- summary(lm(scale(MPS_total_infections) ~  scale(total_infections) + m_age + GENDER + msmoke + EDUCM + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + plate, data =mps_pheno_df))$coefficients[2,] 

lm_tri1 <- summary(lm(scale(MPS_tri1_infections) ~  scale(infection_trimester1) + m_age + GENDER + msmoke + EDUCM + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + plate, data =mps_pheno_df))$coefficients[2,] 

lm_tri2 <- summary(lm(scale(MPS_tri2_infections) ~  scale(infection_trimester2) + m_age + GENDER + msmoke + EDUCM + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + plate, data =mps_pheno_df))$coefficients[2,] 

lm_tri3 <- summary(lm(scale(MPS_tri3_infections) ~  scale(infection_trimester3) + m_age + GENDER + msmoke + EDUCM + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + plate, data =mps_pheno_df))$coefficients[2,] 

library(tibble) # create df with results to later save to excel 
lm_total <- rbind(lm_totalinf, lm_tri1, lm_tri2, lm_tri3)
lm_total <- as.data.frame(lm_total)
lm_total <-  tibble::rownames_to_column(lm_total, "Exposure")

# 2.2. Calculate incremental r2 --> this indicates prediction peformance of MPS 
# mod 1=basic model, including all covariates in your EWAS
mod1 <- lm(scale(total_infections) ~ m_age + GENDER + msmoke + EDUCM + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + plate, mps_pheno_df)
r.mod1 = summary(mod1)$r.squared   

mod2 <- lm(scale(total_infections) ~ scale(MPS_total_infections) + m_age + GENDER + msmoke + EDUCM + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + plate, mps_pheno_df)
r.mod2 = summary(mod2)$r.squared
r2change =  r.mod2 - r.mod1 # 0.0009437336

mod3 <- lm(scale(infection_trimester1) ~ m_age + GENDER + msmoke + EDUCM + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + plate, mps_pheno_df)
r.mod3 = summary(mod3)$r.squared   

mod4 <- lm(scale(infection_trimester1) ~ scale(MPS_tri1_infections) + m_age + GENDER + msmoke + EDUCM + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + plate, mps_pheno_df)
r.mod4 = summary(mod4)$r.squared
r2change2 =  r.mod4 - r.mod3 #0.0002250439

mod5 <- lm(scale(infection_trimester2) ~m_age + GENDER + msmoke + EDUCM + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + plate, mps_pheno_df)
r.mod5 = summary(mod5)$r.squared   

mod6 <- lm(scale(infection_trimester2) ~ scale(MPS_tri2_infections) + m_age + GENDER + msmoke + EDUCM + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + plate, mps_pheno_df)
r.mod6 = summary(mod6)$r.squared
r2change3 =  r.mod6 - r.mod5 # 0.001708174

mod7 <- lm(scale(infection_trimester3) ~ m_age + GENDER + msmoke + EDUCM + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + plate, mps_pheno_df)
r.mod7 = summary(mod7)$r.squared   

mod8 <- lm(scale(infection_trimester3) ~ scale(MPS_tri3_infections) +m_age + GENDER + msmoke + EDUCM + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + plate, mps_pheno_df)
r.mod8 = summary(mod8)$r.squared
r2change4 =  r.mod8 - r.mod7 #2.650397e-05

lm_total$r2_change <- c(r2change, r2change2, r2change3, r2change4) # add results to lm validation df 

setwd('/home/a.suleri/infection_EWAS_project/replication_ALSPAC/results')
write_xlsx(lm_total, 'replication_results_lm_infections.xlsx')

#\ 

### 3.5. Associate MPS-infection with child phenotypes 
## Step 1: SDQ scales 
fl("sdq", data = phenotype.data)

#SDQ total difficulties, SDQ emotional symptoms, SDQ hyperactivity, SDQ conduct problems = tc4025a, b,c,f
sdq_vars <- dplyr::select(phenotype.data, c("cidB3361","tc4025a", 'tc4025b',"tc4025c",'tc4025f'))
sdq_vars$tc4025a <- as.numeric(as.character(sdq_vars$tc4025a))
sdq_vars$tc4025b <- as.numeric(as.character(sdq_vars$tc4025b))
sdq_vars$tc4025c <- as.numeric(as.character(sdq_vars$tc4025c))
sdq_vars$tc4025f <- as.numeric(as.character(sdq_vars$tc4025f))

mps_pheno_df2 <- merge(mps_pheno_df, sdq_vars, by = "cidB3361", all.x =T)

mps_pheno_df2 <- rename(mps_pheno_df2, "SDQ_total_difficulties" = tc4025a, "SDQ_emotional_symptoms" = tc4025b, "SDQ_hyperactivity" = tc4025c, "SDQ_conduct_problems" = tc4025f)

#saveRDS(mps_pheno_df2,'mps_pheno_df2.rds')

child_pheno_sdq <- c('SDQ_total_difficulties', 'SDQ_emotional_symptoms', 'SDQ_hyperactivity', 'SDQ_conduct_problems') 

results_df <- data.frame() 

for(x in child_pheno_sdq) {
  # lm formula
  a <- paste0("scale(",x,")", "~ scale(MPS_total_infections) + m_age + GENDER + msmoke + EDUCM + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + plate")
  
  b <- paste0("scale(",x,")", "~ scale(MPS_tri1_infections) + m_age + GENDER + msmoke + EDUCM + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + plate")
  
  c <- paste0("scale(",x,")", "~ scale(MPS_tri2_infections) + m_age + GENDER + msmoke + EDUCM + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + plate")
  
  d <- paste0("scale(",x,")", "~ scale(MPS_tri3_infections) + m_age + GENDER + msmoke + EDUCM + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + plate")
  
  #extract coefficients 
  coef_totalinf <- summary(lm(as.formula(a), data = mps_pheno_df2))$coefficients[2,]
  coef_tri1inf <- summary(lm(as.formula(b), data = mps_pheno_df2))$coefficients[2,]
  coef_tri2inf <- summary(lm(as.formula(c), data = mps_pheno_df2))$coefficients[2,]
  coef_tri3inf <- summary(lm(as.formula(d), data = mps_pheno_df2))$coefficients[2,]
  
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

results_df2 <- t(results_df)
results_df2 <- as.data.frame(results_df2)

write_xlsx(results_df2, path = "replication_mps_child_phenotypes_sdq.xlsx")

## step 2: Create correlation plot for Alspac variables
cor_vars <- c("GENDER", "bestgest", "m_age","msmoke", "EDUCM", "bmi_f17", "infection_trimester1", "infection_trimester2", "infection_trimester3", "total_infections", "MPS_total_infections", "MPS_tri1_infections", "MPS_tri2_infections", "MPS_tri3_infections", 'SDQ_total_difficulties', 'SDQ_emotional_symptoms', 'SDQ_hyperactivity', 'SDQ_conduct_problems')

selected_df <- mps_pheno_df2[cor_vars]

# Convert variables to numeric
selected_df <- as.data.frame(lapply(selected_df, as.numeric))
selected_df <- subset(selected_df, complete.cases(selected_df))
selected_df <- rename(selected_df, "13. Child sex" = GENDER, "9. Gestational age at birth" = bestgest, "10. Maternal age" = m_age, "18. BMI (age 17)" = bmi_f17, "2. Infection sum score (trimester 1)" = infection_trimester1, "3. Infection sum score (trimester 2)" = infection_trimester2, "4. Infection sum score (trimester 3)" = infection_trimester3, "1. Total infection sum score" = total_infections, "5. MPS (total infections)" = MPS_total_infections, "6. MPS (trimester 1)" = MPS_tri1_infections, "7. MPS (trimester 2)" = MPS_tri2_infections, "8. MPS (trimester 3)" = MPS_tri3_infections, "12. Maternal smoking" = msmoke, "11. Maternal education" = EDUCM, '14. SDQ - total difficulties' = SDQ_total_difficulties, '15. SDQ - emotional symptoms' = SDQ_emotional_symptoms, '16. SDQ - hyperactivity' = SDQ_hyperactivity, '17. SDQ - conduct problems' = SDQ_conduct_problems)

# Create a correlation plot
library(corrplot)

correlation <- cor(selected_df, use="pairwise.complete.obs")

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
  "14. SDQ - total difficulties",
  "15. SDQ - emotional symptoms",
  "16. SDQ - hyperactivity",
  "17. SDQ - conduct problems",
  "18. BMI (age 17)"
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

### 3.6 Associate winsorized MPS-infection with infection & child phenotypes 
#Winsorize all MPS-infection variables 
#Our goal is to transform extreme values with less extreme values to reduce impact of outliers
library(DescTools)
mps_pheno_df$MPS_total_infections_win <- Winsorize(mps_pheno_df$MPS_total_infections, probs = c(0.05, 0.95), na.rm = T)
mps_pheno_df$MPS_tri1_infections_win <- Winsorize(mps_pheno_df$MPS_tri1_infections, probs = c(0.05, 0.95), na.rm = T)
mps_pheno_df$MPS_tri2_infections_win <- Winsorize(mps_pheno_df$MPS_tri2_infections, probs = c(0.05, 0.95), na.rm = T)
mps_pheno_df$MPS_tri3_infections_win <- Winsorize(mps_pheno_df$MPS_tri3_infections, probs = c(0.05, 0.95), na.rm = T)

# plot histograms 
plot_ly(data = mps_pheno_df, x = ~MPS_total_infections_win, type = "histogram", marker = list(color = "purple")) %>%
  layout(title = "A. Histogram of MPS total infection (winsorized)",
         xaxis = list(title = "MPS total infection"),
         yaxis = list(title = "Frequency"))

plot_ly(data = mps_pheno_df, x = ~MPS_tri1_infections_win, type = "histogram", marker = list(color = "purple")) %>%
  layout(title = "B. Histogram of MPS trimester 1 (winsorized)",
         xaxis = list(title = "MPS trimester 1"),
         yaxis = list(title = "Frequency"))

plot_ly(data = mps_pheno_df, x = ~MPS_tri2_infections_win, type = "histogram", marker = list(color = "purple")) %>%
  layout(title = "C. Histogram of MPS trimester 2 (winsorized)",
         xaxis = list(title = "MPS trimester 2"),
         yaxis = list(title = "Frequency"))

plot_ly(data = mps_pheno_df, x = ~MPS_tri3_infections_win, type = "histogram", marker = list(color = "purple")) %>%
  layout(title = "D. Histogram of MPS trimester 3 (winsorized)",
         xaxis = list(title = "MPS trimester 3"),
         yaxis = list(title = "Frequency"))

# Run lm with infections
lm_totalinf_win <- summary(lm(scale(MPS_total_infections_win) ~  scale(total_infections) + m_age + GENDER + msmoke + EDUCM + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + plate, data =mps_pheno_df))$coefficients[2,] 

lm_tri1_win <- summary(lm(scale(MPS_tri1_infections_win) ~  scale(infection_trimester1) + m_age + GENDER + msmoke + EDUCM + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + plate, data =mps_pheno_df))$coefficients[2,] 

lm_tri2_win <- summary(lm(scale(MPS_tri2_infections_win) ~  scale(infection_trimester2) + m_age + GENDER + msmoke + EDUCM + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + plate, data =mps_pheno_df))$coefficients[2,] 

lm_tri3_win <- summary(lm(scale(MPS_tri3_infections_win) ~  scale(infection_trimester3) + m_age + GENDER + msmoke + EDUCM + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + plate, data =mps_pheno_df))$coefficients[2,] 

lm_total_win <- rbind(lm_totalinf_win, lm_tri1_win, lm_tri2_win, lm_tri3_win)

# create plot of lm 
scatter <- ggplot(mps_pheno_df, aes(x = total_infections, y = scale(MPS_total_infections_win))) +
  geom_point(color = "#bc5090") +
  labs(title = "Regression Plot", x = "Prenatal infection sum score", y = "Methylation profile score of infections (winsorized)")

reg_line <- geom_smooth(data = mps_pheno_df, method = "lm", se = FALSE, color = "#58508d")

plot <- ggplotly(scatter + reg_line)
plot <- plot %>% layout(plot_bgcolor = "lightblue")
plot 

# calculate incremental r2
mod1_win <- lm(scale(total_infections) ~ m_age + GENDER + msmoke + EDUCM + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + plate, mps_pheno_df)
r.mod1_win = summary(mod1_win)$r.squared   

mod2_win <- lm(scale(total_infections) ~ scale(MPS_total_infections_win) + m_age + GENDER + msmoke + EDUCM + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + plate, mps_pheno_df)
r.mod2_win = summary(mod2_win)$r.squared
r2change_win =  r.mod2_win - r.mod1_win 

mod3_win <- lm(scale(infection_trimester1) ~ m_age + GENDER + msmoke + EDUCM + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + plate, mps_pheno_df)
r.mod3_win = summary(mod3_win)$r.squared   

mod4_win <- lm(scale(infection_trimester1) ~ scale(MPS_tri1_infections_win) + m_age + GENDER + msmoke + EDUCM + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + plate, mps_pheno_df)
r.mod4_win = summary(mod4_win)$r.squared
r2change2_win =  r.mod4_win - r.mod3_win 

mod5_win <- lm(scale(infection_trimester2) ~m_age + GENDER + msmoke + EDUCM + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + plate, mps_pheno_df)
r.mod5_win = summary(mod5_win)$r.squared   

mod6_win <- lm(scale(infection_trimester2) ~ scale(MPS_tri2_infections_win) + m_age + GENDER + msmoke + EDUCM + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + plate, mps_pheno_df)
r.mod6_win = summary(mod6_win)$r.squared
r2change3_win =  r.mod6_win - r.mod5_win 

mod7_win <- lm(scale(infection_trimester3) ~ m_age + GENDER + msmoke + EDUCM + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + plate, mps_pheno_df)
r.mod7_win = summary(mod7_win)$r.squared   

mod8_win <- lm(scale(infection_trimester3) ~ scale(MPS_tri3_infections_win) +m_age + GENDER + msmoke + EDUCM + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + plate, mps_pheno_df)
r.mod8_win = summary(mod8_win)$r.squared
r2change4_win =  r.mod8_win - r.mod7_win 

r2_change <- c(0.003134761, 3.47514e-06, 0.0004123876, 2.468865e-07)
lm_total_win2 <- cbind(lm_total_win2, r2_change)

setwd('/home/a.suleri/infection_EWAS_project/replication_ALSPAC/results')
lm_total_win2 <- as.data.frame(lm_total_win2)
write_xlsx(lm_total_win2, 'replication_results_lm_infections_winsorized_mps.xlsx')

#\ 

### 3.7 Associate infection with gestational age clocks 
## Step 1) Create gestational age clocks 
## We want: Bohlin & Epic - raw & residual form of EGAA 
rm(list=ls())

set.seed(2023)

# Step 1.1: load libraries 
libraries <- c('Biobase', 'tibble', 'impute', 'ggplot2', 'ggpmisc', 'GEOquery', 'preprocessCore', 'foreign', 'tidyverse', 'meffil', 'devtools', "MASS", "robustbase")
invisible(lapply(libraries, require, character.only = T)) 

setwd('/home/a.suleri/infection_EWAS_project/replication_ALSPAC/')
devtools::load_all('methylclock-main') #of note, this package is downloaded via code > download in zip file > unpack > and then stored on server (https://github.com/KristinaSalontaji/methylclock)

# Step 1.2 Load methylation data
dnam_alspac_orig.data <- meffil.gds.methylation("~/Alspac_Tempo/methylation/data/betas/450.gds")
samples <- read.csv("~/Alspac_Tempo/methylation/data/samplesheet/samplesheet.csv")
samples_birth.data <- samples[samples$time_point == "cord", ]
dnam_alspac_cord.data <- dnam_alspac_orig.data[ ,samples_birth.data$Sample_Name]
dim(dnam_alspac_cord.data) #rows = cpg's and columns = id 
meth <- as.data.frame(dnam_alspac_cord.data)
meth[1:5, 1:5]

meth$cpg <- rownames(meth)
meth2 <- dplyr::select(meth, c(914, 1:913))

# Step 1.3 Check available cpgs
cpgs.missing.GA <- checkClocksGA(meth2) #no missings, so no need to impute 

# Step 1.4 Load imputed final pheno data
pheno <- readRDS('df_final_imp.rds')

# Step 1.5 Select cpg's of interest and align dimensions with pheno dataframe file 
est_Bohlin <- meth[rownames(meth) %in% coefBohlin$CpGmarker[-1],]
est_Bohlin2 <- est_Bohlin[colnames(est_Bohlin) %in% pheno$Sample_Name]
est_Bohlin3 <- as.data.frame(t(est_Bohlin2))
est_Bohlin4 <- as.matrix(est_Bohlin3)

est_EPIC <- meth[rownames(meth) %in% coefEPIC$CpGmarker[-1],]
est_EPIC2 <- est_EPIC[colnames(est_EPIC) %in% pheno$Sample_Name]
est_EPIC3 <- as.data.frame(t(est_EPIC2))
est_EPIC4 <- as.matrix(est_EPIC3)

# Step 1.6 Calculate raw & residual acceleration
# first calculate clock estimate 
min.perc = 0.5
bohlin <- predAge(est_Bohlin4, coefBohlin, intercept = TRUE, min.perc)

bohlin <- data.frame(
  id = rownames(est_Bohlin4),
  Bohlin = bohlin / 7
)

epic <- predAge(est_EPIC4, coefEPIC, intercept = TRUE, min.perc)

epic <- data.frame(
  id = rownames(est_EPIC4),
  EPIC = epic / 7
)

# then calculate raw and residual acceleration 
Bohlin <- ageAcc1(est_Bohlin3, pheno$bestgest, lab = "Bohlin")
Epic <- ageAcc1(est_EPIC3, pheno$bestgest, lab = "EPIC")

Bohlin2 <- dplyr::select(Bohlin, c('ageAcc', 'ageAcc2'))
Epic2 <- dplyr::select(Epic, c('ageAcc', 'ageAcc2'))

setwd('set_path_to_results')
saveRDS(Bohlin,"Bohlin_clock_alspac.rds")
saveRDS(Epic,"Epic_clock_alspac.rds")

# plot clock estimate and residual accelration
# looks good
plotDNAmAge(Bohlin$ageAcc2, bohlin$Bohlin, 
            clock="GA")

plotDNAmAge(Epic$ageAcc2, epic$EPIC, 
            clock="GA")

# Step 1.8 Combine phenotype data and clocks 
pheno_clocks_df_bohlin <- cbind(pheno, Bohlin2)
pheno_clocks_df_epic <- cbind(pheno, Epic2)

## Step 2) lm between infections and clocks (residual estimate of clocks)
# 2.1 create vector of exposure
infection <- c("total_infections", "infection_trimester1", "infection_trimester2", "infection_trimester3")

# 2.2 Create lists to store results
results_list_bohlin <- list()
results_list_epic <- list()

# 2.3 run replication lm for total infections
for (exposure in infection) {
    
    # Model formulas
    m1 <- paste0("scale(ageAcc2) ~ scale(", exposure, ") + m_age + GENDER + msmoke + EDUCM + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + plate")
    m2 <- paste0("scale(ageAcc2) ~ scale(", exposure, ") + m_age + GENDER + msmoke + EDUCM + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + plate + bestgest")
    
    # Extract coefficients
    coef_totalinf_m1b <- summary(lm(as.formula(m1), data = pheno_clocks_df_bohlin))$coefficients[2,]
    coef_totalinf_m2b <- summary(lm(as.formula(m2), data = pheno_clocks_df_bohlin))$coefficients[2,]
    
    coef_totalinf_m1e <- summary(lm(as.formula(m1), data = pheno_clocks_df_epic))$coefficients[2,]
    coef_totalinf_m2e <- summary(lm(as.formula(m2), data = pheno_clocks_df_epic))$coefficients[2,]
    
    # Create a data frame for the current iteration
    result_df_m1 <- data.frame(exposure = exposure, coef_totalinf_m1b, coef_totalinf_m2b)
    result_df_m2 <- data.frame(exposure = exposure, coef_totalinf_m1e, coef_totalinf_m2e)
    
    # Append to the results list
    results_list_bohlin[[length(results_list_bohlin) + 1]] <- result_df_m1
    results_list_epic[[length(results_list_epic) + 1]] <- result_df_m2
  }

# Combine the results data frames
results_clocks_bohlin <- do.call(rbind, results_list_bohlin)
results_clocks_epic <- do.call(rbind, results_list_epic)

setwd('set_path_to_results')
write_xlsx(results_clocks_bohlin, "results_ALSPAC_lm_GA_clocks_bohlin.xlsx")
write_xlsx(results_clocks_epic, "results_ALSPAC_lm_GA_clocks_epic.xlsx")

## Step 3) robust lm between infections and clocks (residual estimate of clocks) 

inc <- lmrob.control() #increase number of iterations as lmrob may fail to converge in rare cases
inc$maxit.scale<-2500 #you may increase this number if a model fails to converge (you will get a warning saying that it fails to converge)
inc$k.max <-2500 #same here
set.seed(2023)

library(broom.mixed)

model0 <- with(pheno_clocks_df_bohlin, lmrob(scale(ageAcc2) ~ scale(total_infections) + m_age + GENDER + msmoke + EDUCM + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + plate, control = inc, trace.lev=4))
summary <- tidy(model0, conf.int = TRUE, conf.level = 0.95)

model <- rlm(scale(ageAcc2) ~ scale(total_infections) + m_age + GENDER + msmoke + EDUCM + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + plate, data = pheno_clocks_df_bohlin)

# similar to lm and because plate leads to problems with lmrob, we stick to lm

#\ END OF SCRIPT 
