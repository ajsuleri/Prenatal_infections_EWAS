Author: Kristina Salontaji (k.salontaji@erasmusmc.nl)

####################################################################################################################################################################
# The following R code will allow you to calculate epigenetic gestational age for the Knight, Bohlin, and EPIC overlap clocks, 
# along with raw ("ageACC") and residual ("ageAcc2") epigenetic age acceleration 
# You shouldn't have to rewrite or add to the following code, unless otherwise stated.
# The current script calculates Generation R data.
####################################################################################################################################################################

#library(methylclock)
library(Biobase)
library(tibble)
library(impute)
library(ggplot2)
library(ggpmisc)
library(GEOquery)
library(preprocessCore)
library(foreign)
library(tidyverse)
library(minfi)
#library(minfiData)
library(dplyr)
library(meffonym)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(devtools)

#download the adapted methylclock package from github and extract it into your working directory: https://github.com/KristinaSalontaji/methylclock
devtools::load_all('~/methylclock-main') #load local adapted methylclock package

#setwd("your working dir") #change to your own working directory 

#load gestational age file. This file contains: participant ID (IDC), sample ID (ID used in our 450k array), gestational age, cell type proportions (bcell, cd4t, cd8t, gran, mono, nk, nRBC)
ageacc <- readRDS(file = "EPICageacc.rds")
#ageacc$sample.id <- gsub(" ", "", ageacc$sample.id)
colnames(ageacc)[colnames(ageacc) == "SampleID"] <- "Sample_ID"
colnames(ageacc)[colnames(ageacc) == "GESTBIR"] <- "GA"
#ageacc$GA <- ageacc$GA 

GA <- ageacc$GA
age = GA
min.perc = 0.5


genR_DNAm = load("data") # matrix with cpgs as rows and subjects as cols ; paste your own file location here; do this for both 450k methylation data as well as epic data seperately 
genR_DNAm <- GENR_EPICv1METH_Norm_Betas_birth_ALL.data
GENR_EPICv1METH_Norm_Betas_birth_ALL.data <- NULL
gc()

backupdata <- genR_DNAm

genR_DNAm[1:5,1:5]



mat <- backupdata
meth <- as.data.frame(mat)

meth[1:5, 1:5]


#epic_numeric <- as.numeric(unlist(genR_DNAm))
epic_matrix <- do.call(cbind, genR_DNAm)
epic_matrix[1:5, 1:5]
rownames(epic_matrix) <- rownames(genR_DNAm)

cpgs.missing.GA <- checkClocksGA(epic_matrix) #no missings, so no need to impute


# Find rows with row names containing "."
rows_to_remove <- grep("\\.", rownames(epic_matrix))

# Remove rows by index
epic_matrix <- epic_matrix[-rows_to_remove, , drop = FALSE]

epic_matrix[1:5, 1:5]

#turn rownames to column
#meth <- rownames_to_column(meth, var = "cpgsite") 
#epic_matrix <- do.call(cbind, meth)
#epic_matrix[1:5, 1:5]



merged_matrix1 <- cbind(imputed_subset$data, imputed_subset2$data, imputed_subset3$data) #,
#imputed_subset4$data,imputed_subset5$data,imputed_subset6$data)
#save(merged_matrix1, file ="imputationoutput_epic1.Rdata") ##########load this and run calculations separately = call it cpgimpute.
cpgs.missing.GA <- checkClocksGA(merged_matrix1) #no missings, so no need to impute
merged1_transposed <- t(merged_matrix1)
merged1_transposed <- t(epic_matrix)


bohlin <- predAge(merged1_transposed, coefBohlin, intercept = TRUE, min.perc)
Bohlin <- data.frame(id = rownames(merged1_transposed), Bohlin = bohlin / 7)

epic <- predAge(merged1_transposed, coefEPIC, intercept = TRUE, min.perc)
EPIC <- data.frame(id = rownames(merged1_transposed), EPIC = epic / 7)


#calculate raw acceleration per clock
Bohlin <- ageAcc1(Bohlin, age, lab = "Bohlin")
EPIC <- ageAcc1(EPIC, age, lab = "EPIC")



#join the files from the different clocks together
out <- Bohlin %>%
  full_join(EPIC, by = "id")

out <- tibble::as_tibble(out)

age = GA
out <- add_column(out, age = age) #add age column

ageAcc_data_epic = out

save(ageAcc_data_epic, file = "ageAcc_data_epic.rda") 



load("ageAcc_data_epic.rda")
out = ageAcc_data_epic

cpg.names <- rownames(epic_matrix)
#checkEPIC <- coefEPIC$CpGmarker[-1][!coefEPIC$CpGmarker[-1] %in% cpg.names]

save(checkEPIC, file = "checkEPIC.rda")

# we plot the clocks against gestational age at birth to see how well they perform in each dataset
# I am running the epigenetic clock calculations on an external server and will import the ageAcc_data.rda file together with the following graphs to a new folder that will allow me to use RStudio. 
# In the next part of this analysis We will create a new folder called 

pdf(file="Knight.pdf") # save the plot 

dev.copy(png,'Knight.png')
plotDNAmAge(out$Knight, out$age, 
            tit="GA Knight", 
            clock="GA") + 
  theme_classic() +
  ggplot2::xlim(25, 45) +
  ggplot2::ylim(25, 45)
dev.off()


pdf(file="Bohlin_epic.pdf") # save the plot 
dev.copy(png,'Bohlin.png')
plotDNAmAge(out$Bohlin, out$age,  
             tit="GA Bohlin",
             clock="GA") + 
  theme_classic() +
  ggplot2::xlim(25, 45) +
  ggplot2::ylim(25, 45)
dev.off()


pdf(file="EPIC_epic.pdf") # save the plot 
dev.copy(png,'EPIC.png')
plotDNAmAge(out$EPIC, out$age, 
            tit="GA EPIC overlap clock", 
            clock="GA") + 
  theme_classic() +
  ggplot2::xlim(25, 45) +
  ggplot2::ylim(25, 45)

dev.off()

#\ END OF THIS SCRIPT.
