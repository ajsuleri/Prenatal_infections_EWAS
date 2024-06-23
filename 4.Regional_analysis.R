#########################################################################
#####################PROJECT: PRENATAL INFECTION EWAS####################
#########################################################################

# Author: Anna Suleri

# In this script we perform the following steps:
# Perform a regional analysis for model 1, first for all individual arrays and then meta-analyze the results. We have four exposures: total infection and trimester 1-3 infections.
# So we will have 4 exposures x 2 arrays = 8 regional analyses that we will later meta-analyze
# Run part 1 and 2 on gena server and run part 3 on local r studio 

### Step 1: preparation 
## set WD
setwd('set_path_to_data')

## Load libraries of interest 
library(dplyr)
library(dmrff)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

## Annotate individual ewases to get chr and pos

# Load annotation file
annotation <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# load files
load('overlap_ewas_450k_model1.rdata')
ewas_450k_total_infection <- overlap_ewas_450k_model1.data
load('overlap_ewas_450k_trimester1_model1.rdata')
overlap_ewas_450k_trimester1 <- overlap_ewas_450k_trimester1_model1.data
load('overlap_ewas_450k_trimester2_model1.rdata')
overlap_ewas_450k_trimester2 <- overlap_ewas_450k_trimester2_model1.data
load('overlap_ewas_450k_trimester3_model1.rdata')
overlap_ewas_450k_trimester3 <- overlap_ewas_450k_trimester3_model1.data

load('overlap_ewas_epic_model1.rdata')
ewas_epic_total_infection <- overlap_ewas_epic_model1.data
load('overlap_ewas_epic_trimester1_model1.rdata')
overlap_ewas_epic_trimester1 <- overlap_ewas_epic_trimester1_model1.data
load('overlap_ewas_epic_trimester2_model1.rdata')
overlap_ewas_epic_trimester2 <- overlap_ewas_epic_trimester2_model1.data
load('overlap_ewas_epic_trimester3_model1.rdata')
overlap_ewas_epic_trimester3 <- overlap_ewas_epic_trimester3_model1.data

# Merge each rdata file with annotation 
an_ewas_450k_total_infection <-merge(ewas_450k_total_infection, annotation,by.x = "cpg", by.y = "Name", all.x = T) 
an_ewas_450k_trimester1_infection <-merge(overlap_ewas_450k_trimester1, annotation,by.x = "cpg", by.y = "Name", all.x = T)
an_ewas_450k_trimester2_infection <-merge(overlap_ewas_450k_trimester2, annotation,by.x = "cpg", by.y = "Name", all.x = T)
an_ewas_450k_trimester3_infection <-merge(overlap_ewas_450k_trimester3, annotation,by.x = "cpg", by.y = "Name", all.x = T)

an_ewas_epic_total_infection <-merge(ewas_epic_total_infection, annotation,by.x = "cpg", by.y = "Name", all.x = T) 
an_ewas_epic_trimester1_infection <-merge(overlap_ewas_epic_trimester1, annotation,by.x = "cpg", by.y = "Name", all.x = T)
an_ewas_epic_trimester2_infection <-merge(overlap_ewas_epic_trimester2, annotation,by.x = "cpg", by.y = "Name", all.x = T)
an_ewas_epic_trimester3_infection <-merge(overlap_ewas_epic_trimester3, annotation,by.x = "cpg", by.y = "Name", all.x = T)

# Save the merged dataframe as a CSV file in the specified output folder
output_folder <- 'path_to_results'

# Create the output folder if it doesn't exist
dir.create(output_folder, showWarnings = FALSE)

# Construct the full path for the CSV file
csv_file_path <- file.path(output_folder, 'an_ewas_450k_total_infection.csv')
csv_file_path2 <- file.path(output_folder, 'an_ewas_450k_trimester1_infection.csv')
csv_file_path3 <- file.path(output_folder, 'an_ewas_450k_trimester2_infection.csv')
csv_file_path4 <- file.path(output_folder, 'an_ewas_450k_trimester3_infection.csv')

csv_file_path5 <- file.path(output_folder, 'an_ewas_epic_total_infection.csv')
csv_file_path6 <- file.path(output_folder, 'an_ewas_epic_trimester1_infection.csv')
csv_file_path7 <- file.path(output_folder, 'an_ewas_epic_trimester2_infection.csv')
csv_file_path8 <- file.path(output_folder, 'an_ewas_epic_trimester3_infection.csv')

# Save the dataframe as a CSV file
write.csv(an_ewas_450k_total_infection, file = csv_file_path, row.names = FALSE)
write.csv(an_ewas_450k_trimester1_infection, file = csv_file_path2, row.names = FALSE)
write.csv(an_ewas_450k_trimester2_infection, file = csv_file_path3, row.names = FALSE)
write.csv(an_ewas_450k_trimester3_infection, file = csv_file_path4, row.names = FALSE)

write.csv(an_ewas_epic_total_infection, file = csv_file_path5, row.names = FALSE)
write.csv(an_ewas_epic_trimester1_infection, file = csv_file_path6, row.names = FALSE)
write.csv(an_ewas_epic_trimester2_infection, file = csv_file_path7, row.names = FALSE)
write.csv(an_ewas_epic_trimester3_infection, file = csv_file_path8, row.names = FALSE)

# Print a message indicating where the file is saved (sanity check)
cat("CSV file saved at:", csv_file_path, "\n")
cat("CSV file saved at:", csv_file_path2, "\n")
cat("CSV file saved at:", csv_file_path3, "\n")
cat("CSV file saved at:", csv_file_path4, "\n")
cat("CSV file saved at:", csv_file_path5, "\n")
cat("CSV file saved at:", csv_file_path6, "\n")
cat("CSV file saved at:", csv_file_path7, "\n")
cat("CSV file saved at:", csv_file_path8, "\n")

#\ 

### Step 2: Run DMR analysis (meta-analyze across assays)
## http://htmlpreview.github.io/?https://github.com/perishky/dmrff/blob/master/docs/meta-analysis.html 
## load and prepare dataframes
setwd("/path_to_results")

csv_files <- list.files(pattern = "\\.csv$")

for (file in csv_files) {
  name <- tools::file_path_sans_ext(file)
  assign(name, read.csv(file))
}

## load methylation data
setwd('path_to_data')
overlap_df_450k <- readRDS("overlap_df_450k.rds")
methylation <- as.matrix(t(overlap_df_450k))
overlap_df_epic <- readRDS("overlap_df_epic.rds")
methylation2 <- as.matrix(t(overlap_df_epic))

## align order of cpg's between ewas output and methylation output  !!!
methylation_450k <- methylation[an_ewas_450k_total_infection$cpg,]
methylation_epic <- methylation2[an_ewas_epic_total_infection$cpg,]

#double check by looking at first 6 cpg's if order is correct
head(methylation_450k)
head(an_ewas_450k_total_infection)

head(methylation_epic)
head(an_ewas_epic_total_infection)

## Run dmrff meta analysis + create and save tables 

# (first run dmrff on individual dataset (create pre objects), then apply meta functino, then select fdr sig regions, then use dmrff.sites function to show CpG sites in each DMR and their meta-analyzed sum stats, and create tables with all info we want)

#' total infections 
total_inf_450k <- dmrff.pre(an_ewas_450k_total_infection$Estimate, an_ewas_450k_total_infection$'Std..Error',methylation_450k, an_ewas_450k_total_infection$chr, an_ewas_450k_total_infection$pos)

total_inf_epic <- dmrff.pre(an_ewas_epic_total_infection$Estimate, an_ewas_epic_total_infection$'Std..Error',methylation_epic, an_ewas_epic_total_infection$chr, an_ewas_epic_total_infection$pos)

options(mc.cores=20)
pre_totalinf <- list(total_inf_450k, total_inf_epic)
meta_totalinf <- dmrff.meta(pre_totalinf)

dmrs_totalinf <- meta_totalinf$dmrs[which(meta_totalinf$dmrs$p.adjust < 0.09 & meta_totalinf$dmrs$n >= 2), ]
dmrs_totalinf2 <- dmrs_totalinf[,c("chr","start","end","n","estimate","se","p.value","p.adjust")]

setwd('path_to_results')
#write_xlsx(dmrs_totalinf2, 'dmrff_total_infections_nocpgsites.xlsx') # 0 findings 

#sites_totalinf <- dmrff.sites(dmrs_totalinf2, meta_totalinf$ewas$chr, meta_totalinf$ewas$pos) # here we use the dmrff.sites function to show the CpG sites in each DMR and their meta-analysed summary statistics.

#library(tibble)
#meta_totalinf$ewas <- tibble::rownames_to_column(meta_totalinf$ewas, "cpg")

#sites_totalinf2 <- cbind(sites_totalinf[,c("region", "site" ,"chr","pos")], meta_totalinf$ewas[sites_totalinf$site,c("cpg", "estimate","se","p.value")]) 

#write_xlsx(sites_totalinf2, 'dmrff-MA_totalinfections.xlsx')

#' trimester 1
tri1_inf_450k <- dmrff.pre(an_ewas_450k_trimester1_infection$Estimate, an_ewas_450k_trimester1_infection$'Std..Error',methylation_450k, an_ewas_450k_trimester1_infection$chr, an_ewas_450k_trimester1_infection$pos)

tri1_inf_epic <- dmrff.pre(an_ewas_epic_trimester1_infection$Estimate, an_ewas_epic_trimester1_infection$'Std..Error',methylation_epic, an_ewas_epic_trimester1_infection$chr, an_ewas_epic_trimester1_infection$pos)

options(mc.cores=20)
pre_tri1_inf <- list(tri1_inf_450k, tri1_inf_epic)
meta_tri1_inf <- dmrff.meta(pre_tri1_inf)

dmrs_tri1_inf <- meta_tri1_inf$dmrs[which(meta_tri1_inf$dmrs$p.adjust < 0.05 & meta_tri1_inf$dmrs$n >= 2), ]
dmrs_tri1_inf <- dmrs_tri1_inf[,c("chr","start","end","n","estimate","se", "p.value","p.adjust")]

#write_xlsx(dmrs_tri1_inf, 'dmrff_trimester1_infections_nocpgsites.xlsx') #0 findings

#sites_tri1_inf <- dmrff.sites(dmrs_tri1_inf, meta_tri1_inf$ewas$chr, meta_tri1_inf$ewas$pos) 

#meta_tri1_inf$ewas <- tibble::rownames_to_column(meta_tri1_inf$ewas, "cpg")

#sites_tri1_inf <- cbind(sites_tri1_inf[,c("region", "site" ,"chr","pos")], meta_tri1_inf$ewas[sites_tri1_inf$site,c("cpg", "estimate","se","p.value")]) 

#write_xlsx(sites_tri1_inf, 'dmrff-MA_tri1_infections.xlsx')

#' trimester 2
tri2_inf_450k <- dmrff.pre(an_ewas_450k_trimester2_infection$Estimate, an_ewas_450k_trimester2_infection$'Std..Error',methylation_450k, an_ewas_450k_trimester2_infection$chr, an_ewas_450k_trimester2_infection$pos)

tri2_inf_epic <- dmrff.pre(an_ewas_epic_trimester2_infection$Estimate, an_ewas_epic_trimester2_infection$'Std..Error',methylation_epic, an_ewas_epic_trimester2_infection$chr, an_ewas_epic_trimester2_infection$pos)

options(mc.cores=20)
pre_tri2_inf <- list(tri2_inf_450k, tri2_inf_epic)
meta_tri2_inf <- dmrff.meta(pre_tri2_inf)

dmrs_tri2_inf <- meta_tri2_inf$dmrs[which(meta_tri2_inf$dmrs$p.adjust < 0.05 & meta_tri2_inf$dmrs$n >= 2), ]
dmrs_tri2_inf <- dmrs_tri2_inf[,c("chr","start","end","n","estimate","se", "p.value","p.adjust")]

#write_xlsx(dmrs_tri2_inf, 'dmrff_trimester2_infections_nocpgsites.xlsx') #0 findings

#sites_tri2_inf <- dmrff.sites(dmrs_tri2_inf, meta_tri2_inf$ewas$chr, meta_tri2_inf$ewas$pos) 

#meta_tri2_inf$ewas <- tibble::rownames_to_column(meta_tri2_inf$ewas, "cpg")

#sites_tri2_inf <- cbind(sites_tri2_inf[,c("region", "site" ,"chr","pos")], meta_tri2_inf$ewas[sites_tri2_inf$site,c("cpg", "estimate","se","p.value")]) 

#write_xlsx(sites_tri2_inf, 'dmrff-MA_tri2_infections.xlsx')

#' trimester 3
tri3_inf_450k <- dmrff.pre(an_ewas_450k_trimester3_infection$Estimate, an_ewas_450k_trimester3_infection$'Std..Error',methylation_450k, an_ewas_450k_trimester3_infection$chr, an_ewas_450k_trimester3_infection$pos)

tri3_inf_epic <- dmrff.pre(an_ewas_epic_trimester3_infection$Estimate, an_ewas_epic_trimester3_infection$'Std..Error',methylation_epic, an_ewas_epic_trimester3_infection$chr, an_ewas_epic_trimester3_infection$pos)

options(mc.cores=20)
pre_tri3_inf <- list(tri3_inf_450k, tri3_inf_epic)
meta_tri3_inf <- dmrff.meta(pre_tri3_inf)

dmrs_tri3_inf <- meta_tri3_inf$dmrs[which(meta_tri3_inf$dmrs$p.adjust < 0.05 & meta_tri3_inf$dmrs$n >= 2), ]
dmrs_tri3_inf <- dmrs_tri3_inf[,c("chr","start","end","n","estimate","se","p.value","p.adjust")]

#write_xlsx(dmrs_tri3_inf, 'dmrff_trimester3_infections_nocpgsites.xlsx') #0 findings

#sites_tri3_inf <- dmrff.sites(dmrs_tri3_inf, meta_tri3_inf$ewas$chr, meta_tri3_inf$ewas$pos) 

#meta_tri3_inf$ewas <- tibble::rownames_to_column(meta_tri3_inf$ewas, "cpg")

#sites_tri3_inf <- cbind(sites_tri3_inf[,c("region", "site" ,"chr","pos")], meta_tri3_inf$ewas[sites_tri3_inf$site,c("cpg", "estimate","se","p.value")]) 

#write_xlsx(sites_tri3_inf, 'dmrff-MA_tri3_infections.xlsx')

#\ END OF SCRIPT
