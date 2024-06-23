#########################################################################
#####################PROJECT: PRENATAL INFECTION EWAS####################
#########################################################################

# Author: Anna Suleri

# This script is part 3 for this project in which it is the goal to do run the eQTM and gene ontology analysis as functional analyses.
# We do both these analyses for the suggestive hits 
# We do that for each exposure but only for the main model, namely model 1 

##########------------Part 1: eQTM------------##########
## First load significant results and only select suggestive significance hits (this is for the post processing dataset)
setwd("set_path_to_cleaned_results_folder")

library(readxl)

df_tri1 <- read_excel("final_cleaned_sign_results_meta_analyzed_results_trimester1_infection_model1.xlsx")
df_tri2 <- read_excel("final_cleaned_sign_results_meta_analyzed_results_trimester2_infection_model1.xlsx")
df_tri3 <- read_excel("final_cleaned_sign_results_meta_analyzed_results_trimester3_infection_model1.xlsx")
df_tot <- read_excel("final_sign_results_meta_analyzed_results_total_infection_model1.xlsx")

## Load in helix database for eQTM
setwd("set_path_to_functional_list")
eqtms_uniq <- readRDS("eQTM_autosome_adj.cells_SIG_unique.RData")

## Merge datasets with helix databse to identify eQTM's
df_tri1$eqtm <- ifelse(df_tri1$cpg %in% eqtms_uniq$CpG, "yes", "no")
df_tri2$eqtm <- ifelse(df_tri2$cpg %in% eqtms_uniq$CpG, "yes", "no")
df_tri3$eqtm <- ifelse(df_tri3$cpg %in% eqtms_uniq$CpG, "yes", "no")
df_tot$eqtm <- ifelse(df_tot$cpg %in% eqtms_uniq$CpG, "yes", "no")

## Select only suggestive significance hits and save datasets
setwd("set_path_to_results_folder")

# First create dataframe with only hits that are below p=5e-05 (suggestive significance threshold)
library(tidyverse)

df_tri1_eqtm <- subset(df_tri1, meta_pval < 0.00005)
df_tri2_eqtm <- subset(df_tri2, meta_pval < 0.00005)
df_tri3_eqtm <- subset(df_tri3, meta_pval < 0.00005)
df_tot_eqtm <- subset(df_tot, meta_pval < 0.00005)

# look up for the sig cpgs which eqtm gene is associated
subset(eqtms_uniq, CpG == "cg00702872")
subset(eqtms_uniq, CpG == "cg01304814")
subset(eqtms_uniq, CpG == "cg08337633")
subset(eqtms_uniq, CpG == "cg03987884")

# save datasets
library(xlsx)

write.xlsx(df_tri1_eqtm, 'df_tri1_eqtm.xlsx')
write.xlsx(df_tri2_eqtm, 'df_tri2_eqtm.xlsx')
write.xlsx(df_tri3_eqtm, 'df_tri3_eqtm.xlsx')
write.xlsx(df_tot_eqtm, 'df_tot_eqtm.xlsx')

#\

##########------------Part 2: Gene ontology analysis------------##########
# run this part of the analysis on the gena server ! 

# missMethyl: Analysing Illumina HumanMethylation BeadChip Data
# https://bioconductor.org/packages/devel/bioc/vignettes/missMethyl/inst/doc/missMethyl.html

# The GO terms of the BP, CC and MF categories enrichment if testing GO pathways
# Biological Process (BP) 	
# Cellular Component (CC)
# Molecular Function (MF)
# N	= number of genes in the GO or KEGG term
# DE	= number of genes that are differentially methylated
# P.DE	 p-value for over-representation of the GO or KEGG term term

## load libraries
library(missMethyl)  
library(data.table)
library(tibble)

## set wd
setwd('set_path_to_data')

## load & prepare data 
# load all data

# List all CSV files in the current directory
csv_files <- list.files(pattern = "\\.csv$")

# Loop through each CSV file, read it, and assign it to an object with the original name
for (file in csv_files) {
  # Extract the file name without the file extension
  name <- tools::file_path_sans_ext(file)
  
  # Read the CSV file and assign it to an object with the original name
  assign(name, read.csv(file))
}

# create vector for all cpgs 
# also check length 
all_tot <- merged_meta_analyzed_results_total_infection_model1$cpg
length(all_tot)
all_tri1 <- merged_meta_analyzed_results_trimester1_infection_model1$cpg
length(all_tri1)
all_tri2 <- merged_meta_analyzed_results_trimester2_infection_model1$cpg
length(all_tri2)
all_tri3 <- merged_meta_analyzed_results_trimester3_infection_model1$cpg
length(all_tri3)

# create vector for suggestive significant hits
# also check length 
sig_cpgs_tot <- merged_meta_analyzed_results_total_infection_model1$cpg[merged_meta_analyzed_results_total_infection_model1$meta_pval<0.00005]
length(sig_cpgs_tot)
sig_cpgs_tri1 <- merged_meta_analyzed_results_trimester1_infection_model1$cpg[merged_meta_analyzed_results_trimester1_infection_model1$meta_pval<0.00005]
length(sig_cpgs_tri1)
sig_cpgs_tri2 <- merged_meta_analyzed_results_trimester2_infection_model1$cpg[merged_meta_analyzed_results_trimester2_infection_model1$meta_pval<0.00005]
length(sig_cpgs_tri2)
sig_cpgs_tri3 <- merged_meta_analyzed_results_trimester3_infection_model1$cpg[merged_meta_analyzed_results_trimester3_infection_model1$meta_pval<0.00005]
length(sig_cpgs_tri3)

## apply GO analysis
#This code is performing a gene ontology analysis using the gometh function, considering the differentially methylated CpG sites (sigCpGs) compared to all CpG sites (all). The analysis is based on the Gene Ontology database (collection = "GO"), and it includes the option to generate a plot to visualize potential bias (plot.bias = TRUE). The results of the analysis are stored in the variable gst.
gst_tot <- gometh(sig.cpg=sig_cpgs_tot, all.cpg= all_tot, collection="GO", plot.bias=TRUE)
gst_tri1 <- gometh(sig.cpg=sig_cpgs_tri1, all.cpg= all_tri1, collection="GO", plot.bias=TRUE)
gst_tri2 <- gometh(sig.cpg=sig_cpgs_tri2, all.cpg= all_tri2, collection="GO", plot.bias=TRUE)
gst_tri3 <- gometh(sig.cpg=sig_cpgs_tri3, all.cpg= all_tri3, collection="GO", plot.bias=TRUE)

# check results of significant enrichment after fdr correction
# result is no enrichment 
table(gst_tot$FDR < .05)
table(gst_tri1$FDR < .05)
table(gst_tri2$FDR < .05)
table(gst_tri3$FDR < .05)

# hence, create a table for top 20 hits 
setwd('set_path_to_results')
library(writexl)

tot_top20<-topGSA(gst_tot, n=20)
write_xlsx(tot_top20, "Go_tot_top20_5x10-5.xlsx")

tri1_top20<-topGSA(gst_tri1, n=20)
write_xlsx(tri1_top20, "Go_tri1_top20_5x10-5.xlsx")

tri2_top20<-topGSA(gst_tri2, n=20)
write_xlsx(tri2_top20, "Go_tri2_top20_5x10-5.xlsx")

tri3_top20<-topGSA(gst_tri3, n=20) # interesting enrichment --> vascular processes 
write_xlsx(tri3_top20, "Go_tri3_top20_5x10-5.xlsx")

## apply KEGG to interpret biological meaning of a set of genes 
# apply keg, then check fdr sig results, and because there are none, select top 20 and save table
tot_gst.kegg <- gometh(sig.cpg = sig_cpgs_tot, all.cpg = all_tot, collection = "KEGG")
table(tot_gst.kegg$FDR < .05)
tot_top20_kegg<- topGSA(tot_gst.kegg, n=20)
write_xlsx(tot_top20_kegg, "Kegg_tot_top20_5x10-5.xlsx")

tri1_gst.kegg <- gometh(sig.cpg = sig_cpgs_tri1, all.cpg = all_tri1, collection = "KEGG")
table(tri1_gst.kegg$FDR < .05)
tri1_top20_kegg<- topGSA(tri1_gst.kegg, n=20)
write_xlsx(tri1_top20_kegg, "Kegg_tri1_top20_5x10-5.xlsx")

tri2_gst.kegg <- gometh(sig.cpg = sig_cpgs_tri2, all.cpg = all_tri2, collection = "KEGG")
table(tri2_gst.kegg$FDR < .05)
tri2_top20_kegg<- topGSA(tri2_gst.kegg, n=20)
write_xlsx(tri2_top20_kegg, "Kegg_tri2_top20_5x10-5.xlsx")

tri3_gst.kegg <- gometh(sig.cpg = sig_cpgs_tri3, all.cpg = all_tri3, collection = "KEGG")
table(tri3_gst.kegg$FDR < .05)
tri3_top20_kegg<- topGSA(tri3_gst.kegg, n=20)
write_xlsx(tri3_top20_kegg, "Kegg_tri3_top20_5x10-5.xlsx")

#\ END OF THIS SCRIPT
