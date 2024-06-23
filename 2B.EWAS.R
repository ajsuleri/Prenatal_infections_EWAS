#########################################################################
#####################PROJECT: PRENATAL INFECTION EWAS####################
#########################################################################

#Author: Anna Suleri

#This script is to repeat the EWAS but now without adjust for cell type covariates ast a post-hoc analysis 
#1. We create a dataframe of CpGs that overlap in 450k and in epic for both 450k and epic
#2. We run the following EWASES on each dataframe:
######a) total infection - CpG at birth

#3. We run these four EWASES for the following two models:
######a. Model 1: child sex, maternal age at delivery, maternal education, maternal smoking, parity, batch effects, cel type proportions
######b. Model 2: model 1 + additionally gestational age at birth and birth weight 
#4.We check the diagnostics of each EWAS with a Q-Q plot of the pvalues, lambda inflation factor   
#5.We check the quality of the ewas  + apply post processing
# Run this script on gena server 

###Set working directory and load packages###
library(foreign)
library(haven)
library(dplyr)
library(pbmcapply) 

setwd('set_path_to_data')

## source functions file
source('0.Functions_script.R')

## run this script on the GENA server 

###----------------STEP 1: Prepare dataframes----------------###
## load methylation overlap data, for 450k and then for epic
overlap_df_450k <- readRDS('overlap_df_450k.rds')
overlap_df_epic <- readRDS('overlap_df_epic.rds') 

## load in phenotypic data
df_phenotype_450k <- readRDS('df_final_450k.rds') 
df_phenotype_epic <- readRDS('df_final_epic.rds') 

#check str of variables and recode if necessary 
df_phenotype_450k$Sample_Plate <- as.factor(df_phenotype_450k$Sample_Plate)
df_phenotype_epic$Sample_Plate <- as.factor(df_phenotype_epic$Sample_Plate)

# correlation plot of infections + cell type proportions 
df_cor_450k <- dplyr::select(df_phenotype_450k, c('sumscore_inf_tot', 'Bcell', 'CD4T', 'CD8T', 'Gran', 'Mono', 'NK', 'nRBC'))
df_cor_epic <- dplyr::select(df_phenotype_epic, c('sumscore_inf_tot', 'Bcell', 'CD4T', 'CD8T', 'Gran', 'Mono', 'NK', 'nRBC'))

library(corrplot)
correlation <- cor(df_cor_450k, use="pairwise.complete.obs")

png(file='correlation_infections_celltypes.png', width = 300, height = 400)

corrplot::corrplot(
  correlation, 
  method = 'color',
  addCoef.col = "black", 
  number.cex = 0.6,
  type = 'lower', 
  diag = FALSE,  
  tl.col = 'black', 
  tl.cex = 0.7, 
  sig.level = 0.05,
  col = colorRampPalette(c("midnightblue", "white", "darkred"))(100),
  colnames = colnames(correlation)
)

dev.off()

#\

###----------------STEP 2: RUN EWAS----------------###
## functions to run ewas models for total infection score
#function for model 1, 450k
regress1 <- function(cpg){
  cpg_phenotype_data <- cbind(cpg, df_phenotype_450k)
  coef(summary(lm(cpg ~ sumscore_inf_tot_standardized + GENDER + SMOKE_ALL + AGE_M_v2 + PARITY + EDUCM_3groups + Sample_Plate, data = cpg_phenotype_data)))[2,]
}

#function for model 2, 450k
regress2 <- function(cpg){
  cpg_phenotype_data <- cbind(cpg, df_phenotype_450k)
  coef(summary(lm(cpg ~ sumscore_inf_tot_standardized + GENDER + SMOKE_ALL + AGE_M_v2 + PARITY + EDUCM_3groups + GESTBIR + Sample_Plate, data = cpg_phenotype_data)))[2,]
}

#function for model 1, epic 
regress3 <- function(cpg){
  cpg_phenotype_data <- cbind(cpg, df_phenotype_epic)
  coef(summary(lm(cpg ~ sumscore_inf_tot_standardized + GENDER + SMOKE_ALL + AGE_M_v2 + PARITY + EDUCM_3groups + Sample_Plate, data = cpg_phenotype_data)))[2,]
}

#function for model 2, epic
regress4 <- function(cpg){
  cpg_phenotype_data <- cbind(cpg, df_phenotype_epic)
  coef(summary(lm(cpg ~ sumscore_inf_tot_standardized + GENDER + SMOKE_ALL + AGE_M_v2 + PARITY + EDUCM_3groups + GESTBIR + Sample_Plate, data = cpg_phenotype_data)))[2,]
}


## run all the ewases and save results 
setwd('set_path_to_results_folder')

# first for 450k
#ewas 1: total infection - 450k data - model 1
overlap_ewas_450k_model1.list <- pbmclapply(overlap_df_450k, regress1, mc.cores = 12)
overlap_ewas_450k_model1.data <- as.data.frame(do.call(rbind, overlap_ewas_450k_model1.list))
overlap_ewas_450k_model1.data$cpg <- names(overlap_df_450k)
save(overlap_ewas_450k_model1.data, file = 'overlap_ewas_450k_model1_nocell.rdata')

#ewas 2: total infection - 450k data - model 2
overlap_ewas_450k_model2.list <- pbmclapply(overlap_df_450k, regress2, mc.cores = 12)
overlap_ewas_450k_model2.data <- as.data.frame(do.call(rbind, overlap_ewas_450k_model2.list))
overlap_ewas_450k_model2.data$cpg <- names(overlap_df_450k)
save(overlap_ewas_450k_model2.data, file = 'overlap_ewas_450k_model2_nocell.rdata')

# then for epic

#ewas 9: total infection - epic data - model 1
overlap_ewas_epic_model1.list <- pbmclapply(overlap_df_epic, regress3, mc.cores = 12)
overlap_ewas_epic_model1.data <- as.data.frame(do.call(rbind, overlap_ewas_epic_model1.list))
overlap_ewas_epic_model1.data$cpg <- names(overlap_df_epic)
save(overlap_ewas_epic_model1.data, file = 'overlap_ewas_epic_model1_nocell.rdata')

#ewas 10: total infection - epic data- model 2
overlap_ewas_epic_model2.list <- pbmclapply(overlap_df_epic, regress4, mc.cores = 12)
overlap_ewas_epic_model2.data <- as.data.frame(do.call(rbind, overlap_ewas_epic_model2.list))
overlap_ewas_epic_model2.data$cpg <- names(overlap_df_epic)
save(overlap_ewas_epic_model2.data, file = 'overlap_ewas_epic_model2_nocell.rdata')

#\ 

###----------------STEP 4: Metafor meta-analyses----------------###

# load library for meta analysis and to save results
library(metafor)

# we will run an inverse weighted fixed-effects meta-analysis for total infection, trimester 1 infection, trimester 2 infection and trimester 3 infection. 
# For each exposure type we have a model 1 and 2. 
# We will meta-analyze the output of 16 ewases of the 450k and epic array together
# This will result in a final of 8 ewases
# because the vectors in each dataframes are too big, we need to do this with parallel processing 

setwd('set_path_to_meta_analysis_results_folder')

# combine datasets for each model and exposure together that we want to meta-analyze together 
total_infection_model1 <- rbind(overlap_ewas_450k_model1.data, overlap_ewas_epic_model1.data)
total_infection_model2 <- rbind(overlap_ewas_450k_model2.data, overlap_ewas_epic_model2.data)

# Use pbmclapply to run the models in parallel
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

# Obtain the summary of the meta-analysis results
result_total_infection_model1_df <- as.data.frame(matrix(unlist(results_total_infection_model1), ncol = 4, byrow = TRUE))
colnames(result_total_infection_model1_df) <- c("cpg","meta_beta", "meta_se", "meta_pval")
result_total_infection_model1_df_unique <- subset(result_total_infection_model1_df, !duplicated(cpg))

result_total_infection_model2_df <- as.data.frame(matrix(unlist(results_total_infection_model2), ncol = 4, byrow = TRUE))
colnames(result_total_infection_model2_df) <- c("cpg","meta_beta", "meta_se", "meta_pval")
result_total_infection_model2_df_unique <- subset(result_total_infection_model2_df, !duplicated(cpg))

# Save the results to an Excel file
write.csv(result_total_infection_model1_df_unique, "meta_analyzed_results_total_infection_model1_nocell.csv", row.names = FALSE)
write.csv(result_total_infection_model2_df_unique, "meta_analyzed_results_total_infection_model2_nocell.csv", row.names = FALSE)

#\ 

###----------------STEP 5: Quality checks meta-analysis EWAS----------------###
meta_analyzed_results_total_infection_model1_nocell <- read.csv("meta_analyzed_results_total_infection_model1_nocell.csv")
meta_analyzed_results_total_infection_model2_nocell <- read.csv("meta_analyzed_results_total_infection_model2_nocell.csv")

# calculate lambda for meta analyzed results
Lambda<-function(P){
  chisq <- qchisq(1-P,1)
  median(chisq,na.rm=T)/qchisq(0.5,1)
}

Lambda(meta_analyzed_results_total_infection_model1_nocell$meta_pval) #0.883
Lambda(meta_analyzed_results_total_infection_model2_nocell$meta_pval) #0.883

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
 
png(file='qq_plot_revisions_ewas_nocell_model1.png')
qqunif.plot(meta_analyzed_results_total_infection_model1_nocell$meta_pval)
dev.off()

png(file='qq_plot_revisions_ewas_nocell_model2.png')
qqunif.plot(meta_analyzed_results_total_infection_model2_nocell$meta_pval)
dev.off()

#\

###----------------STEP 6: Annotation, creating basic tables, postprocessing----------------###
# set wd 
setwd('path_to_cleaned_results_folder')

## now we will annotate the files 
library(IlluminaHumanMethylation450kanno.ilmn12.hg19) 

# first annotate dataframes with chromosome and position 
annotation <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
annotation$cpg <- annotation$Name # get consistent names 

# merge annotation with df s
meta_analyzed_results_total_infection_model1_nocell_annotated <- merge(meta_analyzed_results_total_infection_model1_nocell, annotation, by = 'cpg', all.x = TRUE)
meta_analyzed_results_total_infection_model2_nocell_annotated <- merge(meta_analyzed_results_total_infection_model2_nocell, annotation, by = 'cpg', all.x = TRUE)

write.csv(meta_analyzed_results_total_infection_model1_nocell_annotated, 'meta_analyzed_results_total_infection_model1_nocell_annotated.csv')
write.csv(meta_analyzed_results_total_infection_model2_nocell_annotated, 'meta_analyzed_results_total_infection_model2_nocell_annotated.csv')

# first make tables with only p < 0.05 (nominal significance)
sign_meta_analyzed_results_total_infection_model1_nocell_annotated <- meta_analyzed_results_total_infection_model1_nocell_annotated[meta_analyzed_results_total_infection_model1_nocell_annotated$meta_pval < 0.05,]
sign_meta_analyzed_results_total_infection_model2_nocell_annotated <- meta_analyzed_results_total_infection_model2_nocell_annotated[meta_analyzed_results_total_infection_model2_nocell_annotated$meta_pval < 0.05,]

## Post processing 
library(writexl)

## read crossreactive & flagged probes: combination of Naeem and Chen 
cross_reactive_probes<-read.csv("crossreactiveprobes.csv")
flagged_probes<-read.csv("flaggedprobes.csv")

## exclude cross reactive probes 
cleaned_sign_results_meta_analyzed_results_total_infection_model1_nocell <- subset(sign_meta_analyzed_results_total_infection_model1_nocell_annotated, !(sign_meta_analyzed_results_total_infection_model1_nocell_annotated$cpg %in% cross_reactive_probes$MarkerName))
cleaned_sign_results_meta_analyzed_results_total_infection_model2_nocell <- subset(sign_meta_analyzed_results_total_infection_model2_nocell_annotated, !(sign_meta_analyzed_results_total_infection_model2_nocell_annotated$cpg %in% cross_reactive_probes$MarkerName))

## flag probes containing snps
# first align name with cpg name of our meta analysis results file 
colnames(flagged_probes)[1] <- "cpg"

final_cleaned_sign_results_meta_analyzed_results_total_infection_model1_nocell <- merge(cleaned_sign_results_meta_analyzed_results_total_infection_model1_nocell, flagged_probes, by.x = "cpg", by.y = "cpg", all.x = T, all.y = F)
final_cleaned_sign_results_meta_analyzed_results_total_infection_model2_nocell <- merge(cleaned_sign_results_meta_analyzed_results_total_infection_model2_nocell, flagged_probes, by.x = "cpg", by.y = "cpg", all.x = T, all.y = F)

## save new datasets
final_cleaned_sign_results_meta_analyzed_results_total_infection_model1_nocell <- as.data.frame(final_cleaned_sign_results_meta_analyzed_results_total_infection_model1_nocell)
final_cleaned_sign_results_meta_analyzed_results_total_infection_model2_nocell <- as.data.frame(final_cleaned_sign_results_meta_analyzed_results_total_infection_model2_nocell)

write_xlsx(final_cleaned_sign_results_meta_analyzed_results_total_infection_model1_nocell, 'final_cleaned_sign_results_meta_analyzed_results_total_infection_model1_nocell.xlsx')
write_xlsx(final_cleaned_sign_results_meta_analyzed_results_total_infection_model2_nocell, 'final_cleaned_sign_results_meta_analyzed_results_total_infection_model2_nocell.xlsx')

#\ END OF THIS SCRIPT 
