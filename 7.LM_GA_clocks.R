# Author: Kristina Salontaji, k.salontaji@erasmusmc.nl 

# create an output folder that we will save all our data to within the working directory we just specified
dir.create("output_folder") 

# change to your own working directory 
setwd("path_to_data") 

data - readRDS(file = "data.rds")

library(dplyr)
library(mice)
library(broom.mixed)
library(lattice)
library(MASS)
library(robustbase)

master_clock <- data

clocks <- c("EPIC_ffk", "Bohlin_ffk", "ageAcc2.Bohlin_ffk", "ageAcc2.EPIC_ffk")

inc <- lmrob.control() #increase number of iterations as lmrob may fail to converge in rare cases
inc$maxit.scale<-2500 #you may increase this number if a model fails to converge (you will get a warning saying that it fails to converge)
inc$k.max <-2500 #same here
set.seed(2023)

#using non-imputed data
for (j in clocks){
  clock = j
  print(clock)
  model0 <- with(master_clock, lmrob(scale(get(clock)) ~ scale(sumscore_inf_tot) + sex + agemother_birthchild + education_mother + smoke_mother + PARITY + bcell_ffk + cd4t_ffk + cd8t_ffk + gran_ffk + mono_ffk + nk_ffk + nRBC_ffk + Sample_Plate_ffk, control = inc))
  summary <- tidy(model0, conf.int = TRUE, conf.level = 0.95)
  write.csv(x = summary, file = paste0("./annaproj_output/", clock, "_mod1_tot.csv"), row.names = FALSE) # save to our output folder
  print("finished mod1")
  
  model0 <- with(master_clock, lmrob(scale(get(clock)) ~ scale(sumscore_inf_tot) + sex + agemother_birthchild + education_mother + smoke_mother + PARITY + bcell_ffk + cd4t_ffk + cd8t_ffk + gran_ffk + mono_ffk + nk_ffk + nRBC_ffk + Sample_Plate_ffk + GESTBIR, control = inc))
  summary <- tidy(model0, conf.int = TRUE, conf.level = 0.95)
  write.csv(x = summary, file = paste0("./annaproj_output/", clock, "_mod2_tot.csv"), row.names = FALSE) # save to our output folder
  print("finished mod2")
  
  
  #trimester 1 specific
  model0 <- with(master_clock, lmrob(scale(get(clock)) ~ scale(sumscore_inf_tri1) + sex + agemother_birthchild + education_mother + smoke_mother + PARITY + bcell_ffk + cd4t_ffk + cd8t_ffk + gran_ffk + mono_ffk + nk_ffk + nRBC_ffk + Sample_Plate_ffk, control = inc))
  summary <- tidy(model0, conf.int = TRUE, conf.level = 0.95)
  write.csv(x = summary, file = paste0("./annaproj_output/", clock, "_mod1_tri1.csv"), row.names = FALSE) # save to our output folder
  print("finished mod1")
  
  model0 <- with(master_clock, lmrob(scale(get(clock)) ~ scale(sumscore_inf_tri1) + sex + agemother_birthchild + education_mother + smoke_mother + PARITY + bcell_ffk + cd4t_ffk + cd8t_ffk + gran_ffk + mono_ffk + nk_ffk + nRBC_ffk + Sample_Plate_ffk + GESTBIR, control = inc))
  summary <- tidy(model0, conf.int = TRUE, conf.level = 0.95)
  write.csv(x = summary, file = paste0("./annaproj_output/", clock, "_mod2_tri1.csv"), row.names = FALSE) # save to our output folder
  print("finished mod2")
  
  
  #trimester 2 specific
  model0 <- with(master_clock, lmrob(scale(get(clock)) ~ scale(sumscore_inf_tri2) + sex + agemother_birthchild + education_mother + smoke_mother + PARITY + bcell_ffk + cd4t_ffk + cd8t_ffk + gran_ffk + mono_ffk + nk_ffk + nRBC_ffk + Sample_Plate_ffk, control = inc))
  summary <- tidy(model0, conf.int = TRUE, conf.level = 0.95)
  write.csv(x = summary, file = paste0("./annaproj_output/", clock, "_mod1_tri2.csv"), row.names = FALSE) # save to our output folder
  print("finished mod1")
  
  model0 <- with(master_clock, lmrob(scale(get(clock)) ~ scale(sumscore_inf_tri2) + sex + agemother_birthchild + education_mother + smoke_mother + PARITY + bcell_ffk + cd4t_ffk + cd8t_ffk + gran_ffk + mono_ffk + nk_ffk + nRBC_ffk + Sample_Plate_ffk + GESTBIR, control = inc))
  summary <- tidy(model0, conf.int = TRUE, conf.level = 0.95)
  write.csv(x = summary, file = paste0("./annaproj_output/", clock, "_mod2_tri2.csv"), row.names = FALSE) # save to our output folder
  print("finished mod2")
  
  #trimester 3 specific
  model0 <- with(master_clock, lmrob(scale(get(clock)) ~ scale(sumscore_inf_tri3) + sex + agemother_birthchild + education_mother + smoke_mother + PARITY + bcell_ffk + cd4t_ffk + cd8t_ffk + gran_ffk + mono_ffk + nk_ffk + nRBC_ffk + Sample_Plate_ffk, control = inc))
  summary <- tidy(model0, conf.int = TRUE, conf.level = 0.95)
  write.csv(x = summary, file = paste0("./annaproj_output/", clock, "_mod1_tri3.csv"), row.names = FALSE) # save to our output folder
  print("finished mod1")
  
  model0 <- with(master_clock, lmrob(scale(get(clock)) ~ scale(sumscore_inf_tri3) + sex + agemother_birthchild + education_mother + smoke_mother + PARITY + bcell_ffk + cd4t_ffk + cd8t_ffk + gran_ffk + mono_ffk + nk_ffk + nRBC_ffk + Sample_Plate_ffk + GESTBIR, control = inc))
  summary <- tidy(model0, conf.int = TRUE, conf.level = 0.95)
  write.csv(x = summary, file = paste0("./annaproj_output/", clock, "_mod2_tri3.csv"), row.names = FALSE) # save to our output folder
  print("finished mod2")
  }



master_clock <- data

clocks <- c("ageAcc2.Bohlin_epic", "ageAcc2.EPIC_epic", "EPIC_epic", "Bohlin_epic")


inc <- lmrob.control() #increase number of iterations as lmrob may fail to converge in rare cases
inc$maxit.scale<-2500 #you may increase this number if a model fails to converge (you will get a warning saying that it fails to converge)
inc$k.max <-2500 #same here
set.seed(2023)

#using non-imputed data
for (j in clocks){
  clock = j
  print(clock)
  model0 <- with(master_clock, lmrob(scale(get(clock)) ~ scale(sumscore_inf_tot) + sex + agemother_birthchild + education_mother + smoke_mother + PARITY + bcell_epic + cd4t_epic + cd8t_epic + gran_epic + mono_epic + nk_epic + nRBC_epic + Sample_Plate_epic, control = inc))
  summary <- tidy(model0, conf.int = TRUE, conf.level = 0.95)
  write.csv(x = summary, file = paste0("./annaproj_output/", clock, "_mod1_tot.csv"), row.names = FALSE) # save to our output folder
  print("finished mod1")
  
  model0 <- with(master_clock, lmrob(scale(get(clock)) ~ scale(sumscore_inf_tot) + sex + agemother_birthchild + education_mother + smoke_mother + PARITY + bcell_epic + cd4t_epic + cd8t_epic + gran_epic + mono_epic + nk_epic + nRBC_epic + Sample_Plate_epic + GESTBIR.y, control = inc))
  summary <- tidy(model0, conf.int = TRUE, conf.level = 0.95)
  write.csv(x = summary, file = paste0("./annaproj_output/", clock, "_mod2_tot.csv"), row.names = FALSE) # save to our output folder
  print("finished mod2")
  
  
  #trimester 1 specific
  model0 <- with(master_clock, lmrob(scale(get(clock)) ~ scale(sumscore_inf_tri1) + sex + agemother_birthchild + education_mother + smoke_mother + PARITY + bcell_epic + cd4t_epic + cd8t_epic + gran_epic + mono_epic + nk_epic + nRBC_epic + Sample_Plate_epic, control = inc))
  summary <- tidy(model0, conf.int = TRUE, conf.level = 0.95)
  write.csv(x = summary, file = paste0("./annaproj_output/", clock, "_mod1_tri1.csv"), row.names = FALSE) # save to our output folder
  print("finished mod1")
  
  model0 <- with(master_clock, lmrob(scale(get(clock)) ~ scale(sumscore_inf_tri1) + sex + agemother_birthchild + education_mother + smoke_mother + PARITY + bcell_epic + cd4t_epic + cd8t_epic + gran_epic + mono_epic + nk_epic + nRBC_epic + Sample_Plate_epic + GESTBIR.y, control = inc))
  summary <- tidy(model0, conf.int = TRUE, conf.level = 0.95)
  write.csv(x = summary, file = paste0("./annaproj_output/", clock, "_mod2_tri1.csv"), row.names = FALSE) # save to our output folder
  print("finished mod2")
  
  
  #trimester 2 specific
  model0 <- with(master_clock, lmrob(scale(get(clock)) ~ scale(sumscore_inf_tri2) + sex + agemother_birthchild + education_mother + smoke_mother + PARITY + bcell_epic + cd4t_epic + cd8t_epic + gran_epic + mono_epic + nk_epic + nRBC_epic + Sample_Plate_epic, control = inc))
  summary <- tidy(model0, conf.int = TRUE, conf.level = 0.95)
  write.csv(x = summary, file = paste0("./annaproj_output/", clock, "_mod1_tri2.csv"), row.names = FALSE) # save to our output folder
  print("finished mod1")
  
  model0 <- with(master_clock, lmrob(scale(get(clock)) ~ scale(sumscore_inf_tri2) + sex + agemother_birthchild + education_mother + smoke_mother + PARITY + bcell_epic + cd4t_epic + cd8t_epic + gran_epic + mono_epic + nk_epic + nRBC_epic + Sample_Plate_epic + GESTBIR.y, control = inc))
  summary <- tidy(model0, conf.int = TRUE, conf.level = 0.95)
  write.csv(x = summary, file = paste0("./annaproj_output/", clock, "_mod2_tri2.csv"), row.names = FALSE) # save to our output folder
  print("finished mod2")
  
  #trimester 3 specific
  model0 <- with(master_clock, lmrob(scale(get(clock)) ~ scale(sumscore_inf_tri3) + sex + agemother_birthchild + education_mother + smoke_mother + PARITY + bcell_epic + cd4t_epic + cd8t_epic + gran_epic + mono_epic + nk_epic + nRBC_epic + Sample_Plate_epic, control = inc))
  summary <- tidy(model0, conf.int = TRUE, conf.level = 0.95)
  write.csv(x = summary, file = paste0("./annaproj_output/", clock, "_mod1_tri3.csv"), row.names = FALSE) # save to our output folder
  print("finished mod1")
  
  model0 <- with(master_clock, lmrob(scale(get(clock)) ~ scale(sumscore_inf_tri3) + sex + agemother_birthchild + education_mother + smoke_mother + PARITY + bcell_epic + cd4t_epic + cd8t_epic + gran_epic + mono_epic + nk_epic + nRBC_epic + Sample_Plate_epic + GESTBIR.y, control = inc))
  summary <- tidy(model0, conf.int = TRUE, conf.level = 0.95)
  write.csv(x = summary, file = paste0("./annaproj_output/", clock, "_mod2_tri3.csv"), row.names = FALSE) # save to our output folder
  print("finished mod2")
}

#\ END OF INDIVIDUAL REGRESSIONS

#meta analysis

# https://bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/pooling-es.html
# Linear Regression — Linear Regression - F-Squared. The effect size measure of choice for (simple and multiple) linear regression is f2.
# What effect size do you use for regression?
# Cohen's ƒ2 is a measure of effect size used for a (multiple) regression. Effect size measures for ƒ2are 0.02, 0.15, and 0.35, 
# indicating small, medium, and large, respectively.

#load packages
library(metafor)
library(readr)
library(purrr)
library(stringr)
library(meta)
library(dmetar)

#load the datasets
table(is.na(master_clock$sumscore_inf_tot),is.na(master_clock$Bohlin_epic))

N_ffk = 1106
N_epic = 871

dir_genr = "annaproj_output/ffk/"
filenames <- list.files(dir_genr, pattern="*.csv", full.names=TRUE)
filenames_short <- list.files(dir_genr, pattern="*.csv", full.names=FALSE)
modelnames <- gsub("\\.csv*$","",filenames_short) 
ldf <- lapply(filenames, read.csv)
names(ldf) <- modelnames

ffk <- map_df(ldf, ~as.data.frame(.x), .id="id")
ffk$N <- N_ffk
ffk$author <- "ffk"

dir_genr = "annaproj_output/epic/"
filenames <- list.files(dir_genr, pattern="*.csv", full.names=TRUE)
filenames_short <- list.files(dir_genr, pattern="*.csv", full.names=FALSE)
modelnames <- gsub("\\.csv*$","",filenames_short) 
ldf <- lapply(filenames, read.csv)
names(ldf) <- modelnames


epic <- map_df(ldf, ~as.data.frame(.x), .id="id")
epic$N <- N_epic
epic$author <- "epic"


#append all to one data frame for all cohorts
data_metaanalysis = rbind(ffk, epic)
data_metaanalysis$low <- data_metaanalysis$conf.low
data_metaanalysis$hi <- data_metaanalysis$conf.high
data_metaanalysis$est <- data_metaanalysis$estimate


#keep only row names that are a clock variable - partial string match ageA

#data_metaanalysis <- data_metaanalysis[str_detect(data$term, "ageA"), ]

#data$SD_val= sqrt(data$N)*data$estimate

#data$cohens_d = data$estimate/data$SD_val #no, I cannot use it like this!!!!!!!!!!

#data$weight = 1/data$std.error^2

saveRDS(data_metaanalysis, file = "data_metaanalysis_annaproj.rds")

#\ END OF THIS SCRIPT.
