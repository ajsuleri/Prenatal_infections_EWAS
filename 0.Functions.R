#########################################################################
#####################PROJECT: PRENATAL INFECTION EWAS####################
#########################################################################

#Author: Anna Suleri

# In this script all functions are noted that we use in the project. We will then source this script in the other scripts.
# Make sure right wd is loaded before saving results or data of interest 

########################################
#------------Baseline table------------#
########################################
baseline_table_function <- function(baselinevars, df){ #feed baselinevars vector and dataframe 
  
  #create empty dataframe to store results 
  results <- data.frame(Variable = character(0), Output = character(0))
  
  #for each baseline variable as specified in baselinevars decide whether continous or categorical and give output accordingly 
  for(i in baselinevars){ 
    
    #x = i vars that are columns in the dataframe df
    x <- df[, i] 
    
    #show column name as heading per output 
    message(i) 
    
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
    
    #if else to apply correct function for vars type 
    if (class(x) == 'numeric') {
      output <- summary_continuous(x)
    } 
    else 
    {
      output <- summary_categorical(x) 
    }
    results <- rbind(results, data.frame(Variable = i, Output = output))
  }
  write_xlsx(results, 'baseline_table_results.xlsx') #load write_xl package
}

#\ 

########################################
#------Summary of missing values-------#
########################################
miss_values_function <- function(df){
  missvalues <- cbind("# NA" = sort(colSums(is.na(df))),
                      "% NA" = round(sort(colMeans(is.na(df))) * 100, 2))
  print(missvalues)
}

#\ 

########################################
#---------Multiple imputation----------#
########################################
# first load mice library 
imputation_function <- function(data, exclude_imp_vars, exclude_predictors, method = "default") {
  
  # Running setup imputation run
  if (method == "rf") { #e.g. in case of brain to allow any type of combination between variables 
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
  saveRDS(imp.test, file = paste0(imputedDataName, ".rds"))
  
  # Output 
  return(list(imp0 = imp0, imp.test = imp.test))
}

#\ 

########################################
#------------BaCON ESTIMATES-----------#
########################################
# first load bacon library

performBaconAnalysis <- function(dataframe) {
  bc <- bacon(dataframe[["t value"]])
  
  # Display the results
  print(inflation(bc))
  print(bias(bc))
  
  # Create a histogram plot
  plot(bc, type = "hist")
  dev.off()
  
  # Create a QQ plot
  plot(bc, type = "qq")
  dev.off()
}

#\ 

########################################
#---------EWAS META-ANALYSIS-----------#
########################################
# Define a function to create a fixed-effets meta-analysis model for each dataset
meta_analyze <- function(data) {
  rma(yi = data$Estimate, sei = data$`Std. Error`, method = "FE")
}

#\ 
