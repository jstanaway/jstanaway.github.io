# Project: iNTS Forecasting
#
# Purpose: 
# This script runs crossvalidation on the meta-model to determine the combination of hyperparameter values that 
# result in the smallest out-of-sample prediction errors.  There are two relevant hyperparameters: alpha, which
# governs the mix of L1 and L2 regularization, and lambda, which governs the strength of regularization.  This
# script takes as input a single value each form alpha and lambda, runs the meta-model with those hyperparameter
# values on each fold of the crossvalidation, and calculates the out-of-sample prediction errors.  The results
# are saved to disk and will be used in the next step to select the best models to run on the full dataset.
#
# Author: Jeff Stanaway (stanaway@uw.edu)
# Date created: February 2024
# Date last modified: 25 March 2024
#
# Libraries: data.table, feather, glmnet, parallel, pracma
# User inputs: ROOT_DIR -- everything else is passed as arguments from the launch script

# Clear the environment
rm(list = ls())

# Load libraries
library(data.table)
library(feather)
library(glmnet)
library(parallel)
library(pracma)


# Determine if this is an interactive session (i.e. development or debugging) or not
INTERACTIVE = commandArgs()[2]=='--interactive'

print(commandArgs())

# If not interactive, get run arguments from the command line
if (INTERACTIVE == FALSE) {  
  args <- commandArgs(trailingOnly = T)
  alpha <- as.integer(args[1])
  lambda <- as.integer(args[2])
  xval_version <- args[3]
  
# If interactive hardcode testing parameters  
} else { 
  alpha <- 50
  lambda <- -1
  xval_version <- 3
}

# Set directories
ROOT_DIR <- 'FILEPATH_REMOVED_FOR_SECURITY'
xval_dir <- file.path(ROOT_DIR, 'xvalidation', xval_version)

# Load file with x, y, pf, train_indices, and test_indices (created by launch script)
load(file.path(xval_dir, 'xval_glmnet_data.RData'))


# Define a function to run the model and calculate fit statistics
xval_glmnet <- function(alpha, lambda, fold, std) {
  # If fold is 999, then we are running on the full dataset
  if (fold == 999) {
    x_train <- x
    y_train <- y
    x_test <- x
    y_test <- y
    
    # Otherwise, we are running on a fold of the crossvalidation
  } else {
    x_train <- x[train_indices[[fold]], ]
    y_train <- y[train_indices[[fold]]]
    x_test <- x[test_indices[[fold]], ]
    y_test <- y[test_indices[[fold]]]
  }
  
  # Run the model on the training data
  mdl <- glmnet(x_train, y_train, alpha = alpha, lambda = lambda, 
                penalty.factor = pf, standardize = std)
  
  # Predict on the test data
  pred <- data.table(truth = y_test, pred = predict(mdl, newx = x_test, s = "lambda.min"))
  
  # Get the included predictors
  mdl_coefs <- coef(mdl)
  weak_learners  <- mdl_coefs@Dimnames[[1]][-1][mdl_coefs@i]
  
  # Calculate fit statistics OOS
  fit_stats <- data.table(t(unlist(rmserr(pred$truth, pred$pred))))
  fit_stats[, `:=`(alpha = alpha, lambda = lambda, fold = fold, std = std, included = list(weak_learners))]
  return(fit_stats)
}


# Run the model on each fold of the crossvalidation
# Note: that alpha and lambda are passed as integers to avoid rounding errors and simplify
# the naming of the output files.  The math on alpha and lambda below are to convert those
# integers to the actual desired alpha and lambda values

fit_stats <- rbindlist(mclapply(c(1:length(train_indices), 999), function(f) {
  rbindlist(lapply(c(T, F), function(std) xval_glmnet(alpha/100, 10^(lambda/10), f, std)))}, mc.cores = 7))




# Save the resulting fit statistics (need to use RDS format because of the list column)
saveRDS(fit_stats, file.path(xval_dir, paste0('step_2_fit_stats_alpha', alpha, '_lambda', lambda, '.rds')))
