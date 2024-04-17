# Project: iNTS Forecasting
#
# Purpose: 
# This script runs Shape Constrained Additive Models (SCAMs) to predict iNTS incidence. These models are the
# weak-learners in the ensemble model. It is called at three points in the analysis: 1) it is first called by
# 02c_xval_model_run_launch.R to run fold-specific models and assess out-of-sample prediction erros during 
# the crossvalidation, 2) it is called by 02d_xval_select_best_candidates.R to run the best candidate models
# on the full dataset (selected based on OOS performance in the previous step), and 3) it is called by
# 02e_xval_ar1_model_launch.R to determine the residual autocorrelation in the full models, calculate the
# the autocorrelation coefficient, and then run the model on the full dataset with AR1 errors to account 
# for residual autocorrelation.  Note that this extra step is required because, while SCAM can account for 
# autocorrelation, it cannot estimate the autocorrelation coefficient and requires the user to pass that
# in as an argument. 
#
# Author: Jeff Stanaway (stanaway@uw.edu)
# Date created: April 2023
# Date last modified: 25 March 2024
#
# Libraries: data.table, scam, pracma, DescTools, pcse, prais
# User inputs: No hard coded inputs -- everything is passed as arguments from the launch script


rm(list = ls())

library(data.table)
library(scam)
library(pracma)
library(DescTools)
library(pcse)
library(prais)

# Determine if this is an interactive run (i.e. for development and debugging) or not
INTERACTIVE = commandArgs()[2]=='--interactive'

print(commandArgs())

if (!INTERACTIVE) {  
  args <- commandArgs(trailingOnly = T)
  model_id <- args[1]
  fold_id <- args[2]
  step <- args[3]
  xval_version <- args[4]
  is_ar1 <- as.logical(args[5])
  
} else { # running interactively (need to hard code inputs for testing)
  model_id <- 1127
  fold_id <- 999
  step <- 1
  xval_version <- 3
  is_ar1 <- TRUE
}


message(paste0("Running model ", model_id, " on fold ", fold_id, " for step ", step, " with AR1 = ", is_ar1))

root_dir <- 'FILEPATH_REMOVED_FOR_SECURITY'
xval_dir <- file.path(root_dir, 'xvalidation', xval_version)



### READ IN INPUTS ##

# Read in the cross-validation folds #
fold <- fread(file.path(root_dir, 'inputs', 'xvalidation_folds.csv'))
ref_year <- max(fold$year_id)

# Read in the model and cross-validation fold information #
message('Reading in the model information')
model <- readRDS(file.path(xval_dir, paste0('step', step, '_models.rds')))[model_id, ]

# Read in the input dataset (No iNTS in age_group 2, so dropping those rows) #
message('Reading in the dataset')
all_points <- fread(file.path(root_dir, "FILEPATH_REMOVED_FOR_SECURITY"))[age_group_id > 2 & year_id <= ref_year, ]
all_points[, age_group_id := as.factor(age_group_id)]
all_points[, ihme_loc_id:= as.factor(ihme_loc_id)]




### ESTABLISH TRAIN AND TEST SETS ###

# Split data into training and testing -- fold_id 999 indicates a run on the full data set #
message('Creating train and test samples')
if (fold_id == 999) {
  train <- copy(all_points)
  test <- copy(all_points)
} else {
  fold <- fold[fold == fold_id, ]
  
  train <- copy(all_points)[year_id %in% fold[set == 'train', year_id]]
  test  <- copy(all_points)[year_id %in% fold[set == 'test', year_id]]
}



### RUN MODEL ###

if (!is_ar1) {
  # Run the model on the training data and save the model object     
  message('Running the SCAM model')
  mdl <- scam(as.formula(model), data = train)
  
  message('SCAM model complete.  Saving model object')
  saveRDS(mdl, file.path(xval_dir, paste0('step', step, '_model', model_id, '_fold', fold_id, '.rds')))
  
} else {
  # We need to estimate the AR1 coefficient, accounting for predictors.  
  # Since our model is non-linear, and we have panel data we can't do this with a simple ARIMA 
  # Instead, we'll use predictions from the naive model (i.e. the model not accounting for autocorrelation)
  # as the external regressor, and estimate the Prais-Winsten Estimator for AR1 with a panel data structure
  
  # Bring in the naive model and predict the outcome
  naive_mdl <- readRDS(file.path(xval_dir, paste0('step', step, '_model', model_id, '_fold', fold_id, '.rds')))
  all_points[, pred_naive := as.numeric(predict(naive_mdl, all_points))]                     

  # Estimate the AR1 coefficient using the Prais-Winsten Estimator
  message('Estimating the AR1 coefficient')
  ar_coef_mdl <- prais_winsten(iNTS_log ~ pred_naive, data = all_points, index = c('location_id', 'age_group_id', 'sex_id', 'year_id'))
  rho <- rev(ar_coef_mdl$rho)[1] # The Prais-Winsten Estimator uses iteration, and we want the last estimates, so reversing it 
  
  # With panel data, scam needs a logical variable to indicate the starting observation for each panel.  
  # Creating the panel start variable here
  train[, ar_start := (year_id == min(year_id)), by = c('location_id', 'age_group_id', 'sex_id')]
  train <- train[order(location_id, sex_id, age_group_id, year_id), ]

  # Run the model on the training data and save the model object #      
  message('Running the SCAM model')
  mdl <- scam(as.formula(model), data = train,  AR1.rho = rho, AR.start = train$ar_start)
  
  message('SCAM model complete.  Saving model object')
  saveRDS(mdl, file.path(xval_dir, paste0('step', step, '_model', model_id, '_fold', fold_id, '_ar1.rds')))
}



### CALCULATE FIT STATISTICS ###

# If we're running a fold for cross-validation, we need to predict on the test data and calculate fit statistics #
if (fold_id != 999) {
  message('Predicting and estimating OOS fit statistics')
  # Use model to predict in the test data #
  test$pred <- as.numeric(predict(mdl, newdata = test))
  test <- test[order(location_id, age_group_id, sex_id, year_id), ]
  
  # Calculate annual change in observed data and predictions -- need to asses OOS trend prediction accuracy
  test[, paste0(c('iNTS_log', 'pred'), '_trend') := lapply(.SD, function(x) x - shift(x)), by = c('location_id', 'age_group_id', 'sex_id'), .SDcols = c('iNTS_log', 'pred')]
  
  # Determine trend directions in observed data and predictions, and make confusion matrix (for calculating Kappa)
  test[, dir_truth := factor(sign(iNTS_log_trend), levels = -1:1, labels = c('neg', 'zero', 'pos'))]
  test[, dir_pred := factor(sign(pred_trend), levels = -1:1, labels = c('neg', 'zero', 'pos'))]
  cmat <- table(test$dir_pred, test$dir_truth)
  
  # Calculate measures of OOS fit accuracy and place in a neat data.table for export
  fit_stats_est <- data.table(t(unlist(rmserr(test$iNTS_log, test$pred))))
  fit_stats_trend <- data.table(t(unlist(rmserr(na.omit(test$iNTS_log_trend), na.omit(test$pred_trend)))))
  names(fit_stats_trend) <- paste0(names(fit_stats_trend), '_trend')
  
  # Bring all of the fit stats and model information together in preparation for export
  fit_stats <- cbind(fit_stats_est, fit_stats_trend, 
                     pr_correct_sign = mean(na.omit(test$dir_truth) == na.omit(test$dir_pred)), 
                     sign_kappa = CohenKappa(cmat, weights = "Equal-Spacing"),
                     model = model)
  
  # Save the fit stats
  write.csv(fit_stats, file.path(xval_dir, paste0('fit_stats_step', step, '_model', model_id, '_fold', fold_id, '.csv')), row.names = F)
}

message('Done.')