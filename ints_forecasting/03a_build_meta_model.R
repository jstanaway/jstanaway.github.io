# Project: iNTS Forecasting
#
# Purpose: 
# This script builds the meta-model based on weak learner predictions. The meta-model here is an elastic net
# with a penalty factor derived from OOS RMSE, such that the models with smaller OOS errors are given lower
# penalties (and therefore more likely to be included and with larger coefficients in the meta-model), and 
# models with larger OOS errors are given larger penalties (and therefore less likely to be included or have
# smaller coefficients in the meta-model.  Note that only those weak-learners meeting the 1SE threshold in 
# the previous step are eligible for inclusion in the meta model. It also launches crossvalidation on the 
# meta-model to determine the combination of hyperparameter values that result in the smallest out-of-sample 
# prediction errors.  There are two relevant hyperparameters: alpha, which governs the mix of L1 and L2 
# regularization (i.e. the degree to which the model behaves like lasso vs ridge regression), and lambda, 
# which governs the strength of regularization (i.e. the amount of shrinkage).
#
# The overall structure of the script is...
# 1) Read in the cross_validation fit statistics and calculate the penalty factor
# 2) Create some diagnostic plots to evaluate the penalty factor and 1SE threshold relative to RMSE
# 3) Read in predictions from all eligible candidate weak-learners, and create needed objects to run 
#    crossvalidation to select alpha and lambda parameters for the elastic net meta model
# 4) Launch the jobs to run the crossvalidation models (note, this could be done with a pre-build package
#    package like caret, but the run time and memory requirements were excessive with that approach, and
#    the parallelized approach was much faster)
# 5) Read in the crossvalidation results and select and fit the best model
#
# Unlike the previous crossvalidation steps, this one runs pretty quickly: if the cluster is reasonably 
# clear this process can run in well under an hour
#
# Author: Jeff Stanaway (stanaway@uw.edu)
# Date created: February 2024
# Date last modified: 25 March 2024
#
# Libraries: data.table, tidyverse, viridis, pbapply
# User inputs: RELEASE, XVAL_VERSION, ROOT_DIR, THREADS, MEM (you'll likely not need to change the last two)


# Clear the environment
rm(list = ls())

# Load libraries
library(data.table)
library(tidyverse)
library(viridis)
library(pbapply)

# Set user inputs
RELEASE <- 6
XVAL_VERSION <- 3
ROOT_DIR <- 'FILEPATH_REMOVED_FOR_SECURITY'
user <- Sys.getenv('USER')


# Build out filepaths
xval_dir <- file.path(ROOT_DIR, 'xvalidation', XVAL_VERSION)
graph_dir <- file.path(xval_dir, 'graphs')

dir.create(graph_dir, recursive = T, showWarnings = F)


# Read in the compiled fit statistics
complete <- readRDS(file.path(xval_dir, 'step1_fit_stats_complete.rds'))


# Check for the existence of the full and AR1 model fits -- we exclude any model that failed to converge
complete[, base_status := file.exists(file.path(xval_dir, paste0('step', step, '_model', model_id, '_fold999.rds')))]
complete[, ar1_status := file.exists(file.path(xval_dir, paste0('step', step, '_model', model_id, '_fold999_ar1.rds')))]
complete[, ar1_pred_file := file.path(xval_dir, paste0('step', step, '_model', model_id, '_fold999_ar1_preds.feather'))]
complete[, ar1_preds_status := file.exists(ar1_pred_file)]


# Scale RMSE and trend RMSE to 0-1 to make comparable
complete[, rmse_mean_scaled := (rmse_mean - min(rmse_mean)) / (max(rmse_mean) - min(rmse_mean))]
complete[, rmse_trend_mean_scaled := (rmse_trend_mean - min(rmse_trend_mean)) / (max(rmse_trend_mean) - min(rmse_trend_mean))]
complete[, rmse_total := rmse_mean_scaled + rmse_trend_mean_scaled]
complete[, rmse_std_total := sqrt(((mse_mean / mean(mse_mean)) + (mse_trend_mean / mean(mse_trend_mean)))/2)]


# Create a penalty factor as the total scaled RMSE minus the minimum total scaled RMSE
# This gives a zero penalty to the best model
complete[, penalty_factor := rmse_std_total - min(rmse_std_total)]

# Rank the models based on the penalty factor
complete <- complete[order(penalty_factor), ][, model_rank := 1:.N]


# Create a vector of model ids orderd by model rank, and excluding those beyond the 1SE threshold,
# and those that failed to converge. We'll use this vector below as candidate models for building the meta-model
complete[, eligible := within_1se & ar1_status]
complete[, model_step_tag := paste0('step', step, '_model', model_id)]
eligible_model_ids <- complete[eligible == TRUE, ][order(penalty_factor), model_step_tag]


# Some Dx plotting
complete %>% ggplot(aes(x = model_rank, y = penalty_factor, color = eligible)) + 
  geom_point(alpha = 0.3) + theme_minimal()

complete %>% dplyr::filter(eligible==TRUE) %>% ggplot(aes(x = model_rank, y = penalty_factor/sum(penalty_factor))) + 
  geom_point() + theme_minimal()


complete %>% filter(within_1se==TRUE) %>% 
  ggplot(aes(x = rmse_mean, y = rmse_trend_mean, color = penalty_factor)) +
  geom_point(alpha = 0.5) + 
  geom_point(data = complete[within_1se==FALSE & ar1_status==T, ], aes(x = rmse_mean, y = rmse_trend_mean), 
             color = 'darkred', color = 'firebrick', shape = 4) +
  theme_minimal() + scale_color_viridis() + 
  labs(title = 'Model fit by RMSE and trend RMSE', x = 'RMSE', y = 'Trend RMSE', color = 'Scaled RMSE total')


# Read in training data and ensure variables are formatted as they were during model building
allPoints <- fread(file.path(ROOT_DIR, 'FILEPATH_REMOVED_FOR_SECURITY', 'master_df.csv'))[age_group_id>2, ]
allPoints[, age_group_id := as.factor(age_group_id)]
allPoints[, ihme_loc_id:= as.factor(ihme_loc_id)]


# Define a function to pull predictions from the AR1 models and format them for comparison
pull_preds <- function(file, se = FALSE) {
  tmp <- read_feather(file)
  setDT(tmp)
  
  model <- gsub('_fold999_ar1_preds.feather', '', gsub(paste0(xval_dir, '/'), '', file))
  model_id <- as.integer(gsub('step[0-9]+_model', '', model))
  setnames(tmp, paste0(c('pred_', 'se_'), model_id) , c('pred', 'se'))
    
  if (!se) tmp[, se := NULL]
    
  tmp[, model := model]
  return(tmp)
}


# Read in the predictions from the best models 
# Note that we rbind then dcast to wide format as it's faster than merging all files
preds <- rbindlist(pblapply(complete[eligible == T, ar1_pred_file], pull_preds))
preds <- dcast.data.table(preds, demographic_id ~ model, value.var = c('pred'))

# Use that vector to order columns based on OOS fit such that the predictions from the
# best model are first and predictions from the worst model are last
preds <- preds[, .SD, .SDcols = c('demographic_id', eligible_model_ids)]

# Rename prediction columns from model_id to pred_model_id
pred_vars <- paste0('pred_', eligible_model_ids)
names(preds) <- c('demographic_id', pred_vars)

# Merge in the demographic information associated with each demographic id 
# (need to know age and year for subsetting below)
preds <- merge(copy(allPoints)[, .(location_id, ihme_loc_id, age_group_id, sex_id, year_id, demographic_id, iNTS_log)],
               preds, by = 'demographic_id', all = T)


# Split data into training and testing -- fold_id 999 indicates a run on the full data set #
fold  <- fread(file.path(ROOT_DIR, 'inputs', 'xvalidation_folds.csv'))
mdl_dt <- preds[year_id %in% fold$year_id & age_group_id != 2, .SD, .SDcols = c('year_id', 'iNTS_log', pred_vars)]

train_indices <- lapply(unique(fold$fold), function(f) which(mdl_dt$year_id %in% fold[fold == f & set == 'train', year_id]))
test_indices  <- lapply(unique(fold$fold), function(f) which(mdl_dt$year_id %in% fold[fold == f & set == 'test', year_id]))


# Pull the penalty factor for each model in the order of the best models (i.e. the column order of preds)
pf <- as.numeric(complete[sapply(eligible_model_ids, function(i) which(complete$model_step_tag == i)), penalty_factor])

# Pull the response and predictors for the best models (glmnet takes these as a vector and matrix, not a single dataframe)
y <- preds[year_id %in% fold$year_id & age_group_id != 2, iNTS_log]
x <- data.matrix(preds[year_id %in% fold$year_id & age_group_id != 2, .SD, .SDcols = pred_vars])

# Save the predictors, response, penalty factor, and fold indices for use in crossvalidataion
save(x, y, pf, train_indices, test_indices, file = file.path(xval_dir, 'xval_glmnet_data.RData'))


# passing as integers to simplify arguments and file names
# will divide by 100 to get actual alpha values
alphas <- seq(10L, 90L, 5L) 

# passing as integers to simplify arguments and file names
# will convert as lambda = 10^(x/10) in the model
lambdas <- seq(-40L, 0L, 1L) 


# Set parallel job launch parameters
project <- '-A proj_erf '
user <- Sys.getenv('USER')

slurm_output_dir <- paste0('-o FILEPATH_REMOVED_FOR_SECURITY', user, '/output/%x.o%j.out -e FILEPATH_REMOVED_FOR_SECURITY', user, '/errors/%x.e%j.err')
worker_script <- file.path(ROOT_DIR, 'code', '03b_xval_model_run_step2.R')
r_shell <- 'FILEPATH_REMOVED_FOR_SECURITY/execRscript.sh -s'
runtime <- '1:00:00'

THREADS <- 10
MEM <- '32G'

# Launch the jobs
for (alpha in alphas) { 
  for (lambda in lambdas) {
    jname <- paste('xval2', alpha, lambda, sep = '_')
    sys_sub <- paste0('sbatch ', project, slurm_output_dir, ' -J ', jname, ' -c ', THREADS, ' --mem=', MEM,  ' -t ', runtime, ' -p all.q')
    args <- paste(alpha, lambda, XVAL_VERSION)
    system(paste(sys_sub, r_shell, worker_script, args))
  }
}


# Compile all the fit stats from the crossvalidation runs
fit_stats <- rbindlist(lapply(alphas, function(a) 
  rbindlist(lapply(lambdas, function(l) 
    readRDS(file.path(xval_dir, paste0('step_2_fit_stats_alpha', a, '_lambda', l, '.rds')))))))

# Calculate mean fit stats across crossvalidation folds (excluding the full data run here)
fit_stat_means <- copy(fit_stats)[fold!=999, ][, .(mse = mean(mse), se = sd(mse)/7), by = .(alpha, lambda, std)]
fit_stat_means <- fit_stat_means[order(mse), ][, model_rank := 1:.N]
fit_stat_means[, within_1se := mse < min(mse) + se]

# Merge in the varibles included in the models run on the full dataset
fit_stat_means <- merge(fit_stat_means, fit_stats[fold==999, .(alpha, lambda, std, included)], by = c('alpha', 'lambda'), all.x = T)
fit_stat_means$n_included <- sapply(1:nrow(fit_stat_means), function(i) length(fit_stat_means$included[i][[1]]))


# Plot a heat map of the RMSE across all folds for each alpha and lambda
fit_stat_means %>% ggplot(aes(x = lambda, y = alpha, fill = sqrt(mse))) + geom_tile() + 
  scale_fill_viridis(option = 'B', name = 'RMSE') + theme_minimal() + scale_x_log10() +
  facet_wrap(~std) + xlab('Lambda') + ylab('Alpha') + ggtitle('RMSE by alpha and lambda') 


# Plot the relationship between the number of included predictors and the OOS RMSE
fit_stat_means %>% ggplot(aes(x = n_included, y = sqrt(mse), color = model_rank)) + geom_point(alpha = 0.25) +
  scale_color_viridis(option = 'B', name = 'Model Rank', trans = 'log') + theme_minimal() +
  xlab('Number of Included Predictors') + ylab('RMSE') + ggtitle('RMSE by number of included predictors')


# Run and save the best model
best_model <- which(fit_stat_means$mse==min(fit_stat_means$mse))
mdl_best <- glmnet(x, y, alpha = fit_stat_means$alpha[best_model], lambda = fit_stat_means$lambda[best_model], penalty.factor = pf, standardize = T)
saveRDS(mdl_best, file.path(xval_dir, 'glmnet_model.rds'))




