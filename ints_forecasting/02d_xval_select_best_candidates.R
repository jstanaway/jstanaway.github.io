# Project: iNTS Forecasting
#
# Purpose: 
# This script reads in the fit statistics from every model and fold from crossvalidation and assesses models based on errors 
# in OOS predictions in the point estimate (rmse) and the predicted trend (rmse_trend), and retains only those models with 
# RMSE values within 1SE of the best model with regard to both point estimates and trend. It then launches jobs to run those
# selected models on the full dataset. This script should only be run when all jobs launched in the previous step are complete.
#
# Plan for this step to take about a week: it launches far fewer jobs that the first crossvalidation launch, since there are 
# no folds and only only a subset of better models are run, but as these models are run on the full dataset each one takes
# quite a bit longer to complete. If the cluster is reasonably clear all jobs may start running within minutes, and while
# the average model may complete in 12 hours or so, a large number will take a few days, and a handful may take a week or 
# more. 
#
# Author: Jeff Stanaway (stanaway@uw.edu)
# Date created: April 2023
# Date last modified: 25 March 2024
#
# Libraries: data.table, tidyverse
# Dependencies: xval_model_run.R; xvalidation_folds.csv (created by 02a_crossvalidation_fold_maker.R)
# User inputs: XVAL_VERSION, ROOT_DIR, THREADS, MEM (you'll likely not need to change the last two)


rm(list = ls())


library(data.table)
library(tidyverse)


XVAL_VERSION <- 3

ROOT_DIR <- 'FILEPATH_REMOVED_FOR_SECURITY'
xval_dir <- file.path(ROOT_DIR, 'xvalidation', XVAL_VERSION)
graph_dir <- file.path(xval_dir, 'graphs')

dir.create(graph_dir, recursive = T, showWarnings = F)

user <- Sys.getenv("USER")


# Get list of all fit statistics files and read into a single data.table
fit_stat_files <- list.files(xval_dir, '^fit_stats_step.*csv$', full.names = F)
fit_stats <- rbindlist(lapply(fit_stat_files, function(fname) fread(file.path(xval_dir, fname))[, fname := fname]))

# Get the names of the fit statistics columns
stat_names <- setdiff(names(fit_stats), c('model', 'fname'))

# Extract model number and fold number from the file name
fit_stats[, mdl_fname := gsub('.csv', '.rds', gsub('fit_stats_', '', fname))]
fit_stats[, step := gsub('step', '', str_split_i(fit_stats$mdl_fname, "_", 1))]
fit_stats[, model_id := as.integer(gsub('model', '', str_split_i(fit_stats$mdl_fname, "_", 2)))]
fit_stats[, fold_num := as.integer(gsub('fold', '', gsub('.rds', '', str_split_i(fit_stats$mdl_fname, "_", 3))))]

# Calculate the number of complete fits for each model
fit_stats[, n_complete := .N, by = c('step', 'model_id')]
fit_stats[, max_complete := max(n_complete), by = c('step', 'model_id')]
n_folds <- max(fit_stats$n_complete)


# Collapse fold-specific stats to model means, for models in which all folds converged
complete <- fit_stats[n_complete == n_folds, ]
complete <- rbind(complete[, lapply(.SD, mean), by = c('step', 'model_id', 'model'), .SDcols = stat_names][, stat := 'mean'],
                  complete[, lapply(.SD, sd), by = c('step', 'model_id', 'model'), .SDcols = stat_names][, stat := 'sd'],
                  complete[, lapply(.SD, max), by = c('step', 'model_id', 'model'), .SDcols = stat_names][, stat := 'max'])

complete <- dcast.data.table(complete, step + model_id + model ~ stat, value.var = stat_names)                  


# Read in and merge the model details so we can see which variables were included in the best models
model_details <- rbindlist(lapply(unique(complete$step), function(step) {
  data.table(readRDS(file.path(xval_dir, paste0('step', step, '_model_details.rds'))))[, step := step]}))

complete <- merge(complete, model_details, by = c('step', 'model'), all.x = T)


# Drop duplicate models (note: there was a quirk in the algorithm that caused a few variable combinations to be duplicated in the cross-validation)
complete <- complete[, dup_index := 1:.N, by = c('step', 'model')][dup_index == 1, ]


# Determine which models are within 1 standard error of the best model
complete[, rmse_mean_1se := rmse_sd/sqrt(n_folds)]
complete[, rmse_trend_mean_1se := rmse_trend_sd/sqrt(n_folds)]

min_rmse_mean <- min(complete$rmse_mean)
min_rmse_se <- complete[which(complete$rmse_mean == min_rmse_mean), rmse_mean_1se]
rmse_mean_1se_threshold <- min_rmse_mean + min_rmse_se

min_rmse_trend_mean <- min(complete$rmse_trend_mean)
min_rmse_trend_se <- complete[which(complete$rmse_trend_mean == min_rmse_trend_mean), rmse_trend_mean_1se]
rmse_trend_mean_1se_threshold <- min_rmse_trend_mean + min_rmse_trend_se

complete[, within_1se := (rmse_mean < rmse_mean_1se_threshold) & (rmse_trend_mean < rmse_trend_mean_1se_threshold)]


# Write out the complete data.table
saveRDS(complete, file = file.path(xval_dir, 'step1_fit_stats_complete.rds'))



# Slurm Project and output
project <- "-A proj_erf "
slurm_output_dir <- paste0("-o FILEPATH_REMOVED_FOR_SECURITY", user, "/output/%x.o%j.out -e FILEPATH_REMOVED_FOR_SECURITY", user, "/errors/%x.e%j.err")
worker_script <- file.path(ROOT_DIR, 'code', '02c_xval_model_run.R')
r_shell <- "FILEPATH_REMOVED_FOR_SECURITY.sh -s"
runtime <- "380:00:00"

THREADS <- 4
MEM <- '64G'
EST_AR1_COEF <- F # TRUE if you want to estimate the AR1 coefficient, FALSE otherwise (keep F here)
FOLD_ID <- 999

### Determine which models need to be launched ##
# If the output file exists already then we don't need to run it again
complete[, full_fname := file.path(xval_dir, paste0('step', step, '_model', model_id, '_fold999.rds'))]
complete[, full_exists := file.exists(full_fname)]

# If the job is already running then we don't need to launch it again
is_running <- data.table(job_name = grep('^xval_.*_999$', system(paste0("squeue --me -o '%j'"), intern = T), value = T),
                         running = TRUE)
complete[, job_name := paste("xval", model_id, FOLD_ID, sep = '_')]
complete <- merge(complete, is_running, by = 'job_name', all.x = T)

complete[, to_launch := is.na(running) & !full_exists & within_1se]


# Loop through every combination of model and fold, and submit a job to run the cross-validation
for (i in which(complete$to_launch)) { 
  model_id <- complete$model_id[i]
  step <- complete$step[i]
  
  jname <- paste("xval", model_id, FOLD_ID, sep = '_')
  sys_sub <- paste0("sbatch ", project, slurm_output_dir, " -J ", jname, " -c ", THREADS, " --mem=", MEM,  " -t ", runtime, " -p long.q")
  args <- paste(model_id, FOLD_ID, step, XVAL_VERSION, EST_AR1_COEF)
  system(paste(sys_sub, r_shell, worker_script, args))
}


