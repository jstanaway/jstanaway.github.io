# Project: iNTS Forecasting
#
# Purpose: 
# This script determines which models launched in the previous step completed and launches a job to run the 
# corresponding model with AR1 errors to account for residual autocorrelation.  This extra step is required 
# because, while SCAM can account for autocorrelation, it cannot estimate the autocorrelation coefficient 
# and requires the user to pass that in as an argument.The models produced here will be the final models 
# used as weak learners to be stacked in the meta-model or strong learner built in the next steps.  
# Note, this script determines which models are ready to be run (i.e. those in which the full model from the
# previous step completed), and have not yet been launched (i.e. it determines which ar1 models are either 
# currently running or have already completed).  Therefore, you don't have to wait until all full models are 
# complete to run it, and it can be re-run as needed to launch all relevant models -- useful if you're trying 
# to push things through as quickly as possible and don't want to wait for everything from the previous step 
# to finish before moving on. 
#
# Plan for this step to take about a week: if launches far fewer jobs that the first crossvalidation launch, since there are 
# no folds and only only a subset of better models are run, but as these models are run on the full dataset each one takes
# quite a bit longer to complete. If the cluster is reasonably clear all jobs may start running within minutes, and while
# the average model may complete in 12 hours or so, a large number will take a few days, and a handful may take a week or 
# more. 
#
# Author: Jeff Stanaway (stanaway@uw.edu)
# Date created: April 2023
# Date last modified: 25 March 2024
#
# Libraries: data.table, 
# User inputs: XVAL_VERSION, ROOT_DIR, THREADS, MEM, STEP, EST_AR1_COEF (the last four are unlikely to need changing)

rm(list = ls())

library(data.table)

XVAL_VERSION <- 3

ROOT_DIR <- 'FILEPATH_REMOVED_FOR_SECURITY'
xval_dir <- file.path(ROOT_DIR, 'xvalidation', XVAL_VERSION)

user <- Sys.getenv('USER')

# Get the list of models that have completed running on the full dataset (ie fold 999)
complete_models <- data.table(mdl_fname = list.files(xval_dir, 'step1_.*_fold999.rds'))

# Extract model information from the filename
complete_models[, step := gsub('step', '', str_split_i(complete_models$mdl_fname, '_', 1))]
complete_models[, model_id := as.integer(gsub('model', '', str_split_i(complete_models$mdl_fname, '_', 2)))]
complete_models[, fold_num := as.integer(gsub('fold', '', gsub('.rds', '', str_split_i(complete_models$mdl_fname, '_', 3))))]

# Determine if the AR1 model has already been run
complete_models[, ar1_fname := file.path(xval_dir, paste0('step', step, '_model', model_id, '_fold', fold_num, '_ar1.rds'))]
complete_models[, ar1_exists := file.exists(ar1_fname)]
complete_models[, job_name := paste('ar1_xval', model_id, '999', sep = '_')]

is_running <- data.table(job_name = grep('ar1_', system(paste0("squeue --me -o '%j'"), intern = T), value = T),
                         running = TRUE)

complete_models <- merge(complete_models, is_running, by = 'job_name', all.x = T)

# Select only those models that have completed the full model but for which the AR1 model has
# not completed and is not running to be launched
complete_models[, to_launch := is.na(running) & !ar1_exists]

# Slurm Project and output
project <- '-A PROJECT_REMOVED_FOR_SECURITY'
slurm_output_dir <- paste0('-o FILEPATH_REMOVED_FOR_SECURITY', user, '/output/%x.o%j.out -e FILEPATH_REMOVED_FOR_SECURITY', user, '/errors/%x.e%j.err')
worker_script <- file.path(ROOT_DIR, 'code', '02c_xval_model_run.R')
r_shell <- 'FILEPATH_REMOVED_FOR_SECURITY.sh -s'
runtime <- '380:00:00'

THREADS <- 4
MEM <- '64G'
EST_AR1_COEF <- TRUE # TRUE if you want to estimate the AR1 coefficient, FALSE otherwise (keep T here)
FOLD_ID <- 999
  

# Loop through every row with a model to launch and submit the job to run the full AR1 model
for (i in which(complete_models$to_launch)) { 
  model_id <- complete_models$model_id[i]
  step <- complete_models$step[i]
  
  jname <- paste('ar1_xval', model_id, FOLD_ID, sep = '_')
  sys_sub <- paste0('sbatch ', project, slurm_output_dir, ' -J ', jname, ' -c ', THREADS, ' --mem=', MEM,  ' -t ', runtime, ' -p long.q')
  args <- paste(model_id, FOLD_ID, step, XVAL_VERSION, EST_AR1_COEF)
  system(paste(sys_sub, r_shell, worker_script, args))
}

