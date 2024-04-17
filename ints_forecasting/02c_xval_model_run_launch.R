# Project: iNTS Forecasting
#
# Purpose: 
# Given the size of the dataset and complexity of the models, we need to parallelize crossvalidation 
# for this analysis to be tractable in a reasonable amount of time.  This script launches the 
# parallelized crossvalidation jobs.  The script that's doing the actual work is xval_model_run.R, and 
# that is built to be robust to changes between cycles and should require no changes or maintenance 
# unless we're changing our methods. This launch script can also optionally run 02b_xval_model_maker.R 
# which will create an object (and save corresponding .rds files) containing details of every model 
# variant to be tested in crossvalidation. This will need to be run once for a crossvalidation version 
# but can be disabled for later runs to save a bit of computation. Finally, after the launch code,
# you'll find code to monitor the status of launched jobs and check for run errors.
#
# Plan for the crossvalidation run to take a while: the average model will run on a single fold in a 
# few hours, but some may take a few days, and with >8000 model-fold combinations you can plan for this 
# to take a few days when the cluster is reasonably clear or substantially longer when it's busy -- 
# a week is typical. You'll need to wait for all crossvalidation jobs to complete before moving on to 
# the next step.
#
# Author: Jeff Stanaway (stanaway@uw.edu)
# Date created: April 2023
# Date last modified: 25 March 2024
#
# Libraries: data.table
# Dependencies: xval_model_maker.R, xval_model_run.R; xvalidation_folds.csv (created by 02a_crossvalidation_fold_maker.R)
# User inputs: XVAL_VERSION, REF_YEAR, MAKE_MODEL_FILE, ROOT_DIR, THREADS, MEM, STEP, EST_AR1_COEF


rm(list = ls())

library(data.table)


### SET RUN PARAMETERS ###

XVAL_VERSION <- 3    # Set the version here
REF_YEAR <- 2019     # Set the reference year here -- this is the year corresponding to the GBD cycle
MAKE_MODEL_FILE <- F # TRUE if you want to make the model file, FALSE otherwise.  You need to run this once to make the file, then set to FALSE

ROOT_DIR <- 'FILEPATH_REMOVED_FOR_SECURITY'
xval_dir <- file.path(ROOT_DIR, 'xvalidation', XVAL_VERSION)

dir.create(xval_dir, recursive = T)

user <- Sys.getenv("USER")


# Either make the file containing all possible model formulae, or read it in
if (MAKE_MODEL_FILE == T) {
  source(file.path(ROOT_DIR, 'code', 'xval_model_maker.R'))
  step1_models <- make_models(xval_dir)
} else {
  step1_models <- readRDS(file.path(xval_dir, 'step1_models.rds'))
}
 
# Read in the cross-validation folds
folds <- fread(file.path(ROOT_DIR, 'inputs', 'xvalidation_folds.csv'))



### SUBMIT THE CROSSVALIDATION JOBS TO THE CLUSTER ###
# Slurm Project and output
project <- "-A PROJECT_REMOVED_FOR_SECURITY "
slurm_output_dir <- paste0("-o FILEPATH_REMOVED_FOR_SECURITY", user, "/output/%x.o%j.out -e FILEPATH_REMOVED_FOR_SECURITY", user, "/errors/%x.e%j.err")
worker_script <- file.path(ROOT_DIR, 'code', '02c_xval_model_run.R')
r_shell <- "FILEPATH_REMOVED_FOR_SECURITY -s"
runtime <- "380:00:00" 

launch_date <- Sys.Date()

THREADS <- 4
MEM <- '32G'
STEP <- '1'
EST_AR1_COEF <- F # TRUE if you want to estimate the AR1 coefficient, FALSE otherwise


# Loop through every combination of model and fold, and submit a job to run the cross-validation
for (model_id in 1:nrow(step1_models)) { 
  for (fold_id in unique(folds$fold)) {
    jname <- paste("xval", model_id, fold_id, sep = '_')
    sys_sub <- paste0("sbatch ", project, slurm_output_dir, " -J ", jname, " -c ", THREADS, " --mem=", MEM,  " -t ", runtime, " -p long.q")
    args <- paste(model_id, fold_id, STEP, XVAL_VERSION, EST_AR1_COEF)
    system(paste(sys_sub, r_shell, worker_script, args))
  }
}



### CHECK THE STATUS OF THE JOBS ###

setDT(step1_models)
step1_models[, mdl_num := paste0('model', 1:.N)]
step1_models[, dup_count := .N, by = 'model']

# Get the names of the jobs that are currently running
is_running <- data.table(job_name = grep('xval_', system(paste0("squeue --me -o '%j'"), intern = T), value = T),
                         running = TRUE)


# Make a data table will all jobs
status <- data.table(expand.grid(model_id = model_ids, fold_id = unique(folds$fold))) #, complete = F, running = F)
status[, job_name := paste('xval', model_id, fold_id, sep = '_')]

# Merge the running jobs with the jobs that have been submitted
status <- merge(status, is_running, by = 'job_name', all.x = T)
status[is.na(running), running := F]

# Determine which jobs are complete based on existence of output file
status[, file := file.path(xval_dir, paste0('step', STEP, '_model', model_id, '_fold', fold_id, '.rds'))]
status[, complete := file.exists(file)]

# Print the status of the jobs
table(status$complete, status$running, useNA = 'ifany', deparse.level = 2)


