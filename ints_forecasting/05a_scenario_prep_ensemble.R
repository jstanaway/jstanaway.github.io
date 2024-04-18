# Project: iNTS Forecasting
#
# Purpose: 
# This script builds a workspace that can be shared externally and used to make scenario-based 
# iNTS forecasts. It is a stripped down version of the full draw-level forecasting script, and
# can be run outside of IHME's computing environment. To make this tractable, this script compiles
# all necessary functions and data into a single RData file, and the forecasts are run on point
# estimates rather than draws (i.e. no uncertainty).
#
# Author: Jeff Stanaway (stanaway@uw.edu)
# Date created: April 2023
# Date last modified: 2 April 2024


# Clear the workspace
rm(list = ls())

# Load libraries
library(boot)
library(caTools)
library(data.table)
library(glmnet)
library(mgcv)
library(pbapply)
library(scam)
library(tidyverse)
library(zoo)


# Load shared functions
SHARED_FUN_DIR <- 'FILEPATH_REMOVED_FOR_SECURITY'
source(file.path(SHARED_FUN_DIR, 'get_age_metadata.R'))
source(file.path(SHARED_FUN_DIR, 'get_demographics.R'))
source(file.path(SHARED_FUN_DIR, 'get_draws.R'))
source(file.path(SHARED_FUN_DIR, 'get_ids.R'))
source(file.path(SHARED_FUN_DIR, 'get_location_metadata.R'))
source(file.path(SHARED_FUN_DIR, 'get_outputs.R'))
source(file.path(SHARED_FUN_DIR, 'get_population.R'))


# Set release constants
RELEASE <- 6
REF_YEAR  <- 2019
XVAL_VERSION <- 3
fhs_years <- REF_YEAR:2100


# Establish GBD draw directory based on REF_YEAR (REF_YEAR==GBD cycle year)
INTS_DIR <- 'FILEPATH_REMOVED_FOR_SECURITY'
xval_dir <- file.path(INTS_DIR, 'xvalidation', XVAL_VERSION)
raw_dir <- file.path(INTS_DIR, 'raw_extractions')
gbd_draw_dir <- file.path('FILEPATH_REMOVED_FOR_SECURITY', REF_YEAR)

# Load custom functions
source(file.path(INTS_DIR, 'code', 'helper_functions', 'scenario_helper_functions.R'))








############################################################################
#                                INPUT PREP                                #
############################################################################

# Get location metadata
loc_meta <- get_location_metadata(location_set_id = 35, release_id = RELEASE)

# Get age metadata
age_meta <- get_age_metadata(release_id = RELEASE)
age_meta[, age_mid := (age_group_years_start + age_group_years_end)/2]

# Get measures
measures <- get_ids(table = 'measure')[, .(measure_id, measure)]

# Read in inputs
all_points <- fread(file.path(INTS_DIR, 'raw_extractions', 'merged_covariates', 'transformed_cov', 'master_df.csv'))
all_points <- all_points[year_id %in% 1990:2100 & location_id %in% loc_meta[level==3, location_id]]


pop_means  <- rbindlist(pblapply(loc_meta[level<=3, ihme_loc_id], function(loc) {
  get_raw_extraction(file.path(raw_dir, 'population', paste0(loc, '.csv')), c('age_group_weight_value', 'population'), fhs_years)}))


pop <- get_population(location_id = loc_meta[most_detailed==1, location_id], 
                      age_group_id = c(age_meta$age_group_id, 22, 27), 
                      sex_id = 1:3, year_id = REF_YEAR, 
                      release_id = RELEASE)[, run_id := NULL]


# Read in draws of disability weights
dw <- fread(file.path(INTS_DIR, 'inputs', 'weights_draws.csv'))[, lapply(.SD, mean), .SDcols = 'weight']


# Prep HIV prevalence
hiv_prev   <- rbindlist(pblapply(loc_meta[level==3, ihme_loc_id], function(loc) {
  get_raw_extraction(file.path(raw_dir, 'hiv_prev', paste0(loc, '.csv')), 'hiv_prev', fhs_years)}))


# Get reference estimates for baseline shifts 
ref_ests <- rbindlist(pblapply(loc_meta[level==3, location_id], draws_to_means, ref_year = REF_YEAR))

cc_ests  <- get_outputs('cause', cause_id = 959, location_id = loc_meta[level==3, location_id], age_group_id = age_meta$age_group_id, sex_id = 1:2, year_id = REF_YEAR, 
                    measure_id = c(1,6), metric_id = 3, release_id = RELEASE)[, .(year_id, location_id, age_group_id, sex_id, measure, val)]

cc_ests[age_group_id==2, val := 0]
cc_ests[, hiv := ifelse(measure=='incidence', 2, 0)]

# Load reference life table
life_table <- fread(file.path(INTS_DIR, 'inputs', 'standard_lt.csv'))

# Get duration
duration <- ((2/0.21)-1)/365.25  #duration parameters of negative binomial derived from ABC analysis (see GBD code/docs for details)


# Read in models
meta_mdl <- readRDS(file.path(xval_dir, 'glmnet_model.rds'))
meta_mdl_coefs <- coef(meta_mdl)
weak_learners  <- meta_mdl_coefs@Dimnames[[1]][-1][meta_mdl_coefs@i]
weak_learners <- as.integer(gsub('pred_', '', weak_learners))


# The model objects are huge because they contain all of the original training data,
# fitted values, residuals, and other components that are unnecessary for prediction 
# on new data.  The next block of code drops unnecessary components of the model
# object and shrinks others that are needed to a minimum size

# These are model object components to delete altogether
to_drop <- c('y', 'linear.predictors', 'weights', 'fitted.values', 'residuals', 'effects', 
            'prior.weights', 'offset', 'na.action', 'std.rsd', 'R', 'Ve', 'Ve.t', 
            'edf', 'edf1', 'var.summary', 'formula', 'control', 'beta', 'beta.t', 'p.ident',
            'call', 'sp', 'bfgs.info')

# Define the function that shrinks the model objects
shrink_model <- function(model_id) {
  # Read in the model object
  mdl <- readRDS(file.path(xval_dir, paste0('step1_model', model_id, '_fold999_ar1.rds')))
  
  # Drop unnecessary components
  for (x in to_drop) {
    mdl[x] <- NULL
  }
  
  # Drop large and unnecessary subcomponents of the smooths
  for (i in 1:(length(mdl$smooth) - 1)) {
     mdl$smooth[[i]]$Xdf1 <- NULL
     mdl$smooth[[i]]$Xdf2 <- NULL
   }
  
  # Reduce the model component that contains input data to only the minimum necessary
  tmp <- data.table(mdl$model)
  tmp[, touse := 1:.N, by = c('(AR.start)', 'ihme_loc_id', 'age_group_id', 'sex_id')]
  tmp <- tmp[touse == 1, ][, touse := NULL]
  mdl$model <- tmp
  
  # Save and return the shrunk object
  saveRDS(mdl, file.path(xval_dir, paste0('step1_model', model_id, '_fold999_ar1_shrink.rds')))
  return(mdl)
}

# Apply the shrink_model function to all weak learner models
weak_learner_mdls <- pblapply(weak_learners, shrink_model)
saveRDS(weak_learner_mdls, file.path(xval_dir, 'weak_learner_list.rds'))


# Read in the PCA models
pca_model_wash <- readRDS(file.path(INTS_DIR, 'inputs', 'pca_model_wash_sdi.rds'))
pca_model_undernutrition <- readRDS(file.path(INTS_DIR, 'inputs', 'pca_model_undernutrition.rds'))

# Read in the HIV RR model
rr_model <- readRDS(file.path(INTS_DIR, 'inputs', 'hivRrModel.rds'))

# Read in the CFR model
load(file.path(INTS_DIR, 'inputs', 'cfrModel.RData'))
cfr_model <- cfrModel
rm(cfrModel)





### CREATE CFR PREDICTION LOOKUP TABLE ###
sdi_resolution <- 0.001

cfr <- merge(data.frame(sdi = seq(0.2, 0.9, sdi_resolution)), data.frame(age_meta[, .(age_group_id, age_mid)]), all = T)
cfr <- merge(cfr, data.frame(estPrHiv = 0:1), all = T)
setDT(cfr)
setnames(cfr, 'age_mid', 'ageMid')

cfr[, logit_cfr_pred := as.numeric(predict.gam(cfr_model, cfr, se.fit = F))]
cfr[, cfr_pred := inv.logit(logit_cfr_pred)]
cfr[, sdi_int := as.integer(round(sdi / sdi_resolution))]

setnames(cfr, 'estPrHiv', 'est_pr_hiv')
cfr <- cfr[, .(sdi_int, age_group_id, est_pr_hiv, cfr_pred)]


# Clean up full prediction data frame
in_data <- copy(all_points)[year_id >= REF_YEAR, ]
in_data <- in_data[, .SD, .SDcols = c('location_id', 'ihme_loc_id', 'year_id', 'age_group_id', 'sex_id', 
                                    'sdi', 'water', 'sanitation', 'hygiene', 'hiv', 'malaria', 
                                    'underweight', 'wasting','stunting', 'ln_SEV_diarrhea', 'demographic_id')]
setnames(in_data, 'hiv', 'hiv_mort')

in_data <- merge(in_data, hiv_prev, by = c('location_id', 'year_id', 'age_group_id', 'sex_id'), all.x = T)



# Keep only needed rows and columns of loc_meta
loc_meta <- loc_meta[level<=3, ][, .(location_id, parent_id, path_to_top_parent, level, location_name, super_region_id, super_region_name, region_id, region_name, ihme_loc_id)]


# Save the workspace
save(age_meta, cc_ests, cfr, dw, in_data, life_table, loc_meta, meta_mdl, meta_mdl_coefs, 
     pca_model_undernutrition, pca_model_wash, pop_means, ref_ests, rr_model, weak_learner_mdls, 
     duration, fhs_years, REF_YEAR, RELEASE, sdi_resolution, weak_learners, aggregate_locs, 
     baseline_shift, ints_scenario, file = file.path(INTS_DIR, 'code/scenario', 'scenario_workspace_ensemble.RData'), 
     compress = 'gzip')

