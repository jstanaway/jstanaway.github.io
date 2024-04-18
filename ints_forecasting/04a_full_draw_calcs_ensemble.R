# Project: iNTS Forecasting
#
# Purpose: 
# This script makes draw-level predictions of iNTS incidence, prevalence, mortality, YLLs, YLDs, and DALYs 
# for all years (up to 2100), age groups, and sexes for a single location.  It is launched by 
# 04a_full_draw_calcs_ensemble_launch.R.  This is a pretty long piece of code as it's doing a lot of work: 
# processing these draw-level predictions requires reading in a large number of inputs with disperate formatting, 
# processing all forecasts at the draw-level to correctly propogate uncertainty, predicting based on an ensemble
# of models, ensuring continuity of GBD and forecasted estimates (including continuity of uncertainty), and 
# replicating much of the work normally handled by central comp (i.e. gap metric estimation and age/sex aggregation).  
# Broadly, the structure is to read in and process all inputs; make model-based predictions of incidence, HIV PAFs, 
# and case-fatality; calculate prevalence, mortality, YLLs, YLDs, and DALYs based on those model-based predictions; 
# and calculate age and sex aggregates (note that because this is parallelized by location, location aggregation
# is not done here, but in the next step of the code pipeline).
#
# Code outline:
#  1) Set up the environment: read in arguments, load packages and functions, establish directories, 
#     establish estimation years, and set seeds
#  2) Read in and process inputs:
#     a) Get location, age, and measure metadata
#     b) Read in the model objects: we'll read in PCA models, HIV RR model, case fatality model, incidence meta-model, 
#        and information about weak learners (we'll actually read those in only as needed in the incidence prediction 
#        section to save memory, since these model objects are very large)
#     c) Bring in draw-level covariate forecasts from FHS and calculate derived variables at the draw level (e.g. we'll 
#        estimate values of PCA derived variables, transformed variables, etc). We'll use these in conjuction with the 
#        models to make draw-level forecasts in a few steps
#     d) Bring in draw-level GBD estimates produced by central comp. We'll use these as the starting point from which 
#        to run out the forecasts to ensure temporal continuity between GBD and FHS estimates
#     e) Bring in draw-level GBD estimates pre-central comp.  These estimates are HIV-specific (HIV-specificity not 
#        included in central comp runs) and are used to ensure temporal continuity between GBD and FHS estimates related to
#        HIV attribution and HIV-attributable fatal burden that is not included in 'official' GBD estimates (iNTS deaths 
#        that are attributable to HIV are considered HIV deaths within the GBD framework)
#  3) Incidence estimation: we first make draw-level predictions for each weak learner, then for the final meta-model; 
#     then apply the baseline shift
#  4) HIV-attribution estimation: we forecast the proportion of iNTS cases attributable to HIV using a RR-PAF approach, 
#     apply a baseline shift to ensure temporal continuity with GBD estimates, and then apply the proportions to split 
#     incidence forecasts into cases attributable to HIV and cases not attributable to HIV
#  5) Case-fatality estimation: we forecast iNTS case fatality by HIV status and apply a baseline shift to those estimates
#  6) Mortality estimation: we forecast iNTS mortality as the product of incidence and case-fatality by HIV status. 
#     Because we've already applied the baseline shift to incidence, HIV-attribution, and case-fatality we don't need to 
#     apply another shift to ensure continuity with raw (i.e. pre-CoDCorrect) GBD mortality estimates; however, we do need 
#     a baseline shift to ensure continuity with official CodCorrected estimates, so we apply that shift here.  
#     This shift not only ensures temporal continuity, but also simulates the adjustment effect of CoDCorrect on the forecasts
#  7) YLL estimation: we forecast YLLs as the product of mortality and the reference life expecancy at the age of death.  
#     We estimate two versions here: one that is the raw YLL forecast and one that is 'CodCorrected'
#  8) YLD estimation: we forecast prevalence as the product of incidence and duration, and YLDs as the product of prevalence 
#     and disability weights.  Since all inputs have been shifted, no baseline shift is needed here.
#  9) DALY estimates: we forecast DALYs as the sum of YLLs and YLDs and estimate both CoDCorrected and raw versions.  
#     Since all inputs have been shifted, no baseline shift is needed here.
# 10) Age and sex aggregation: we aggregate age- and sex- specific estimates to get estimates for both sexes combined, 
#     all-ages, and age-standardized
# 11) Export results including a draw file and a file with summary statistics (point estimates with 95% UIs)
#
# Author: Jeff Stanaway (stanaway@uw.edu)
# Date created: February 2023
# Date last modified: 29 March 2024
#
# Libraries: data.table, tidyverse, boot, caTools, scam, mgcv, glmnet, zoo, pbapply
# Shared functions: get_age_metadata.R, get_demographics.R, get_draws.R, get_ids.R, get_location_metadata.R,
#                   get_outputs.R, get_population.R
# Custom functions: add_pr_comps.R, agg_subnat.R, baseline_shift_ensemble.R, compile_detailed_covariate_draws.R,
#                   prep_draws.R, rank_draws.R, read_covariate_draws.R (all in the code/helper_function folder) 
# User inputs: MAX_FHS_DRAW, SHARED_FUN_DIR, INTS_DIR 

rm(list = ls())


# Determine if this is an interactive run (i.e. for development and debugging) or not
INTERACTIVE = commandArgs()[2]=='--interactive'

print(commandArgs())

if (!INTERACTIVE) {  
  arg <- commandArgs(trailingOnly = T)
  LOC_ISO <- arg[1]
  RELEASE <- arg[2]
  OUT_DIR <- arg[3]
  XVAL_VERSION <- arg[4]
  DX <- F

} else { # if this is an interactive run, set the parameters here
  LOC_ISO <- 'LSO'
  RELEASE <- 6
  OUT_DIR <- 'FILEPATH_REMOVED_FOR_SECURITY'
  XVAL_VERSION <- 3
  DX <- T
}


MAX_FHS_DRAW <- 499L
SHARED_FUN_DIR <- 'FILEPATH_REMOVED_FOR_SECURITY'
INTS_DIR <- 'FILEPATH_REMOVED_FOR_SECURITY'

# Load libraries
library(data.table)
library(tidyverse)
library(boot)
library(caTools)
library(scam)
library(mgcv)
library(glmnet)
library(zoo)
library(pbapply)

# Load shared functions
source(file.path(SHARED_FUN_DIR, 'get_age_metadata.R'))
source(file.path(SHARED_FUN_DIR, 'get_demographics.R'))
source(file.path(SHARED_FUN_DIR, 'get_draws.R'))
source(file.path(SHARED_FUN_DIR, 'get_ids.R'))
source(file.path(SHARED_FUN_DIR, 'get_location_metadata.R'))
source(file.path(SHARED_FUN_DIR, 'get_outputs.R'))
source(file.path(SHARED_FUN_DIR, 'get_population.R'))

# Load custom functions
helper_fun_dir <- file.path(INTS_DIR, 'code', 'helper_functions')
for (fun in list.files(helper_fun_dir, '.R', full.names = T)) source(fun)

# Get year list
gbd_years <- get_demographics('epi', release_id = RELEASE)$year_id
ref_year  <- max(gbd_years)
fhs_years <- ref_year:2100
all_years <- unique(c(gbd_years, fhs_years))

# Establish GBD draw directory based on ref_year (ref_year==GBD cycle year)
gbd_draw_dir <- paste0('FILEPATH_REMOVED_FOR_SECURITY', ref_year)
extraction_dir <- file.path(INTS_DIR, 'raw_extractions')
xval_dir <- file.path(INTS_DIR, 'xvalidation', XVAL_VERSION)


# We want random seeds for each point in the code where we need to create independent draws, 
# but want them to be consistent across all jobs.  The code below does that for us
set.seed(2065336) # Starting seed taken from serial number of a random bill in my wallet
seed_list <- sample(1:.Machine$integer.max, 99) # randomly select seeds from all possible values that will be used sequentially below



############################################################################
#                             GET META DATA                                #
############################################################################

# Get location metadata
loc_meta <- get_location_metadata(location_set_id = 35, release_id = RELEASE)  #[grep(LOC_ISO, ihme_loc_id), ]
loc_id <- loc_meta[ihme_loc_id == LOC_ISO, location_id]

# Determine if location is most detailed
detailed_loc_ids <- loc_meta[most_detailed == 1 & grepl(paste0(',', loc_id, ','), path_to_top_parent), location_id]
all_loc_ids <- c(loc_id, detailed_loc_ids)
if (length(detailed_loc_ids) == 0) detailed_loc_ids <- all_loc_ids

# Get age metadata
age_meta <- get_age_metadata(release_id = RELEASE)
age_meta[, age_mid := (age_group_years_start + age_group_years_end)/2]

# Get measures
measures <- get_ids(table = 'measure')[, .(measure_id, measure)]



############################################################################
#                           GET MODEL OBJECTS                              #
############################################################################
# Read in meta-model and get weak learner information 
# (we'll read in the weak learners later, as we use them, to save memory)
meta_mdl <- readRDS(file.path(xval_dir, 'glmnet_model.rds'))
meta_mdl_coefs <- coef(meta_mdl)
weak_learners  <- meta_mdl_coefs@Dimnames[[1]][-1][meta_mdl_coefs@i]
weak_learners <- as.integer(gsub('pred_', '', weak_learners))

# Read in the PCA models
pca_model_wash <- readRDS(file.path(INTS_DIR, 'inputs', 'pca_model_wash_sdi.rds'))
pca_model_undernutrition <- readRDS(file.path(INTS_DIR, 'inputs', 'pca_model_undernutrition.rds'))

# Read in the HIV relative risk and case-fatality models
hiv_rr_model <- readRDS(file.path(INTS_DIR, 'inputs', 'hivRrModel.rds'))
load(file.path(INTS_DIR, 'inputs', 'cfrModel.RData'))
cfr_model <- cfrModel
rm(cfrModel)


############################################################################
#                 BRING IN DRAW-LEVEL COVARIATE FORECASTS                  #
############################################################################
# Create the covariate_draws data frame that contains draws for all variables used in incidence models
# Note that, while incrementally building out a data frame in a loop is typically not very efficient,
# because we need to need to determine the common variables between each covariate data frame, it
# works out to be a cleaner solution in this case

covar_list <- c('sdi', 'water', 'sanitation', 'hygiene', 'hiv', 'malaria', 'underweight', 'stunting', 'wasting')

for (i in 1:length(covar_list)) {
  if (i == 1) {
    covariate_draws <- read_covariate_draws(covar_list[i])
  } else {
    draw_tmp  <- read_covariate_draws(covar_list[i])
    by_vars   <- names(covariate_draws)[names(covariate_draws) %in% names(draw_tmp)]
    covariate_draws <- merge(covariate_draws, draw_tmp, by = by_vars, all = T)
  }
}

# We used the same seed to resample all covariates to maintain correlation structure
# between them (this occured in the 'read_covariate_draws' function in the loop above).  
# Now shift to next seed
seed_list <- seed_list[-1]

# Underweight, stunting, and wasting are only relevant for a couple of age groups, 
# so it's missing otherwise.  We set all those missing values to zero
undernutrition_vars <- c('underweight', 'wasting', 'stunting')
covariate_draws[, c(undernutrition_vars) := lapply(.SD, function(x) {ifelse(is.na(x), 0, x)}), .SDcols = undernutrition_vars]

# Calculate age-standardized under-nutrition variables
covariate_draws <- merge(covariate_draws, age_meta[, .(age_group_id, age_group_weight_value)], by = 'age_group_id', all.x = T)
covariate_draws[, paste0(undernutrition_vars, '_age_std') := lapply(.SD, function(x) sum(x * age_group_weight_value)), 
           by = c('location_id', 'year_id', 'sex_id', 'draw'), .SDcols = undernutrition_vars]

# Create 5-year malaria avg with uncertainty
# Note: coercion to numeric necessary for locations with no malaria, as those files 
# read in as integers, causing errors in runmean
covariate_draws <- covariate_draws[order(location_id, age_group_id, sex_id, draw, year_id)]
covariate_draws[, malaria_avg5 := runmean(as.numeric(malaria), k=5, endrule = 'mean', align = 'left', alg = 'exact'), 
          by = c('location_id', "age_group_id", "sex_id", 'draw')]

# Convert hiv and malaria to fourth root
root_vars <- c('malaria_avg5', 'hiv')
covariate_draws[, paste0(root_vars, '_quart_root') := lapply(.SD, nthroot, 4), .SDcols = root_vars]

# Create PCA variables
covariate_draws <- merge(covariate_draws, add_pr_comps(pca_model_wash, 'wash_sdi'), 
                   by = c('location_id', 'year_id', 'age_group_id', 'sex_id', 'draw'), all.x = T)

covariate_draws <- merge(covariate_draws, add_pr_comps(pca_model_undernutrition, 'undernutrition'), 
                   by = c('location_id', 'year_id', 'age_group_id', 'sex_id', 'draw'), all.x = T)


# Read in forecasted HIV prevalence draws
hiv_prev_draws <- fread(file.path(extraction_dir, 'hiv_prev', paste0(LOC_ISO, '.csv')))
hiv_prev_draws <- hiv_prev_draws[, .(location_id, age_group_id, sex_id, year_id, draw, hiv_prev)][year_id %in% all_years, ]

# Read in forecasted population draws
pop_draws <- fread(file.path(extraction_dir, 'population', paste0(LOC_ISO, '.csv')))
pop_draws <- pop_draws[, .(location_id, age_group_id, sex_id, year_id, draw, population, age_group_weight_value)][year_id %in% all_years & sex_id<3,  ]





############################################################################
#                   BRING IN GBD DRAWS FROM CENTRAL COMP                   #
############################################################################
# Pull draws of incidence and mortality from central comp (i.e. official GBD outputs after COMO and CoDCorrect)
mr_draws_cc <- get_draws('cause_id', gbd_id = 959, location_id = loc_id, age_group_id = c(age_meta$age_group_id, 22, 27), 
                         sex_id = 1:3, year_id = ref_year, measure_id = c(1, 4), metric_id = 1, source = 'codcorrect', release_id = RELEASE)

inc_draws_cc <- get_draws('cause_id', gbd_id = 959, location_id = loc_id, age_group_id = c(age_meta$age_group_id, 22, 27), 
                       sex_id = 1:3, year_id = ref_year, measure_id = c(3,5,6), metric_id = 3, source = 'como', release_id = RELEASE)


# Forecasts have fewer draws than GBD; Sample needed number of draws from GBD & reshape to long by draws
max_gbd_draw <- max(as.integer(gsub('draw_', '', grep('draw_', names(mr_draws_cc), value = T))))

# Set seed for reproducibility & shift to next seed
set.seed(seed_list[1])
seed_list <- seed_list[-1]
draw_sample <- data.table(draw = 0:MAX_FHS_DRAW, gbd_draw = sample(0:max_gbd_draw, (MAX_FHS_DRAW+1), replace = F))


# Ensure that the draws are in the correct format 
# (melt to long, clean up variables, sample correct number of draws)
inc_draws_cc <- prep_draws(inc_draws_cc, draw_sample)
mr_draws_cc  <- prep_draws(mr_draws_cc, draw_sample)


# Mortality draws in are absolute count space; merge in pop to convert to rates
mr_draws_cc <- merge(mr_draws_cc,
                     get_population(location_id = all_loc_ids, age_group_id = c(age_meta$age_group_id, 22, 27), 
                                    sex_id = 1:3, year_id = ref_year, release_id = RELEASE)[, run_id := NULL], 
                     by = c('location_id', 'age_group_id', 'sex_id', 'year_id'), all.x = T) 

mr_draws_cc[, `:=` (mr = death / population, yll_rate = yll / population, hiv = 0)]



# Read in draws of disability weights
dw <- fread(file.path(INTS_DIR, 'inputs', 'weights_draws.csv'))
setnames(dw, c('draw', 'weight'), c('gbd_draw', 'dw'))
dw <- merge(dw, draw_sample, by = 'gbd_draw', all.y = T)[, .(draw, dw)]




############################################################################
#             BRING IN HIV-SPECIFIC GBD DRAWS PRE-CENTRAL COMP             #
############################################################################

# Read in the draw-files and compile
inc_draws_by_hiv <- rbind(compile_detailed_covariate_draws('inc_hiv', draw_sample = draw_sample, hiv_level = 1),
                          compile_detailed_covariate_draws('inc_nohiv', draw_sample = draw_sample, hiv_level = 0))

death_draws_by_hiv <- rbind(compile_detailed_covariate_draws('death_hiv', draw_sample = draw_sample, hiv_level = 1),
                             compile_detailed_covariate_draws('death', draw_sample = draw_sample, hiv_level = 0))

inf_sev_draws <- compile_detailed_covariate_draws(var = 'inf_sev', draw_sample = draw_sample)


# If we need to aggregate subnationals to the national estimate do so here
if (length(detailed_loc_ids)>1) {
  inf_sev_draws     <- agg_subnat(inf_sev_draws, agg_var = c('incidence', 'prevalence'))
  inc_draws_by_hiv   <- agg_subnat(inc_draws_by_hiv, agg_var = 'incidence')
  death_draws_by_hiv <- agg_subnat(death_draws_by_hiv, agg_var = 'death')
}

# Back calculate HIV PAFs from incidence draws
inc_draws_by_hiv[, pr_total_by_hiv := incidence / sum(incidence), 
                 by =  c('location_id', 'age_group_id', 'sex_id', 'year_id', 'draw')]

inc_draws_by_hiv[, logit_hiv_paf := logit(pr_total_by_hiv)]


# Merge in HIV PAFs to incidence and prevalence draws and calculate draws of incidence and prevalence by HIV status
inf_sev_draws_by_hiv  <- merge(inf_sev_draws, copy(inc_draws_by_hiv)[, incidence := NULL], 
                               by =  c('location_id', 'age_group_id', 'sex_id', 'year_id', 'draw'), all.x = T)

inf_sev_draws_by_hiv[, c('incidence', 'prevalence') := lapply(.SD, function(x) x * pr_total_by_hiv), .SDcols = c('incidence', 'prevalence')]


# Merge morality and incidence draws by HIV status to back calculate CFR by HIV
cfr_draws_by_hiv   <- merge(death_draws_by_hiv, inc_draws_by_hiv, 
                            by = c('location_id', 'age_group_id', 'sex_id', 'year_id', 'hiv', 'draw'), all = T)

cfr_draws_by_hiv[, cfr := ifelse(age_group_id==2, 0, death / incidence)][, c('death', 'incidence') := NULL]
cfr_draws_by_hiv[, demog_group := .GRP, by = c('location_id', 'age_group_id', 'sex_id', 'year_id', 'hiv')]


# Where there is a zero incidence draw, cfr will be NaN.  For those rows, resample from other non-missing draws of same demographic
miss <- cfr_draws_by_hiv[is.na(cfr), ][, cfr := NULL]

if (nrow(miss)>0) {
  no_miss   <- cfr_draws_by_hiv[!is.na(cfr), ] 
  miss$cfr <-sapply(1:nrow(miss), function(i) {sample(no_miss[demog_group == miss$demog_group[i], cfr], 1)})
  cfr_draws_by_hiv <- rbind(no_miss, miss)
}




############################################################################
#                        INCIDENCE ESTIMATION                              #
############################################################################

# Convert age group and loc to factors and drop age 2 to match the model
covariate_draws[age_group_id==2, age_group_id := NA] # no predictions for age_group 2; convert to NA to avoid error in predict.scam
covariate_draws[, age_group_id := as.factor(age_group_id)]
covariate_draws[, ihme_loc_id := as.factor(ihme_loc_id)]

# Apply forecast model-based uncertainty
set.seed(seed_list[1])
seed_list <- seed_list[-1]
error_draws <- data.table(draw = 0:MAX_FHS_DRAW, random = rnorm(n = (MAX_FHS_DRAW+1), mean = 0, sd = 1))

covariate_draws <- merge(covariate_draws, error_draws, by = 'draw', all.x = T)


# Make predictions from weak learners
weak_learner_preds <- pblapply(weak_learners, function(model_id) {
  do.call(cbind, predict.scam(readRDS(file.path(xval_dir, paste0('step1_model', model_id, '_fold999_ar1.rds'))), 
               covariate_draws, se.fit = TRUE))})

weak_learner_preds <- data.table(cbind(intercept = 1, do.call(cbind, weak_learner_preds)))
weak_learner_preds <- cbind(weak_learner_preds, covariate_draws$random) 

# Change column names of weak_learner_preds to include model_id
names(weak_learner_preds) <- c('intercept', paste0(c('model_', 'se_'), rep(weak_learners, each = 2)), 'random')

# Calculate the SEs of predictions for meta-model (coefficient weighted mean of weak learner errors)
weak_learner_errors <- copy(weak_learner_preds)[, lapply(.SD, function(x) x * random), .SDcols = paste0('se_', weak_learners)]
covariate_draws$inc_log_se  = rowSums(weak_learner_errors[, Map("*", .SD, meta_mdl_coefs@x[-1])])


# Calculate the predictions from the strong learner
strong_learner_preds <- copy(weak_learner_preds)[, .SD, .SDcols = c('intercept', paste0('model_', weak_learners))]
strong_learner_preds <- rowSums(strong_learner_preds[, Map("*", .SD, meta_mdl_coefs@x)])
covariate_draws$inc_total <- exp(strong_learner_preds)

# Revert age_group_id variable back to integer with no NAs
covariate_draws[, age_group_id := as.integer(as.character(age_group_id))]
covariate_draws[is.na(age_group_id), age_group_id := 2] # convert the NAs back to 2 for age_group 2

# Set HIV status to 3 for all draws (this is HIV+ & HIV- combined)
covariate_draws[, hiv := 3]

# Shift the baseline incidence to match GBD
covariate_draws <- baseline_shift(fhs.data = covariate_draws, gbd.data = inc_draws_cc[, hiv := 3], 
                                  fhs.var = 'inc_total', gbd.var = 'incidence', se.var = 'inc_log_se')

# Age group 2 has no incidence, so set to 0 
covariate_draws[age_group_id==2, inc_total_shifted := 0]


# Retain only needed variables #
full <- copy(covariate_draws)[, .SD, .SDcols = c('location_id', 'year_id', 'age_group_id', 'sex_id', 'draw', 'sdi', 'inc_total_shifted')]



############################################################################
#                                HIV SPLIT                                 #
############################################################################

# Create error draws
set.seed(seed_list[1])
seed_list <- seed_list[-1]
error_draws <- data.table(draw = 0:MAX_FHS_DRAW, random = rnorm(n = (MAX_FHS_DRAW+1), mean = 0, sd = 1), merge = 1)

# Read in the diarrhea SEVs
ln_SEV_diarrhea <- fread(file.path(extraction_dir, 'ln_SEV_diarrhea', 'summaries', paste0(LOC_ISO, '.csv')))[ year_id %in% fhs_years]
hiv_rr <- merge(ln_SEV_diarrhea, age_meta[, .(age_group_id, age_mid)], by = 'age_group_id', all.x = T)
hiv_rr[age_mid<=3, age_mid := 3]
setnames(hiv_rr, c('ln_SEV_diarrhea', 'age_mid'), c('lnSevD', 'ageMid'))


predictions <- predict(hiv_rr_model, newdata = hiv_rr, se.fit = TRUE)
hiv_rr <- cbind(hiv_rr, data.table(ln_hiv_rr_mean = predictions$fit, ln_hiv_rr_se = predictions$se.fit))
hiv_rr <- hiv_rr[, .SD, .SDcols = c('location_id', 'year_id', 'age_group_id', 'sex_id', 'ln_hiv_rr_mean', 'ln_hiv_rr_se')]
hiv_rr[, merge := 1]

# Apply error to create draws
hiv_rr <- merge(hiv_rr, error_draws, by = 'merge', all = T, allow.cartesian = T)
hiv_rr[, ln_hiv_rr := ln_hiv_rr_mean + ln_hiv_rr_se * random]
hiv_rr <- hiv_rr[, c('ln_hiv_rr_mean', 'ln_hiv_rr_se', 'merge', 'random') := NULL]


# Merge in HIV prevalence
hiv_rr <- merge(hiv_rr, hiv_prev_draws, by = c('location_id', 'year_id', 'age_group_id', 'sex_id', 'draw'), all.x = T)

# Calculate PAFs from RRs and prevalence
hiv_rr[, hiv_paf := hiv_prev * (exp(ln_hiv_rr) - 1) / (hiv_prev * (exp(ln_hiv_rr) - 1) + 1)]
hiv_rr[hiv_paf < 0, hiv_paf := 0]
hiv_rr[, logit_nohiv_paf := logit(1 - hiv_paf)]
hiv_rr[, hiv := 0]

# Apply baseline shift
hiv_rr <- baseline_shift(fhs.data = hiv_rr, gbd.data = inc_draws_by_hiv[hiv==0, ], fhs.var = 'logit_nohiv_paf', gbd.var = 'logit_hiv_paf')
hiv_rr[, nohiv_paf_shifted := inv.logit(logit_nohiv_paf_shifted)]
hiv_rr[age_group_id==2 | hiv_prev==0, nohiv_paf_shifted := 1]
hiv_rr[is.na(nohiv_paf_shifted)==T, nohiv_paf_shifted := 1 - hiv_paf]
hiv_rr[, hiv_paf := 1 - nohiv_paf_shifted]

# Retain only needed variables
hiv_rr <- hiv_rr[, .(location_id, year_id, age_group_id, sex_id, draw, hiv_paf)]

# Merge into full data table
full <- merge(full, hiv_rr, by = c('location_id', 'year_id', 'age_group_id', 'sex_id', 'draw'), all = T)

# Apply HIV split
full[, inc_hiv_split_shifted := inc_total_shifted * hiv_paf][, hiv := 1]

full <- rbind(full, copy(full)[, hiv := 0])
full[hiv==0, inc_hiv_split_shifted := inc_total_shifted - inc_hiv_split_shifted]

full[, estPrHiv := hiv]


############################################################################
#                              CASE FATALITY                               #
############################################################################

# Read error draws
error_draws <- data.table(draw = 0:MAX_FHS_DRAW, error = as.numeric(fread(file.path(INTS_DIR, 'inputs', 'pred_cfrDrawDeviations.txt'))[1, ]))

full <- merge(full, error_draws, by = 'draw', all = T)
full <- merge(full, copy(age_meta)[, .(age_group_id, age_mid)], by = 'age_group_id', all = T)

# Reconcile variable names with those used to build model
setnames(full, 'age_mid', 'ageMid')

# Truncate SDI to range of data
full[, sdi := ifelse(sdi<0.2, 0.2, ifelse(sdi>0.9, 0.9, sdi))]

# Predict CFR
cfr_pred <- predict.gam(cfr_model, full, se.fit = TRUE)
full$logit_cfr_pred <- as.numeric(cfr_pred$fit)
full$logit_cfr_pred_se <- as.numeric(cfr_pred$se.fit * dfa)

# Smooth standard errors across age groups to correct discontinuities
full <- full[order(location_id, year_id, sex_id, estPrHiv, draw, ageMid), ]
full[, logit_cfr_pred_se_sm := rollmean(logit_cfr_pred_se, k=5, na.pad = TRUE, partial = TRUE), by = c('location_id', 'year_id', 'sex_id', 'estPrHiv', 'draw')]
full[age_group_id<=4 | age_group_id>=31, logit_cfr_pred_se_sm := logit_cfr_pred_se]

# Calculate CFR draws
full[, cfr_pred := inv.logit(logit_cfr_pred + (logit_cfr_pred_se_sm * error))]

# Clean up 
full[, c('logit_cfr_pred', 'logit_cfr_pred_se', 'logit_cfr_pred_se_sm', 'error', 'sdi', 'estPrHiv') := NULL]

# Apply baseline shift
full <- baseline_shift(fhs.data = full, gbd.data = cfr_draws_by_hiv, fhs.var = 'cfr_pred', gbd.var = 'cfr')

full[age_group_id==2 | cfr_pred_shifted<0, cfr_pred_shifted := 0]
full[cfr_pred_shifted>1, cfr_pred_shifted := 1]




############################################################################
#                                MORTALITY                                 #
############################################################################

# Calculate MR draws as the product of incidence and CFR
# Note: Naming 'shifted' for naming consistency, but no baseline shift needed here
# as all inputs have already been shifted
full[, mr_pred_shifted := inc_hiv_split_shifted * cfr_pred_shifted]
full[age_group_id==2 | inc_hiv_split_shifted==0, mr_pred_shifted := 0]

# Create CodCorrected MR draws via baseline shift
full[, mr_cod_corrected := mr_pred_shifted]
full <- baseline_shift(fhs.data = full, gbd.data = mr_draws_cc, fhs.var = 'mr_cod_corrected', gbd.var = 'mr')
full[age_group_id==2 | hiv==1| inc_hiv_split_shifted==0, mr_cod_corrected_shifted := 0]

if (DX==T) {
  # Check the effect of the baseline shift to CodCorrected estimates
  full %>% filter(year_id==ref_year) %>% ggplot(aes(x = mr_pred_shifted, y = mr_cod_corrected_shifted, color = as.factor(hiv))) + 
    geom_point(alpha = 0.1) + geom_abline() + theme_minimal() + facet_wrap(~age_group_id, scales = 'free')
}



############################################################################
#                             YLL CALCULATION                              #
############################################################################

# Read in the reference life table and merge into full data table
standard_lt <- fread(file.path(INTS_DIR, 'inputs', 'standard_lt.csv'))
full <- merge(full, standard_lt, by = 'age_group_id', all.x = T)

# Calculate YLLs as the product of MR and life expectancy
full[, yll_pred := mr_pred_shifted * life_expectancy]
full[, yll_cod_corrected := mr_cod_corrected_shifted * life_expectancy]



############################################################################
#                             YLD CALCULATION                              #
############################################################################



# Generate draws of duration (used to calculate prevalence from incidence)
duration_draws <- copy(inf_sev_draws_by_hiv)[year_id == ref_year, ]
duration_draws <- duration_draws[, duration := prevalence/incidence][, lapply(.SD, mean, na.rm = T), by = 'draw', .SDcols = 'duration']

# Merge duration and disability weight draws into full data table
full <- merge(full, duration_draws, by = 'draw', all.x = T)
full <- merge(full, dw, by = 'draw', all.x = T)

# Calculate prevalence as the product of incidence and duration
# Since all inputs have been baseline shifted, no need to shift here
full[, prev_pred_shifted := inc_hiv_split_shifted * duration]

# Calculate YLDs as the product of prevalence and disability weight
# Since all inputs have been baseline shifted, no need to shift here
full[, yld_pred := prev_pred_shifted * dw]



############################################################################
#                            DALY CALCULATION                              #
############################################################################

# Calculate DALYs as the sum of YLLs and YLDs
# All inputs have been shifted, so no need to shift here
full[, daly_pred := yld_pred + yll_pred]
full[, daly_cod_corrected := yld_pred + yll_cod_corrected]





############################################################################
#                                CLEAN UP                                  #
############################################################################
# Make a backup copy of the full data table if this is a dianostic run
if (DX==T) {
  bkup <- copy(full)
}

# Retain only needed variables
full <- full[, .(location_id, year_id, age_group_id, sex_id, hiv, draw, inc_hiv_split_shifted, 
                 prev_pred_shifted, mr_pred_shifted, mr_cod_corrected_shifted,
                 yll_pred, yll_cod_corrected, yld_pred, daly_pred, daly_cod_corrected)]

# Clean up variable names
names(full) <- gsub('_pred', '', gsub('_shifted', '', names(full)))
setnames(full, 'inc_hiv_split', 'inc')



############################################################################
#            CALCULATE TOTAL ESTIMATES COMBINING HIV+ & HIV-               #
############################################################################

# We currently have estimates for HIV+ and HIV- separately, but we need to combine them to get total estimates
by_vars <- c('location_id', 'year_id', 'age_group_id', 'sex_id')
measure_vars <- setdiff(names(full), c(by_vars, 'hiv', 'draw'))
full <- rbind(full, copy(full)[, lapply(.SD, sum), by = c(by_vars, 'draw'), .SDcols = measure_vars][, hiv := 2])
full[, hiv := factor(hiv, levels = 0:2, labels = c('No HIV', 'HIV', 'Total'))]




############################################################################
#                 CALCULATE ALL-AGE & BOTH-SEX ESTIMATES                   #
############################################################################

merge_vars   <- c('location_id', 'year_id', 'age_group_id', 'sex_id', 'draw')
measure_vars <- c('inc', 'prev', 'mr', 'mr_cod_corrected', 'yll', 'yll_cod_corrected', 'yld', 'daly', 'daly_cod_corrected')

# Merge in population draws and calculate counts from rates
counts <- merge(full, pop_draws[year_id>=ref_year, ], by = merge_vars, all.x = T)
counts <- counts[, c(measure_vars) := lapply(.SD, function(x) {x*population}), .SDcols = measure_vars]


# Estimate both-sex counts
by_vars <- setdiff(c(merge_vars, 'hiv', 'age_group_weight_value'), 'sex_id')
counts <- rbind(counts,
                copy(counts)[, lapply(.SD, sum), by = by_vars, .SDcols = c(measure_vars, 'population')][, sex_id := 3],
                fill = T)

# Estimate all-age counts
by_vars <- setdiff(c(merge_vars, 'hiv'), 'age_group_id')
counts <- rbind(counts, 
                copy(counts)[, lapply(.SD, sum), by = by_vars, .SDcols = c(measure_vars, 'population')][, age_group_id := 22],
                fill = T)

counts[, `:=` (metric_id = 1, metric_name = 'Number')]


# Convert counts to rates
rates <- copy(counts)[, c(measure_vars) := lapply(.SD, function(x) {x / population}), .SDcols = measure_vars]

# Estimate age-standardized rates
age_std <- copy(rates)[age_group_id!=22, ]

by_vars <- setdiff(c(merge_vars, 'hiv'), 'age_group_id')
age_std[,  lapply(.SD, function(x) sum(age_group_weight_value * x)), by = by_vars, .SDcols = measure_vars][, age_group_id := 27]
rates <- rbind(rates, age_std, fill = T)[, `:=` (metric_id = 3, metric_name = 'Rate')]

# Append counts and rates
full <- rbind(counts, rates)

# Back calculate case fatality rates from mortality and incidence
# Note: we already had CFR estimates, but this allows to get CFR for age & sex aggregates too
full[, cfr := ifelse(inc==0, 0, mr / inc)]
full[, cfr_cod_corrected := ifelse(inc==0, 0, mr_cod_corrected / inc)]


# Calculate summary statistics (means and UIs from draws)
by_vars <- c('location_id', 'year_id', 'age_group_id', 'sex_id', 'metric_id', 'hiv')

means <- full[, lapply(.SD, mean), by = by_vars, .SDcols = measure_vars]
setnames(means, measure_vars, paste0(measure_vars, '_mean_pred'))
lowers <- full[, lapply(.SD, quantile, 0.025), by = by_vars, .SDcols = measure_vars]
setnames(lowers, measure_vars, paste0(measure_vars, '_lower_pred'))
uppers <- full[, lapply(.SD, quantile, 0.975), by = by_vars, .SDcols = measure_vars]
setnames(uppers, measure_vars, paste0(measure_vars, '_upper_pred'))

summary <- merge(means, lowers, by = by_vars, all = T)
summary <- merge(summary, uppers, by = by_vars, all = T)


# Export results (if not a diagnostic run)
if (!DX) {
  for (y in fhs_years) {
    print(paste0('saving ', y))
    write.csv(full[year_id==y, ], paste0(OUT_DIR, '/draws/', LOC_ISO, '_', y, '.csv'), row.names = F)
  }
  write.csv(summary, paste0(OUT_DIR, '/summaries/', LOC_ISO, '.csv'), row.names = F)

# Merge with GBD estimates and make diagnostic plots (if a diagnostic run)    
} else {
  gbd_est  <- get_outputs('cause', cause_id = 959, location_id = loc_id, age_group_id = unique(summary$age_group_id), 
                         year_id = 1990:2019,  sex_id = 1:3, measure_id = 1:6, metric_id = c(1,3), release_id = RELEASE)
  
  gbd_est <- gbd_est[, .(location_id, age_group_id, sex_id, year_id, measure, metric_id, val, lower, upper)]
  setnames(gbd_est, c('val', 'lower', 'upper'), c('_mean', '_lower', '_upper'))
  
  gbd_est <- melt.data.table(gbd_est, id.vars = c('location_id', 'age_group_id', 'sex_id', 'year_id', 'measure', 'metric_id'), 
                            measure.vars = c('_mean', '_lower', '_upper'))
  
  gbd_est[, measure := paste0(measure, as.character(variable))][, variable := NULL]
  gbd_est[is.na(value)==T & age_group_id==2, value := 0]
  
  gbd_est <- dcast.data.table(gbd_est, location_id + age_group_id + sex_id + year_id + metric_id ~ measure, value.var = 'value')
  
  compare <- merge(summary, gbd_est, by = c('location_id', 'age_group_id', 'sex_id', 'year_id', 'metric_id'), all = T)

  
  compare %>% dplyr::filter((is.na(hiv)==T | hiv=='Total') & sex_id==1 & age_group_id==18) %>% 
    ggplot(aes(x = year_id)) + 
    geom_vline(xintercept = 2019, linetype = 2, linewidth = 0.1) +
    geom_ribbon(aes(ymin = incidence_lower, ymax = incidence_upper), alpha = 0.2) +
    geom_ribbon(aes(ymin = inc_lower_pred, ymax = inc_upper_pred), alpha = 0.2) +
    geom_line(aes(y = incidence_mean)) + geom_line(aes(y = inc_mean_pred)) +
    theme_minimal() + facet_wrap(~metric_id, scales = 'free') + 
    scale_y_continuous(trans='log2')
  
  
  compare %>% filter((is.na(hiv)==T | hiv=='Total') & sex_id==3 & metric_id==3) %>% 
    ggplot(aes(x = year_id)) + 
    geom_vline(xintercept = 2019, linetype = 2, linewidth = 0.1) +
    geom_ribbon(aes(ymin = incidence_lower, ymax = incidence_upper), alpha = 0.2) +
    geom_ribbon(aes(ymin = inc_lower_pred, ymax = inc_upper_pred), alpha = 0.2) +
    geom_line(aes(y = incidence_mean)) + geom_line(aes(y = inc_mean_pred)) +
    theme_minimal() + facet_wrap(~age_group_id, scales = 'free') 
  
  compare %>% filter((is.na(hiv)==T | hiv=='Total') & metric_id==3) %>% 
    ggplot(aes(x = year_id, fill = as.factor(sex_id))) + 
    geom_vline(xintercept = 2019, linetype = 2, linewidth = 0.1) +
    geom_ribbon(aes(ymin = incidence_lower, ymax = incidence_upper), alpha = 0.2) +
    geom_ribbon(aes(ymin = inc_lower_pred, ymax = inc_upper_pred), alpha = 0.2) +
    geom_line(aes(y = incidence_mean, color = as.factor(sex_id))) + 
    geom_line(aes(y = inc_mean_pred, color = as.factor(sex_id))) +
    theme_minimal() + facet_wrap(~age_group_id, scales = 'free') 
  
  compare %>% filter((is.na(hiv)==T | hiv=='Total') & sex_id==1 & age_group_id==5) %>% ggplot(aes(x = year_id)) + 
    geom_ribbon(aes(ymin = prevalence_lower, ymax = prevalence_upper), alpha = 0.2) +
    geom_ribbon(aes(ymin = prev_lower_pred, ymax = prev_upper_pred), alpha = 0.2) +
    geom_line(aes(y = prevalence_mean)) + geom_line(aes(y = prev_mean_pred)) +
    geom_vline(xintercept = 2019) + theme_minimal() + facet_wrap(~metric_id, scales = 'free')
  
  compare %>% filter((is.na(hiv)==T | hiv=='Total') & sex_id==3 & metric_id==1) %>% ggplot(aes(x = year_id)) + 
    geom_ribbon(aes(ymin = prevalence_lower, ymax = prevalence_upper), alpha = 0.2) +
    geom_ribbon(aes(ymin = prev_lower_pred, ymax = prev_upper_pred), alpha = 0.2) +
    geom_line(aes(y = prevalence_mean)) + geom_line(aes(y = prev_mean_pred)) +
    geom_vline(xintercept = 2019) + theme_minimal() + facet_wrap(~age_group_id, scales = 'free')
  
  compare %>% filter((is.na(hiv)==T | hiv=='No HIV') & sex_id==1 & age_group_id==5) %>% ggplot(aes(x = year_id)) + 
    geom_ribbon(aes(ymin = yll_lower, ymax = yll_upper), alpha = 0.2) +
    geom_ribbon(aes(ymin = yll_cod_corrected_lower_pred, ymax = yll_cod_corrected_upper_pred), alpha = 0.2) +
    geom_line(aes(y = yll_mean)) + geom_line(aes(y = yll_cod_corrected_mean_pred)) +
    geom_vline(xintercept = 2019) + theme_minimal() + facet_wrap(~metric_id, scales = 'free')
  
  compare %>% filter((is.na(hiv)==T | hiv=='Total') & sex_id==1 & age_group_id==5) %>% ggplot(aes(x = year_id)) + 
    geom_ribbon(aes(ymin = death_lower, ymax = death_upper), alpha = 0.2) +
    geom_ribbon(aes(ymin = mr_cod_corrected_lower_pred, ymax = mr_cod_corrected_upper_pred), alpha = 0.2) +
    geom_line(aes(y = death_mean)) + geom_line(aes(y = mr_cod_corrected_mean_pred)) +
    geom_vline(xintercept = 2019) + theme_minimal() + facet_wrap(~metric_id, scales = 'free')
  
  compare %>% filter((is.na(hiv)==T | hiv=='Total') & sex_id==1 & age_group_id==5) %>% ggplot(aes(x = year_id)) + 
    geom_ribbon(aes(ymin = daly_lower, ymax = daly_upper), alpha = 0.2) +
    geom_ribbon(aes(ymin = daly_cod_corrected_lower_pred, ymax = daly_cod_corrected_upper_pred), alpha = 0.2) +
    geom_line(aes(y = daly_mean)) + geom_line(aes(y = daly_cod_corrected_mean_pred)) +
    geom_vline(xintercept = 2019) + theme_minimal()  + facet_wrap(~metric_id, scales = 'free')
  
  
  compare %>% filter((is.na(hiv)==T | hiv=='Total') & sex_id==3 & metric_id==3) %>% ggplot(aes(x = year_id)) + 
    geom_vline(xintercept = 2019, linetype = 2, linewidth = 0.1) +
    geom_ribbon(aes(ymin = daly_lower, ymax = daly_upper), alpha = 0.2) +
    geom_ribbon(aes(ymin = daly_cod_corrected_lower_pred, ymax = daly_cod_corrected_upper_pred), alpha = 0.2) +
    geom_line(aes(y = daly_mean)) + geom_line(aes(y = daly_cod_corrected_mean_pred)) +
    theme_minimal() + facet_wrap(~age_group_id, scales = 'free')

  compare %>% filter((is.na(hiv)==T | hiv=='Total') & sex_id==3 & metric_id==1) %>% ggplot(aes(x = year_id)) + 
    geom_vline(xintercept = 2019, linetype = 2, linewidth = 0.1) +
    geom_ribbon(aes(ymin = death_lower, ymax = death_upper), alpha = 0.2) +
    geom_ribbon(aes(ymin = mr_cod_corrected_lower_pred, ymax = mr_cod_corrected_upper_pred), alpha = 0.2) +
    geom_line(aes(y = death_mean)) + geom_line(aes(y = mr_cod_corrected_mean_pred)) +
    theme_minimal() + facet_wrap(~age_group_id, scales = 'free')
  
  compare %>% filter((is.na(hiv)==T | hiv=='Total') & sex_id==3 & metric_id==1) %>% ggplot(aes(x = year_id)) + 
    geom_ribbon(aes(ymin = daly_lower, ymax = daly_upper), alpha = 0.2) +
    geom_ribbon(aes(ymin = daly_cod_corrected_lower_pred, ymax = daly_cod_corrected_upper_pred), alpha = 0.2) +
    geom_line(aes(y = daly_mean)) + geom_line(aes(y = daly_cod_corrected_mean_pred)) +
    geom_vline(xintercept = 2019) + theme_minimal() + facet_wrap(~age_group_id, scales = 'free')
  
  }

