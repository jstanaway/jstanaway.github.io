# Project: iNTS Forecasting
#
# Purpose: 
# This script aggregates draws from child locations to a parent location (i.e. aggregating
# country-level estimates to the region-level estimate).  It is a worker script that is launched
# by 04b_location_aggregator.R, and parallelized by location (i.e. one job for each parent location).
# It aggregates in count space, then calculates rates (including age-standardized) and case-fatality, 
# then saves both a draw file and summary file.
#
# Author: Jeff Stanaway (stanaway@uw.edu)
# Date created: April 2023
# Date last modified: 29 March 2024
#
# Libraries: data.table
# User inputs: No hard coded inputs -- everything is passed as arguments from the launch script

# Clear the workspace
rm(list = ls())

# Load libraries and functions
SHARED_FUN_DIR <- 'FILEPATH_REMOVED_FOR_SECURITY'

library('data.table')
source(file.path(SHARED_FUN_DIR, 'get_location_metadata.R'))
source(file.path(SHARED_FUN_DIR, 'get_demographics.R'))


# Determine if this is an interactive run (i.e. for development and debugging) or not
INTERACTIVE = commandArgs()[2]=='--interactive'

if (!INTERACTIVE) {  
  arg     <- commandArgs(trailingOnly = T)
  parent  <- arg[1]
  release <- arg[2]
  version_dir <- arg[3]
  pop_dir <- arg[4]
  
} else { # running interactively (need to hard code inputs for testing)
  parent  <- 'BRA'
  release <- 6
  version_dir <- 'FILEPATH_REMOVED_FOR_SECURITY'
  pop_dir <- 'FILEPATH_REMOVED_FOR_SECURITY'
}




# Build subdirectories
draw_dir <- file.path(version_dir, 'draws')
summ_dir <- file.path(version_dir, 'summaries')


# Get the reference year based on the release
ref_year <- max(get_demographics('epi', release_id = release)$year_id)


# Get the location metadata
loc_meta <- get_location_metadata(location_set_id = 35, release_id = release)  


# Identify the level of the parent location
parent_level <- loc_meta[ihme_loc_id==parent, level]


# Get the children of the parent location
if (parent_level==3) {
  pop_files <- list.files(pop_dir, pattern = '.csv')
  subnats  <- pop_files[nchar(pop_files)>7]
  children <- gsub('.csv', '', grep(parent, subnats, value = T))
} else {
  parent_loc_id <- loc_meta[ihme_loc_id==parent, location_id]
  children <- loc_meta[parent_id==parent_loc_id, ihme_loc_id]
  children <- setdiff(children, parent)
}


# Define the function to aggregate location draws
agg.locs <- function(parent, children, year, save = T) {
  message(paste0('Aggregating ', parent, ', ', year))
  
  merge_vars   <- c('location_id', 'year_id', 'age_group_id', 'sex_id', 'hiv', 'draw')
  measure_vars <- c('inc', 'prev', 'mr', 'mr_cod_corrected', 'yll', 'yll_cod_corrected', 'yld', 'daly', 'daly_cod_corrected')
  
  parent_loc_id <- loc_meta[ihme_loc_id==parent, location_id]
  
  # Read in files with child location draws, retain only counts 
  dt <- do.call(rbind, lapply(file.path(draw_dir, paste0(children, '_', year, '.csv')), fread))[metric_id==1, ]
  
  # Sum draws of counts to get location aggregates (easier to sum counts than to take pop weighted mean of rates)
  dt <- dt[, lapply(.SD, sum), by = setdiff(c(merge_vars, 'age_group_weight_value', 'metric_id', 'metric_name'), 'location_id'), .SDcols = c(measure_vars, 'population')]
  dt[, location_id := parent_loc_id]
  
  # Convert counts to rates
  dt <- dt[, c(measure_vars) := lapply(.SD, function(x) {x/population}), .SDcols = measure_vars]
  
  # Estimate age-standardized rates
  age_std <- copy(dt)[age_group_id!=22, ][,  lapply(.SD, function(x) {sum(age_group_weight_value*x)}), by = setdiff(merge_vars, 'age_group_id'), .SDcols = measure_vars]
  age_std[, age_group_id := 27]
  rates <- rbind(dt, age_std, fill = T)
  
  rates[, `:=` (metric_id = 3, metric_name = 'Rate', population = NULL, age_group_weight_value = NULL)]
  
  
  # Merge in parent population estimates
  pop_draws <- fread(paste0(pop_dir, parent, '.csv'))[year_id==year, ][, .SD, .SDcols = c(setdiff(merge_vars, 'hiv'), 'population', 'age_group_weight_value')]
  rates <- merge(rates, pop_draws, by = setdiff(merge_vars, c('hiv')), all.x = T)

  
  # Calculate parent counts (we use parent pops with aggregated child rates to account for aggregate location population scalars)
  counts <- copy(rates)[age_group_id!=27, ][, c(measure_vars) := lapply(.SD, function(x) {x*population}), .SDcols = measure_vars]
  counts[, `:=` (metric_id = 1, metric_name = 'Number')]
  
  
  # Combine counts and rates
  full <- rbind(counts, rates)
  
  # Recalculate cfr for the aggregate location
  full[, cfr := ifelse(inc==0, 0, mr / inc)]
  full[, cfr_cod_corrected := ifelse(inc==0, 0, mr_cod_corrected / inc)]
  
  # Save the draws
  if (save==T) {
    write.csv(full[year_id==year, ], file.path(draw_dir, paste0(parent, '_', year, '.csv')), row.names = F)
  }
  
  return(full)
}


# Run the aggregation function for each year
start <- Sys.time()
full  <- rbindlist(lapply(ref_year:2100, function(y) agg.locs(parent, children, y, save = T)))
difftime(Sys.time(), start)


# Calculate the mean and 95% CI for each measure
by_vars <- c('location_id', 'year_id', 'age_group_id', 'sex_id', 'hiv', 'metric_id')
measure_vars <- c('inc', 'prev', 'mr', 'mr_cod_corrected', 'yll', 'yll_cod_corrected', 'yld', 'daly', 'daly_cod_corrected')

means <- full[, lapply(.SD, mean), by = by_vars, .SDcols = measure_vars]
  setnames(means, measure_vars, paste0(measure_vars, '_mean_pred'))
lowers <- full[, lapply(.SD, quantile, 0.025), by = by_vars, .SDcols = measure_vars]
  setnames(lowers, measure_vars, paste0(measure_vars, '_lower_pred'))
uppers <- full[, lapply(.SD, quantile, 0.975), by = by_vars, .SDcols = measure_vars]
  setnames(uppers, measure_vars, paste0(measure_vars, '_upper_pred'))

summary <- merge(means,   lowers, by = by_vars, all = T)
summary <- merge(summary, uppers, by = by_vars, all = T)

# Save the summary file
write.csv(summary, file.path(summ_dir, paste0(parent, '.csv')), row.names = F)


