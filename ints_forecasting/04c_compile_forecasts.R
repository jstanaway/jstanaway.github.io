# Project: iNTS Forecasting
#
# Purpose: 
# This script compiles the location-specific forecasts into a single file, merges in metadata, cleans
# up the variables, and saves the final forecasts to a single file for sharing with external collaborators.
#
# Author: Jeff Stanaway (stanaway@uw.edu)
# Date created: February 2024
# Date last modified: 01 April 2024
#
# Libraries: data.table
# User inputs: RELEASE, XVAL_VERSION, FORECAST_VERSION, ROOT_DIR, SHARED_FUN_DIR, COMPILE_SUMMARY_FILES

rm(list = ls())

# Set run parameters
RELEASE <- 6
XVAL_VERSION <- 3
FORECAST_VERSION <- 'enet_meta_20240410'
ROOT_DIR <- 'FILEPATH_REMOVED_FOR_SECURITY'
SHARED_FUN_DIR <- 'FILEPATH_REMOVED_FOR_SECURITY'
COMPILE_SUMMARY_FILES <- TRUE

# Load libraries and functions
library(data.table)

source(file.path(SHARED_FUN_DIR, 'get_age_metadata.R'))
source(file.path(SHARED_FUN_DIR, 'get_demographics.R'))
source(file.path(SHARED_FUN_DIR, 'get_ids.R'))
source(file.path(SHARED_FUN_DIR, 'get_location_metadata.R'))
source(file.path(SHARED_FUN_DIR, 'get_outputs.R'))

# Find the reference year for the release
ref_year <- max(get_demographics('epi', release_id = RELEASE)$year_id)

# Set up the directories
version_dir <- file.path(root_dir, paste0('gbd', ref_year), 'forecasting', XVAL_VERSION, FORECAST_VERSION)
draw_dir <- paste0(version_dir, '/draws/')
summ_dir  <- paste0(version_dir, '/summaries/')
dir.create(graph_dir, showWarnings = F)


# Load the location metadata
loc_meta <- get_location_metadata(location_set_id = 35, release_id = RELEASE)
loc_meta <- loc_meta[, .(location_id, location_name, location_type, level, region_name, super_region_name)]
setnames(loc_meta, 'level', 'location_level')

# Load the age metadata
age_meta <- get_age_metadata(release_id = RELEASE)
age_meta[order(age_group_years_start), age_index := 1:.N]
age_meta <- age_meta[, .(age_group_id, age_group_name, age_group_alternative_name, age_index)]

# Sort by age start, but put all age and age standardized first
age_meta[, age_sort := 3:(.N+2)] # Start at 3 so that all age and age standardized are first
age_meta <- rbind(age_meta, data.table(age_group_id = c(22, 27), age_group_name = c('All age', 'Age standardized'), 
                                     age_group_alternative_name = c('All age', 'Age standardized'), age_index = 1:2+nrow(age_meta), age_sort = c(1, 2)))

# Load the sex and metric metadata
sex_meta <- get_ids('sex')[sex_id <= 3, ]
metric_meta <- get_ids('metric')[metric_id %in% c(1, 3), ][, metric_description := NULL]



# If we are compiling the summary files, do so (this needs to be run once per version, but takes
# a while, so we save the results to a file for future use and avoid rerunning)

if (COMPILE_SUMMARY_FILES) {
  # Bring in forecast summaries
  files <- list.files(summ_dir, '.csv', full.names = T)
  files <- setdiff(files, paste0(summ_dir, '/compare_summaries_full.csv'))
  summary <- rbindlist(lapply(files, fread), fill = T)

  # Bring in GBD estimate summaries
  gbd_est  <- get_outputs('cause', cause_id = 959, location_id = unique(summary$location_id), age_group_id = unique(summary$age_group_id), 
                         year_id = 1990:ref_year,  sex_id = 1:3, measure_id = 1:6, metric_id = c(1,3), release_id = RELEASE)
  
  gbd_est <- gbd_est[, .(location_id, age_group_id, sex_id, year_id, measure, metric_id, val, lower, upper)]
  setnames(gbd_est, c('val', 'lower', 'upper'), c('_mean', '_lower', '_upper'))
  
  gbd_est <- melt.data.table(gbd_est, id.vars = c('location_id', 'age_group_id', 'sex_id', 'year_id', 'measure', 'metric_id'), 
                            measure.vars = c('_mean', '_lower', '_upper'))
  
  gbd_est[, measure := paste0(measure, as.character(variable))][, variable := NULL]
  gbd_est[is.na(value)==T & age_group_id==2, value := 0]
  
  gbd_est <- dcast.data.table(gbd_est, location_id + age_group_id + sex_id + year_id + metric_id ~ measure, value.var = 'value')
  
  
  # Merge GBD estimates and forecasts
  full <- merge(summary, gbd_est, by = c('location_id', 'age_group_id', 'sex_id', 'year_id', 'metric_id'), all = T)
  
  write.csv(full, paste0(summ_dir, '/compare_summaries_full.csv'), row.names = F)
  

# If we are not compiling the summary files, load the previously compiled file  
} else {
  full <- fread(paste0(summ_dir, '/compare_summaries_full.csv'))
}

# Create a vector of id variables and a vector of prediction variables
id_list   <- c(grep('_id', names(full), value = T), 'hiv')
pred_list <- grep('_pred', names(full), value = T)
pred_list <- pred_list[order(pred_list)]

# Clean up the data and merge in metadata to make it comprehensible to external audiences
clean <- copy(full)[year_id >= ref_year, ][, .SD, .SDcols = c(id_list, pred_list)]

clean <- merge(clean, loc_meta, by = 'location_id', all.x = T)
clean <- merge(clean, sex_meta, by = 'sex_id', all.x = T)
clean <- merge(clean, metric_meta, by = 'metric_id', all.x = T)
clean <- merge(clean, copy(age_meta)[, .(age_group_id, age_sort)], by = 'age_group_id', all.x = T)

clean$age_group_name <- factor(clean$age_group_id, levels = age_meta$age_group_id, labels = age_meta$age_group_name)
clean$age_group_alternative_name <- factor(clean$age_group_id, levels = age_meta$age_group_id, labels = age_meta$age_group_alternative_name)

# Clean up variable names
names(clean) <- gsub('_pred', '', names(clean))
names(clean) <- gsub('_cod_corrected', '_cc', names(clean))
names(clean) <- gsub('mr', 'mort', names(clean))
names(clean) <- gsub('_mean', '', names(clean))

# Order the data
clean <- clean[order(location_id, year_id, age_sort, -sex_id, -hiv, metric_id), ]
clean <- clean[, .(location_id, location_name, location_level, region_name, super_region_name, year_id, age_group_id, age_group_name, age_group_alternative_name,
                   sex_id, sex, hiv, metric_id, metric_name, inc, inc_lower, inc_upper, prev, prev_lower, prev_upper, mort, mort_lower, mort_upper, yld, yld_lower, yld_upper,
                   yll, yll_lower, yll_upper, daly, daly_lower, daly_upper, mort_cc, mort_cc_lower, mort_cc_upper, yll_cc, yll_cc_lower, yll_cc_upper, daly_cc, daly_cc_lower, daly_cc_upper)]

# Save the file as a csv for widespread compatibility
write.csv(clean, file.path(summ_dir, paste0('intsForecasting_all_estimates_', format(Sys.Date(), "%Y%m%d"), '.csv')), row.names = F)
