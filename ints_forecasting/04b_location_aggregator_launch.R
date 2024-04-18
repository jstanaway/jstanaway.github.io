# Project: iNTS Forecasting
#
# Purpose: 
# This script launches jobs to make draw-level location aggregates of incidence, prevalence, mortality,
# YLLs, YLDs, and DALYs for the iNTS forecasting project. It moves up levels of the location hierarchy
# in turn, first aggregating any subnational estimates to the national level, then aggregating national
# estimates to the region-level, then aggregating region estimates to the super-region level, and 
# finally aggregating super-region estimates to the global level. It monitors progress of each step
# and automatically launches the next step when the previous step is complete.
#
# Author: Jeff Stanaway (stanaway@uw.edu)
# Date created: February 2023
# Date last modified: 29 March 2024
#
# Libraries: data.table
# User inputs: RELEASE, XVAL_VERSION, FORECAST_VERSION, INTS_SCRATCH_DIR, SHARED_FUN_DIR, POP_DIR,
#              slots, mem (you'll likely not need to change the last two)

#clear memory
rm(list=ls())

# Set the release, version, and directories
RELEASE <- 6
XVAL_VERSION <- 3
FORECAST_VERSION <- 'enet_meta_20240410'

INTS_SCRATCH_DIR <- 'FILENAME_REMOVED_FOR_SECURITY'
SHARED_FUN_DIR <- 'FILENAME_REMOVED_FOR_SECURITY'
POP_DIR  <- 'FILENAME_REMOVED_FOR_SECURITY'


# Load libraries and functions
library('data.table')

source(file.path(SHARED_FUN_DIR, 'get_demographics.R'))
source(file.path(SHARED_FUN_DIR, 'get_location_metadata.R'))
source('~/R/submitExperiment.R')


# Find the reference year based on the release
ref_year <- max(get_demographics('epi', release_id = RELEASE)$year_id)

# Set the directories
version_dir <- file.path(INTS_SCRATCH_DIR, paste0('gbd', ref_year), 'forecasting', XVAL_VERSION, FORECAST_VERSION)
draw_dir <- paste0(version_dir, '/draws/')
out_dir  <- paste0(version_dir, '/summaries/')


# Get the location metadata
loc_meta <- get_location_metadata(location_set_id = 35, release_id = RELEASE)  

# Get the list of locations (based on locations for which population files exist)
# Note, we're not pulling this from loc_meta because FHS uses a different location set
loc_list     <- gsub('.csv', '', list.files(POP_DIR, pattern = '.csv'))
subnats     <- loc_list[nchar(loc_list)>3]
parent_list3 <- unique(substr(subnats, 1, 3))  



# Set up the launch parameters
user <- Sys.getenv('USER')
project  <- 'proj_erf'

sge.output.dir <- paste0(' -o FILENAME_REMOVED_FOR_SECURITY', user, '/output/%x.o%j.out -e FILENAME_REMOVED_FOR_SECURITY', user, '/errors/%x.e%j.err ')
r.shell <- 'FILENAME_REMOVED_FOR_SECURITY/execRscript.sh'
save.script <- 'FILENAME_REMOVED_FOR_SECURITY/04b_location_aggregator.R'

mem <- '50G'
slots <- 8

name.stub <- 'ints'
name.args <- c('parent')

# Launch jobs to aggregate subnational estimates to national level
arg_list3  <- list(parent = parent_list3, release_id = RELEASE, version_dir = version_dir, pop_dir = POP_DIR)
arg_table3 <- launch.jobs(arg.list = arg_list3, name.stub = name.stub, name.args = name.args, outfile.suffix = '.csv',
                          outfile.dir = out_dir, pause = 1, max.attempts = 1, return = T, relaunch = T, queue = 'all.q')

# Launch jobs to aggregate national estimates to region level
arg_list2  <- list(parent = loc_meta[level==2, ihme_loc_id], release_id = RELEASE, version_dir = version_dir, pop_dir = POP_DIR)
arg_table2 <- launch.jobs(arg.list = arg_list2, name.stub = name.stub, name.args = name.args, outfile.suffix = '.csv',
                          outfile.dir = out_dir, pause = 1, max.attempts = 1, return = T, relaunch = T, queue = 'all.q')

# Launch jobs to aggregate region estimates to super-region level
arg_list1  <- list(parent = loc_meta[level==1, ihme_loc_id], release_id = RELEASE, version_dir = version_dir, pop_dir = POP_DIR)
arg_table1 <- launch.jobs(arg.list = arg_list1, name.stub = name.stub, name.args = name.args, outfile.suffix = '.csv',
                          outfile.dir = out_dir, pause = 1, max.attempts = 1, return = T, relaunch = T, queue = 'all.q')

# Launch jobs to aggregate super-region estimates to global level
arg_list0  <- list(parent = loc_meta[level==0, ihme_loc_id], release_id = RELEASE, version_dir = version_dir, pop_dir = POP_DIR)
arg_table0 <- launch.jobs(arg.list = arg_list0, name.stub = name.stub, name.args = name.args, outfile.suffix = '.csv',
                          outfile.dir = out_dir, pause = 1, max.attempts = 1, return = T, relaunch = T, queue = 'all.q')



# Check status of jobs
status <- data.table(location_id = loc_list)
status[, status := file.exists(file.path(out_dir, paste0(loc, '.csv')))]

table(status$status)



