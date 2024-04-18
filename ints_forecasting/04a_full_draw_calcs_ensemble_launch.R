# Project: iNTS Forecasting
#
# Purpose: 
# This script launches jobs to make draw-level predictions of incidence, prevalence, mortality,
# YLLs, YLDs, and DALYs for the iNTS forecasting project. It's a pretty standard script for launching
# jobs on the cluster, and parallelizes the work across locations. Most of the jobs will run in
# under an hour, however countries with subnational locations that are modelled in GBD but not
# FHS will take longer.
#
# Author: Jeff Stanaway (stanaway@uw.edu)
# Date created: February 2023
# Date last modified: 29 March 2024
#
# Libraries: data.table
# User inputs: RELEASE, XVAL_VERSION, FORECAST_VERSION, EXTRACTION_DIR, INTS_SCRATCH_DIR

# Clear the workspace
rm(list=ls())

# Load libraries and functions
SHARED_FUN_DIR <- 'FILEPATH_REMOVED_FOR_SECURITY'

library('data.table')
source(file.path(SHARED_FUN_DIR, 'get_location_metadata.R'))
source(file.path(SHARED_FUN_DIR, 'get_demographics.R'))
source('FILEPATH_REMOVED_FOR_SECURITY/submitExperiment.R')

# Set user inputs
RELEASE <- 6
XVAL_VERSION <- 3
FORECAST_VERSION <- 'enet_meta_20240410'
EXTRACTION_DIR <- 'FILEPATH_REMOVED_FOR_SECURITY'
INTS_SCRATCH_DIR <- 'FILEPATH_REMOVED_FOR_SECURITY'


# Determine reference year based on the release
ref_year <- max(get_demographics('epi', release_id = RELEASE)$year_id)

# Set up directories
version_dir <- file.path(INTS_SCRATCH_DIR, paste0('gbd', ref_year), 'forecasting', XVAL_VERSION, FORECAST_VERSION)
draw_dir <- paste0(version_dir, '/draws/')
out_dir  <- paste0(version_dir, '/summaries/')

dir.create(draw_dir, recursive = T, showWarnings = F)
dir.create(out_dir,  showWarnings = F)

# Set up arguments for the jobs
# Need more memory for GBR since it has so many subnationals, so we'll run it separately
loc_list_hi_mem <- 'GBR'
arg.list_hi_mem <- list(loc = loc_list_hi_mem, release_id = RELEASE, version_dir = version_dir, xval_version = XVAL_VERSION)

# All other locations
loc_list <- gsub('.csv', '', list.files(file.path(EXTRACTION_DIR, 'sdi'), pattern = '.csv'))
loc_list <- setdiff(loc_list, loc_list_hi_mem)
arg.list <- list(loc = loc_list, release_id = RELEASE, version_dir = version_dir, xval_version = XVAL_VERSION)

name.stub <- 'ints'
name.args <- c('loc')



# set up run environment
project  <- 'proj_erf'
sge.output.dir <- paste0(' -o FILEPATH_REMOVED_FOR_SECURITY', user, '/output/%x.o%j.out -e FILEPATH_REMOVED_FOR_SECURITY', user, '/errors/%x.e%j.err ')
r.shell <- 'FILEPATH_REMOVED_FOR_SECURITY/execRscript.sh'
save.script <- 'FILEPATH_REMOVED_FOR_SECURITY/full_draw_calcs_ensemble.R'
slots <- 8

user <- Sys.getenv('USER')

# Launch high memory jobs
mem <- '200G'
arg.table_hi_mem <- launch.jobs(arg.list = arg.list_hi_mem, name.stub = name.stub, name.args = name.args, outfile.suffix = '.csv',
                         outfile.dir = out_dir, pause = 1, max.attempts = 1, return = T, relaunch = F, queue = 'all.q')

# Launch all other jobs
mem <- '50G'
arg.table <- launch.jobs(arg.list = arg.list, name.stub = name.stub, name.args = name.args, outfile.suffix = '.csv',
                         outfile.dir = out_dir, pause = 1, max.attempts = 1, return = T, relaunch = F, queue = 'all.q')


# Check status of jobs
status <- data.table(location_id = c(loc_list, loc_list_hi_mem))
status[, status := file.exists(file.path(out_dir, paste0(loc, '.csv')))]

table(status$status)










