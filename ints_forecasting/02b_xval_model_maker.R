# Project: iNTS Forecasting
#
# Purpose: 
# This script defines a function that creates the models to be assessed in cross-validation.
#
#
# Author: Jeff Stanaway (stanaway@uw.edu)
# Date created: April 2023
# Date last modified: 25 March 2024
#
# Libraries: data.table


library(data.table)

make_models <- function(xval_dir) {

  dir.create(xval_dir, recursive = TRUE, showWarnings = FALSE)
  
  ## GROUP 1 VARIABLES: malaria, hiv, and their tensor product
  # At least one group 1 variable is included in every model
  group1 <-  c("s(malaria_avg5_quart_root, bs = 'mpi') + ", "s(hiv_quart_root, bs = 'mpi') + ", 
               "s(malaria_avg5_quart_root, bs = 'mpi') + s(hiv_quart_root, bs = 'mpi') + ", 
               "s(malaria_avg5_quart_root, hiv_quart_root, bs = 'tedmi') + ") 
  
  
  
  ## GROUP 2 VARIABLES: water, sanitation, hygiene, sdi, and their principal components;
  # a given model may contain any combination of either the raw variables OR the principal components
  # but no model may contain both raw and principal component variables
  group2a_pc <- c("s(pc1_wash_sdi, bs = 'mpd')", "s(pc1_wash_sdi, bs = 'mpd', by = age_group_id)")
  group2b_pc <- c("s(pc2_wash_sdi, bs = 'cv')", "s(pc3_wash_sdi, bs = 'cv')", "s(pc4_wash_sdi, bs = 'cv')")
  
  # A model may only include the 2nd principal component if it includes the 1st principal component, and
  # it may only include the 3rd principal component if it includes the 2nd principal component, and so on
  group2_pc <- unlist(lapply(1:length(group2a_pc), function(a) 
    sapply(1:length(group2b_pc), function(b) paste(c(group2a_pc[a], paste(group2b_pc[1:b], collapse = " + ")), collapse = " + "))))
  
  # Add an empty string to the beginning of the group2_pc vector to allow for models that include no principal components
  group2_pc <- c("", paste(c(group2a_pc, group2_pc), " + "))
  
  # Make all combinations of the raw variables, including an empty string to allow for models that include no raw variables
  group2_raw <- c("s(water, bs = 'mpi')", "s(sanitation, bs = 'mpi')", "s(hygiene, bs = 'mpi')", "s(sdi, bs = 'mpd')")
  group2_raw <- c("", paste0(unlist(lapply(1:length(group2_raw), function(n) apply(combn(group2_raw, n), 2, paste, collapse = " + "))), " + "))
  
  
  
  ## GROUP 3 VARIABLES: undernutrition, stunting, wasting and their principal components
  # a given model may contain any combination of either the raw variables OR the principal components
  # but no model may contain both raw and principal component variables
  group3_pc <- c("s(pc1_undernutrition_age_std, bs = 'mpi')", "s(pc2_undernutrition_age_std, bs = 'cv')", "s(pc3_undernutrition_age_std, bs = 'cv')")
  group3_pc <-sapply(1:length(group3_pc), function(i) paste(group3_pc[1:i], collapse = " + "))
  group3_pc <- c("", paste(group3_pc, " + "))
  
  group3_raw <- c("s(underweight_age_std, bs = 'mpi')", "s(stunting_age_std, bs = 'mpi')", "s(wasting_age_std, bs = 'mpi')")
  group3_raw <- c("", paste0(unlist(lapply(1:length(group3_raw), function(n) apply(combn(group3_raw, n), 2, paste, collapse = " + "))), " + "))
  
  
  
  ## FIXED VARIABLES: location, age, and sex and are included in all models
  fixed <- "sex_id + age_group_id + s(ihme_loc_id, bs = 're')"
  
  
  ## COMBINE GROUPS 1, 2, AND 3 TO MAKE ALL MODELS
  models <- rbindlist(list(expand.grid(group1, group2_pc, group3_pc, stringsAsFactors = F),
                           expand.grid(group1, group2_raw, group3_raw, stringsAsFactors = F),
                           expand.grid(group1, group2_pc, group3_raw, stringsAsFactors = F),
                           expand.grid(group1, group2_raw, group3_pc, stringsAsFactors = F)), 
                      fill = T)
  
  ## Add a column with the outcome variable and collapse all the columns into a single string
  #  that contains the full model formula as a character string
  models <- cbind('iNTS_log ~ ', models, fixed)
  step1_models <- data.frame(model = apply(models, 1, paste0, collapse = ""))
  
  # Save all the model formulae to an RDS to be read in during cross-validation
  saveRDS(step1_models, file.path(xval_dir, 'step1_models.rds'))
  
  # Save the details of the models to an RDS for use checking details of best fitting models
  step1_model_details <- cbind(step1_models, models)
  names(step1_model_details) <- c('model', 'outcome', 'malaria_hiv', 'development', 'undernutrition', 'fixed')
  
  saveRDS(step1_model_details, file.path(xval_dir, 'step1_model_details.rds'))
  
  return(step1_models)
}