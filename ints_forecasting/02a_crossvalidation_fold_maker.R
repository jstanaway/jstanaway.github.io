# Project: iNTS Forecasting
#
# Purpose: 
# This script creates the folds to be used in cross-validation.  
# Because we're using a very specific variant of rolling basis cross-validation we need to create the 
# folds manually instead of using a pre-existing function 
#
# Author: Jeff Stanaway (stanaway@uw.edu)
# Date created: April 2023
# Date last modified: 10 October 2023
#
# Libraries: data.table, tidyverse
# Dependencies: get_demographics.R
# User inputs: SHARED_FUN_DIR, INTS_DIR, RELEASE, TRAIN_N, GAP_N, SHIFT

rm(list = ls())

library(data.table)
library(tidyverse)


SHARED_FUN_DIR <- 'FILEPATH_REMOVED_FOR_SECURITY'
INTS_DIR <- 'FILEPATH_REMOVED_FOR_SECURITY'
RELEASE <- 6


# Pull the years for which we have data for the current release
source(file.path(SHARED_FUN_DIR, 'get_demographics.R'))
demog <- get_demographics('epi', release_id = RELEASE)

years = min(demog$year_id):max(demog$year_id)
n_years = length(years)


# Number of years in each training set
TRAIN_N = 15 

# Number of years in the gap between the training and test sets
GAP_N = 5

# Number of years to shift the training and test sets between each fold
SHIFT = 5

# Number of years in the test set
test_n = n_years - TRAIN_N - GAP_N

# Because we're most interested in predictive accuracy to years outside our data
# we're using rolling-basis cross-validation, with a gap between each training and test set.
# To ensure that each year of input data occurs the same number of times in each the train, test, and gap
# years, we're going to effectively wrap sets around the time period, such that, if a training set
# would extend beyond the end of the data, it will instead start at the beginning of the data.
# We could accomplish this several different ways, but the simplest seems to be to 
# repeat a sequence of train, gap, and test sets 3 times, and then shift the start of the sequence by 5 years each time.

# Lets create the repeated gap, train, gap, test sequence 
set_seq <- rep(rev(c(rep('gap', GAP_N), rep('train', TRAIN_N), rep('gap', GAP_N), rep('test', test_n))), times = 3)

# Now move along that sequence by the value of SHIFT, and create a each fold 
sets = data.table(year_id = years, do.call(cbind, lapply(seq(1, n_years+4, SHIFT), function(i) {rev(set_seq[i:(i + n_years - 1)])})))
names(sets) <- gsub('^V', 'fold', names(sets))

# Check to make sure that every year is in each set the same number of times
sets$n_train <- rowSums(sets == 'train')
sets$n_test <- rowSums(sets == 'test')
sets$n_gap <- rowSums(sets == 'gap')

# Melt the data to long format and write to file
long <- melt.data.table(sets, id.vars = 'year_id', measure.vars = grep('^fold', names(sets), value = T),
                        variable.name = 'fold', value.name = 'set')

long[, fold := as.integer(gsub('fold', '', fold))]
write.csv(long, file.path(INTS_DIR, 'inputs', 'xvalidation_folds.csv'), row.names = F)


# Plot the folds
ggplot(long, aes(x = year_id, y = fold, fill = factor(set))) + 
  geom_tile(color = 'black') + theme_minimal() + scale_y_reverse() +
  labs(x = 'Year', y = 'Fold', fill = 'Set') + 
  scale_fill_brewer(palette = 'YlGnBu', labels = c('Gap', 'Test', 'Train'))

ggsave(file.path(INTS_DIR, 'plots', 'xvalidation_folds.png'), width = 6, bg = 'white')

