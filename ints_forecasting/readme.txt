~ iNTS forecasting code pipeline ~

This folder contains all code for the iNTS forecasting project (Note if you are reading this on my personal GitHub, then it will only contain code that I personally wrote; scripts developed entirely or in part by others are excluded to ensure that all work presented here is my own). As of April 2024 this work is complete except for a few loose ends, but not yet published.  I presented an earlier version of this analysis to the BMGF and at the 13th International Conference on Typhoid and Other Invasive Salmonelloses (December 2023, Kigali, Rwanda) and made several improvements based on feedback from those presentations.  The most notable changes include the addition of predictors related to undernutrition, changing transformation of some predictors from log to fourth-root (done to avoid needing an offset with zero values which led to some instability in scenario-based forecasts), and moving from making predictions based on on only the best model to making predictions based on an ensemble/stacking approach in which I now fit a meta-model on the predictions of candidate models.  I'm currently drafting the manuscript and will include that in this repo when it's complete.

A few practical notes: this code is all written to run on the IHME computing cluster and requires access to IHME-specific resources (e.g. shared functions and databases) and can not be run externally.  Moreover, as per IHME policy, all filepaths have been removed from the code for security purposes.


The overall process is as follows:
1) Prepare the data
2) Develop pool of weak-learner models
   a) Create folds for crossvalidation
   b) Define all potential models to be tested in crossvalidation
   c) Run step 1 crossvalidation and calculate fit statistics to quantify each candidate model's out-of-sample predictive performance
   d) Remove any models with RMSE >1SE above best model from the pool of candidate models
   e) Run each remaining candidate model on the full dataset
   f) Assess unexplained temporal autocorrelation for each remaining candidate model run on the full dataset
   g) Run each remaining candidate model on the full dataset with the correction for autocorrelation and make predictions
3) Develop strong-learner/meta-model
   a) Run step 2 crossvalidation to choose hyperparameters for the meta-model (strong learner) -- using elastic net here
   b) Choose the best model as that with minimum OOS RMSE and run on the full dataset
4) Create full suite of predictions with draw-based uncertainty
   a) Run draw-based prediction code (parallelized by location)
   b) Calculate aggregate location estimates
   c) Process prediction draws to produce tables and graphs of final estimates
5) Create scenario-based forecasting workspace


The code to run the pipeline closely follows the above outline and most script filenames are prefaced by the step number to make the run order clear (note that step numbers in file names don't perfectly match those in the above outline as some scripts contain multiple steps). To run everything from start to finish, follow these steps...

1) 01_data_prep.R
This script prepares a dataset for both modelling and prediction (without uncertainty, more below) that includes values of all predictive covariates for every combination of age, sex, location, and year for 1990 through 2100, and values of the outcome (iNTS incidence) for all demographics for the years 1990 through 2019.

2) 02a_crossvalidation_fold_maker.R
This script creates the folds to be used in cross-validation.  Because we're using a very specific variant of rolling basis cross-validation we need to create the folds manually instead of using a pre-existing function 

3) 02c_xval_model_run_launch.R
Given the size of the dataset and complexity of the models, we need to parallelize crossvalidation for this analysis to be tractable in a reasonable amount of time.  This script launches the parallelized crossvalidation jobs.  The script that's doing the actual work is xval_model_run.R, but it is built to be robust to changes between cycles should require no changes or maintainance unless we're changing our methods. The launch script can also optionally run 02b_xval_model_maker.R which will create an object (and save corresponding .rds files) containing details of every model variant to be tested in crossvalidation. This will need to be run once for a crossvalidation version but can be disabled for later runs to save a bit of computation. Plan for the crossvalidation run to take a while: the average model will run on a single fold in a few hours, but some may take a few days, and with >8000 model-fold combinations you can plan for this to take a few days when the cluster is reasonably clear or substantially longer when it's busy -- a week is typical. You'll need to wait for all crossvalidation jobs to complete before moving on to the next step.

Calls: 02b_xval_model_maker.R (optionally), 02c_xval_model_run.R

4) 02d_xval_select_best_candidates.R
This script should only be run when all jobs launched in the previous step are complete.  It reads in the fit statistics from every model and fold and assesses models based on errors in OOS predictions in the point estimate (rmse) and the predicted trend (rmse_trend), and retains only those models with RMSE values within 1SE of the best model with regard to both point estimates and trend.  It then launches jobs to run those selected models on the full dataset.

Calls: 02c_xval_model_run.R

5) 02e_xval_ar1_model_launch.R  
This script determines which models launched in the previous step completed and launches a job to run the corresponding model with AR1 errors to account for residual autocorrelation.  This extra step is required because, while SCAM can account for autocorrelation, it cannot estimate the autocorrelation coefficient and requires the user to pass that in as an argument.The models produced here will be the final models used as weak learners to be stacked in the meta-model or strong learner built in the next steps.  Note, this script determines which models are ready to be run (i.e. those in which the full model from the previous step completed), and have not yet been launched (i.e. it determines which ar1 models are either currently running or have already completed).  Therefore, you don't have to wait until all full models are complete to run it, and it can be re-run as needed to launch all relevant models -- useful if you're trying to push things through as quickly as possible and don't want to wait for everything from the previous step to finish before moving on. 

Calls: 02c_xval_model_run.R

6) 03a_build_meta_model.R
This script builds the meta-model based on weak learner predictions. The meta-model here is an elastic net with a penalty factor derived from OOS RMSE, such that the models with smaller OOS errors are given lower penalties (and therefore more likely to be included and with larger coefficients in the meta-model), and models with larger OOS errors are given larger penalties (and therefore less likely to be included or have smaller coefficients in the meta-model.  Note that only those weak-learners meeting the 1SE threshold in the previous step are eligible for inclusion in the meta model. It also launches crossvalidation on the meta-model to determine the combination of hyperparameter values that result in the smallest out-of-sample prediction errors.  There are two relevant hyperparameters: alpha, which governs the mix of L1 and L2 regularization (i.e. the degree to which the model behaves like lasso vs ridge regression), and lambda, which governs the strength of regularization (i.e. the amount of shrinkage).

Calls: 03b_xval_model_run_step2.R


7) 04a_full_draw_calcs_ensemble_launch.R
This script launches location-specific jobs to estimate draw-level forecasts.  The real work in this step is done by the worker script 04a_full_draw_calcs_ensemble.R

Calls: 04a_full_draw_calcs_ensemble.R
This script makes draw-level predictions of iNTS incidence, prevalence, mortality, YLLs, YLDs, and DALYs for all years (up to 2100), age groups, and sexes for a single location.  It is launched by 04a_full_draw_calcs_ensemble_launch.R.  This is a pretty long piece of code as it's doing a lot of work: processing these draw-level predictions requires reading in a large number of inputs with disperate formatting, processing all forecasts at the draw-level to correctly propogate uncertainty, predicting based on an ensemble of models, ensuring continuity of GBD and forecasted estimates (including continuity of uncertainty), and replicating much of the work normally handled by central comp (i.e. gap metric estimation and age/sex aggregation).  Broadly, the structure is to read in and process all inputs; make model-based predictions of incidence, HIV PAFs, and case-fatality; calculating prevalence, mortality, YLLs, YLDs, and DALYs based on those model-based predictions; and calculating age and sex aggregates.

8) 04b_location_aggregator_launch.R
This script launches jobs to make draw-level location aggregates of incidence, prevalence, mortality, YLLs, YLDs, and DALYs for the iNTS forecasting project. It moves up levels of the location hierarchy in turn, first aggregating any subnational estimates to the national level, then aggregating national estimates to the region-level, then aggregating region estimates to the super-region level, and finally aggregating super-region estimates to the global level. It monitors progress of each step and automatically launches the next step when the previous step is complete.

Calls: 04b_location_aggregator.R
 
9) 04c_compile_forecasts.R
This script compiles the location-specific forecasts into a single file, merges in metadata, cleans up the variables, and saves the final forecasts to a single file for sharing with external collaborators.

10) 05a_scenario_prep_ensemble.R
This script builds a workspace that can be shared externally and used to make scenario-based iNTS forecasts. It is a stripped down version of the full draw-level forecasting script, and can be run outside of IHME's computing environment. To make this tractable, this script compiles all necessary functions and data into a single RData file, and the forecasts are run on point estimates rather than draws (i.e. no uncertainty).
