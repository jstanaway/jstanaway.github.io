#!/usr/bin/env python
# coding: utf-8


import numpy as np
import pandas as pd
import pymysql 
from scipy.special import logit, expit
from db_queries import get_covariate_estimates, get_demographics, get_demographics_template, get_location_metadata, get_model_results


def fill_dismod_trend(MEID, MODEL_VERSION, RELEASE, LOC):
    print(f'..Filling time series for MEID {MEID}, model version {MODEL_VERSION}, release {RELEASE}, location {LOC}')
    # Get epi years -- Need this to know what years to interpolate between and extrapolate from
    epi_years = get_demographics('epi', release_id=RELEASE)['year_id']
    
    # Make template data frame of full demographics for prediction + epi estimates
    template =  make_template(MEID, RELEASE, LOC) # Need to add India (163) to template as we want to deal with it's subnationals differently in the split code 
    
    # Get the metadata with covariate fixed effects from the database (this tells us which covariates were included and gives us the coefficients for making predictions)
    model_data, measure_id = get_model_metadata(MODEL_VERSION)
    
    # Get the covariate metadata so we know if covariates are age/sex-specific (need this to know how to merge everything together)
    if len(model_data) > 0:
        covar_data = get_covariate_metadata(model_data)
    
    # Combine the covariate estimates with the template and calculate the linear combination of fixed effects
    if len(model_data) > 0:
        template = combine_covariate_estimates(template, model_data, covar_data, RELEASE)
    
    # Group the template by age, location, and sex, and fill in the time-series for each
    grouped = template.groupby(['age_group_id', 'location_id', 'sex_id'])
    preds, status = zip(*[fill_time_series(group, epi_years, measure_id) for _, group in grouped])
    
    # Clean up the prediction data frame for return
    preds = pd.concat(preds, ignore_index=True)
    preds = pd.merge(template, preds, on=['age_group_id', 'location_id', 'sex_id', 'year_id'], how = 'outer')
    
    preds['pred'] = preds['mean'].fillna(preds['pred'])
    preds['ref_year'] = preds['ref_year'].fillna(preds['year_id'])
    preds = preds[['location_id', 'year_id', 'ref_year', 'age_group_id', 'sex_id', 'mean', 'pred']]
    
    status = pd.concat(status, ignore_index=True)
    
    return preds, status
    


"""
NOTE: this is currently configured to only accept a single location -- 
consider updating to take a list of locations
This isn't urgent, but potentially useful and more generalizable
"""

def make_template(meid, release, loc):
    # Get template with all combinations of all years, ages, and sexes
    template = get_demographics_template('cod', release_id=release)
    if loc in np.unique(template['location_id']):
        template = template.query(f'location_id == {loc}')
    else:
        tmp_loc = template['location_id'].iloc[0]
        template = template.query(f'location_id == {tmp_loc}').drop(columns = 'location_id')
        template['location_id'] = loc
    
    
    
    # Get the model estimates and keep only necessary variables
    est = get_model_results('epi', meid, release_id=release, 
                            location_id = list(template['location_id'].unique()),
                            sex_id = list(template['sex_id'].unique()), 
                            age_group_id = list(template['age_group_id'].unique()))
    
    est = est[['location_id', 'year_id', 'age_group_id', 'sex_id', 'mean']]
    
    # Merge the template with the model estimates and return
    template = pd.merge(template, est, on=['location_id', 'year_id', 'age_group_id', 'sex_id'], how = 'outer')
    
    return template



def get_model_metadata(MODEL_VERSION):
    # We're going to pull the fixed effects from the DisMod model and use them to extrapolate
    # These are the variables we want to get re: fixed effects
    varlist = ['measure_id', 'country_covariate_id', 'mean_effect', 'upper_effect', 'lower_effect']
    
    # Connect to the epi database, open the cursor, and execute the SQL query
    db = pymysql.connect(host = 'HOST_REMOVED_FOR_SECURITY', user = 'USER_REMOVED_FOR_SECURITY', password = 'PASS_REMOVED_FOR_SECURITY', database = 'DB_REMOVED_FOR_SECURITY') 
    
    with db:
        with db.cursor() as cursor:
            cursor.execute('SELECT ' + ', '.join(varlist) + ' FROM model_effect WHERE model_version_id = ' + str(MODEL_VERSION) + 
                           ' AND (country_covariate_id IS NOT NULL OR study_covariate_id IS NOT NULL)')

            # Fetch all rows of data, put them in a data frame and add column names
            model_data = pd.DataFrame(cursor.fetchall(), columns = varlist)
            measure_id = model_data['measure_id'][0]
            model_data.dropna(subset = ['country_covariate_id'], inplace=True)

    return model_data, measure_id




def get_covariate_metadata(model_data):
    # Now we know what covariates were used in the model, let's get metadata for those covariates.
    # Steps below are essentially identical to the SQL query above
    # Get list of necessary covariates from the model_data 
    covar_str = ', '.join([str(x) for x in model_data['country_covariate_id']])
    covar_varlist = ['covariate_id', 'covariate_name_short', 'by_sex', 'by_age']
    
    # Connect to the shared database, open the cursor, and execute the SQL query
    db = pymysql.connect(host = 'HOST_REMOVED_FOR_SECURITY', user = 'USER_REMOVED_FOR_SECURITY', password = 'PASS_REMOVED_FOR_SECURITY', database = 'DB_REMOVED_FOR_SECURITY') 
    
    with db:
        with db.cursor() as cursor:
            cursor.execute('SELECT ' + ','.join(covar_varlist) + ' FROM covariate WHERE covariate_id IN (' + covar_str + ')')
    
            # Fetch all rows of data, put them in a data frame and add column names
            covar_data = pd.DataFrame(cursor.fetchall(), columns = covar_varlist)
    
    return covar_data



def combine_covariate_estimates(template, model_data, covar_data, release):
    # Now that we have the covariate metadata, we can pull the covariate value.
    # Since the necessary shared function exists, we'll use that below instead of a SQL query
    
    # We're going to add up the linear combination of effects 
    # (ie. coeffcient * covariate value) for each covariate in turn.
    # Intialize lincom variable to zero here.
    template['lincom'] = 0
    
    # Loop through the covariates
    for i in range(len(covar_data)):
        
        # Pull covariate and model metadata from covar_data and model_data, respectively, and move to individual objects
        id, name, by_sex, by_age = list(covar_data.iloc[i])
        effect = float(model_data.loc[model_data.country_covariate_id == id, 'mean_effect'])
    
        # The id variables on which we'll merge will depend on whether or not the covariate is age- and sex-specific.
        # Location and year are the minimum necessary merge variables; 
        # we'll add sex for sex-specific covariates, and age for age-specific covariates below.
        merge_vars = ['location_id', 'year_id']
    
        if by_sex==1:
            merge_vars.append('sex_id')
            
        if by_age==1:
            merge_vars.append('age_group_id')
        
        # Get_covariates returns variables we don't need here; 
        # We only want to keep the merge variables and the mean value of the covariate ('mean_value').
        keep = merge_vars[::]
        keep.append('mean_value')  
        
        # Pull the covariate estimates
        covar = get_covariate_estimates(covariate_id = int(id), location_id = list(template['location_id'].unique()), release_id = release)
        covar = covar[keep]
        
        # Rename 'mean_value' to the short name of the covariate
        covar.rename(columns = {'mean_value': name}, inplace = True)
        
        # Merge the covariate into the template
        template = pd.merge(template, covar, on = merge_vars, how = 'left')
        
        # Calculate linear combination
        template['lincom'] = template['lincom'] + (template[name] * effect)
        
        # Display numbers of missing values as diagnostic:
        # if all is good, I'd expect this to be all zero missing values
        print('....Missing values of ' + name + ': ' + str(template[name].isna().sum()))
        print('....Missing values of pred: ' + str(template['lincom'].isna().sum()))
        
    return template




"""
Here we're going to define the function to actually do the interpolation and extrapolation.
It takes the template data frame and a list of reference years as inputs, 
where reference years are those between which you want to interpolate.

N.B. IMPORTANT!!! the data frame that you feed into this function should contain estimates for only a single
demographic group (i.e. a single location, age, sex combination) with one row per year.  The main function
handles this for you, but take care if you're going to build your own wrapper around this!

It will return a data frame with location, year, age, sex, and the filled time series of predictions

"""

def fill_time_series(df, ref_years, measure_id):
    
    # fts_main is the main function that manages the data and calls the interpolation and extrpolation functions
    def fts_main(df, ref_years, measure_id):
        # Determine year pairs between which we need to interpolate
        year_pairs = zip(ref_years, ref_years[1:])
        year_pairs = [pair for pair in year_pairs if pair[1] - pair[0] != 1]
        
        # Establish variables to be returned by interpolation/extrapolation code
        return_vars = ['age_group_id', 'location_id', 'sex_id', 'year_id', 'ref_year', 'pred']
        
        # Determine if we have a linear effect to use in the predictions; if not fill a lincom variable with zeros
        if 'lincom' not in df.keys():
            #if to_print:
            status_lincom = False    
            df['lincom'] = 0
        else:
            status_lincom = True
    
    
        # Determine the correct transformation for the measure (i.e. log for rates, logit for proportions)
        if measure_id == None:
            xform = None
        elif int(measure_id) in [5, 17, 18]:
            xform = 'logit'
        else:
            xform = 'log'
        
        # Backcast the estimates
        back = covar_extrapolate(df, min(ref_years), return_vars, xform)

        # Interpolate the estimates
        inter = []
        for year_pair in year_pairs:
            result = covar_interpolate(df, year_pair, return_vars, xform)
            inter.append(result)

        preds = pd.concat([back, pd.concat(inter)])  
        
        status = pd.DataFrame({'lincom': [status_lincom], 'xform': [xform]})
        return preds, status
        
    
    # covar_interpolate will interpolate between two years; it takes the single-demographic template dataframe and list with two years 
    def covar_interpolate(df, year_pair, return_vars, xform = None):
        # Restrict the data to only relevant years and sort by year
        df = df.loc[df.year_id.between(*year_pair)]
        df = df.sort_values('year_id')
        
        # Find the log-transformed DisMod estimates for the reference years, and change between them
        if xform=='log':
            start, end = np.log(df['mean'].iloc[[0,-1]]) 
        elif xform=='logit':
            start, end = logit(df['mean'].iloc[[0,-1]]) 
        else:
            start, end = df['mean'].iloc[[0,-1]]
            
        change = end - start
        
        # Get the linear predictions, start and end values, and change between them
        preds = df['lincom']
        pred_start, pred_end = preds.iloc[[0,-1]]
        pred_change = pred_end - pred_start
        
        # Find the difference in the average annual rates of change in model estimates and linear predictions
        # If those two are different then we will get discontinuities in the predicitions between each interpolated period
        # we avoid that by shifting the overall slope of the linear predictions to match the model estimates
        # Think of this as maintaing the shape of the trend line in the linear predictions, fixing one end of the line to match
        # the corresponding model estimate value, and then rotating the trend line until it matches the model estimate at the
        # other end of the line.
        annual_change = (change - pred_change) / (len(preds)-1)
        shift = [float(x)*annual_change for x in range(len(preds))]
        
        # Make the final predictions here
        if xform=='log':
            df['pred'] = np.exp(start - pred_start + preds + shift)
        elif xform=='logit':
            df['pred'] = expit(start - pred_start + preds + shift)
        else:
            df['pred'] = np.exp(start - pred_start + preds + shift)
        
        df['ref_year'] = int(year_pair[0])
        
        return df[return_vars].iloc[:-1]
    
    
    
    def covar_extrapolate(df, ref_year, return_vars, xform = None, direction = 'backwards'): 
        # Determine if we're extrapolating backwards or foward in time and set up accordingly
        if direction=='backwards':
            df = df.loc[df.year_id <= ref_year]
            df = df.sort_values('year_id')
        elif direction=='forward':
            df = df.loc[df.year_id >= ref_year]
            df = df.sort_values('year_id', ascending=False)    
        else:
            print('direction argument must be either "forward" or "backwards"')
            return
        
        # Process here is really similar to interpolation but simplified, as we don't need
        # to shift the slope to match model estimates at both ends
        lincoms = df['lincom']
        lincom_ref = lincoms.iloc[-1]
        
        if xform=='log':
            ref = np.log(df['mean'].iloc[-1])
            df['pred'] = np.exp(lincoms + ref - lincom_ref)
        elif xform=='logit':
            ref = logit(df['mean'].iloc[-1])
            df['pred'] = expit(lincoms + ref - lincom_ref)
        else:
            ref = df['mean'].iloc[-1]
            df['pred'] = lincoms + ref - lincom_ref        
       
        df['ref_year'] = int(ref_year)
        
        if direction=='backwards':
            return df[return_vars].iloc[:-1]
        else:
            df = df.sort_values('year_id')
            return df[return_vars]
    
    return fts_main(df, ref_years, measure_id)

