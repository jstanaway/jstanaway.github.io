#!/usr/bin/env python
# coding: utf-8

import getpass
import numpy as np
import os
import pandas as pd
import polars as pl
import pymysql 
import re
import shutil
import subprocess
import time

from datetime import datetime
from db_queries import get_demographics, get_ids, get_location_metadata, get_population, get_sequela_metadata
from save_results import save_results_epi, save_results_cod

RELEASE = 16
LEVEL_3 = 'ints' # set to either 'ints' or 'intest'
CLEAR_OUTPUT_DIR = False
BEST = True

bundles = {'intest': 556, 'ints':3785}
xwalks = {'intest': 17891, 'ints': 7298}
    
input_dir = f'FILEPATH_REMOVED_FOR_SECURITY/{LEVEL_3}/release_{RELEASE}/inputs/'
output_dir = f'FILEPATH_REMOVED_FOR_SECURITY/{LEVEL_3}/release_{RELEASE}/draws/'

# Make a list of all of the cause and modelable entity ids for which we'll produce output
seq_meta = get_sequela_metadata(sequela_set_id=2, release_id=RELEASE)

if LEVEL_3 == 'intest':
    cid_list = [319, 320]
    seq_meta = seq_meta.loc[seq_meta.cause_id.isin(cid_list), 'modelable_entity_id'].tolist()
    upload_ids = cid_list + seq_meta
    output_ids = upload_ids + [2523, 23991, 23992]
elif LEVEL_3 == 'ints':
    cid_list = [959]
    seq_meta = seq_meta.loc[seq_meta.cause_id.isin(cid_list), 'modelable_entity_id'].tolist()
    upload_ids = cid_list + seq_meta + [27540, 28000]
    intermediate_ids = [9999, 9959, 196800, 196801]
    output_ids = upload_ids + intermediate_ids
else:
    print('Value of LEVEL_3 must be either "intest" or "ints".')

for output_id in output_ids + ['model_info']:
    if CLEAR_OUTPUT_DIR == True:
        shutil.rmtree(os.path.join(output_dir, str(output_id)))
    os.makedirs(os.path.join(output_dir, str(output_id)), exist_ok = True)
    
    
cod_demog = get_demographics(gbd_team="cod", release_id = RELEASE)
locs = cod_demog['location_id']

if LEVEL_3 == 'intest':
    loc_meta = get_location_metadata(35, release_id = RELEASE)
    ind_locs = loc_meta[loc_meta.path_to_top_parent.str.contains(',163,')]['location_id'].tolist()

            
user = getpass.getuser()



launch_time = datetime.now()

for loc in locs:
    if LEVEL_3 == 'intest' and loc in ind_locs:
        mem = '100G'
    else:
        mem = '20G'

    submission_list = ['sbatch', '-J', f'{LEVEL_3}_{loc}', '-e', f'/FILEPATH_REMOVED_FOR_SECURITY/{user}/errors/%x.e%j.txt',
                       '-o', f'FILEPATH_REMOVED_FOR_SECURITY/{user}/output/%x.o%j.txt', '-A', 'proj_erf',
                       f'--mem={mem}', '-c', '4', '-t', '600', '-p', 'all.q', 
                       'FILEPATH_REMOVED_FOR_SECURITY/py_shell.sh ' +  'FILEPATH_REMOVED_FOR_SECURITY/enteric_split.py ' + str(loc) + ' ' + str(RELEASE) + ' ' + str(LEVEL_3)]
    
    submission_str = " ".join(submission_list)
    os.system(submission_str)
    
while True:
    time.sleep(30)
    
    all_jobs = subprocess.check_output('squeue --me -o "%j"', shell = True)
    all_jobs = all_jobs.decode().split("\n")[1:]

    matching_jobs = [job for job in all_jobs if re.match(f"^{LEVEL_3}", job)]
    if len(matching_jobs) == 0:
        print('All jobs done running. Checking output files.')
        break
    else:
        print('Jobs still running.  Will check again in a minute.')
        



def file_checker(file, launch):
    exists = os.path.isfile(file)
    if exists: 
        mtime = datetime.fromtimestamp(os.path.getmtime(file))
        if mtime > launch_time:
            return 'new'
        else:
            return 'old'
    else:
        return 'missing'

    
checks = [[id, loc, file_checker(os.path.join(output_dir, str(id), f"{loc}.csv"), launch_time)] for id in output_ids for loc in locs]
checks = pd.DataFrame(checks, columns = ['meid', 'location_id', 'status'])
checks['complete'] = checks['status'] == 'new'

ready = checks.groupby(['meid'])['complete'].mean().reset_index()

meids_to_upload = ready.loc[ready.complete==1, 'meid']
meids_to_upload = [x for x in meids_to_upload if x in upload_ids]
print(checks.groupby('meid')['status'].value_counts())



print(meids_to_upload)
print(upload_ids)
missing_upload_ids = set(upload_ids) - set(meids_to_upload)
if len(missing_upload_ids) == 0:
    print("Results are complete for all upload IDs")
else:
    print("The following upload ids do not have complete results: " + str(missing_upload_ids))



model_info = (pl.scan_csv(os.path.join(output_dir, 'model_info', '*.csv')).collect(streaming = True)).unique()
label = []
for tool, group in model_info.group_by('tool'):
    label.append(tool + ' = ' + ', '.join([str(id) for id in group['model_version_id']]))

label = '; and '.join(label)
description = f'Natural hx / CODEm hybrid using {label}, with python pipeline'
print(description)




for id in meids_to_upload:
    if id in cid_list:
        type = 'cod'
        mem = '200G'
        measures = 1

    else:
        type = 'epi'
        mem = '100G'

        m_test = pl.read_csv(os.path.join(output_dir, str(id), '161.csv'))
        measures = ' '.join(map(str, m_test['measure_id'].unique()))
'

    submission_list = ['sbatch', '-J', f'upload_{id}', '-e', f'FILEPATH_REMOVED_FOR_SECURITY/{user}/errors/%x.e%j.txt',
                       '-o', f'FILEPATH_REMOVED_FOR_SECURITY/{user}/output/%x.o%j.txt', '-A', 'proj_erf',
                       f'--mem={mem}', '-c', '4', '-t', '1000', '-p', 'all.q', 
                       'FILEPATH_REMOVED_FOR_SECURITY/py_shell.sh ' +  'FILEPATH_REMOVED_FOR_SECURITY/bulk_uploader.py ' + 
                       f'--type {type} --id {id} --path {os.path.join(output_dir, str(id))} --description "{description}" ' +
                       f'--measure {measures} --best {BEST} --release {RELEASE} --bundle {bundles[LEVEL_3]} --xwalk {xwalks[LEVEL_3]}']

    submission_str = " ".join(submission_list)
    print(submission_str)
    os.system(submission_str)


