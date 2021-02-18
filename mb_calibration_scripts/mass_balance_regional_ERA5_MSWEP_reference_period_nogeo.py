# Python imports
import json
import os
import sys

# Libs
import geopandas as gpd
import numpy as np
import pandas as pd
import shapely.geometry as shpg

# Locals
import oggm
from oggm import cfg, utils, tasks
from oggm.workflow import execute_entity_task
from oggm.core.massbalance import (ConstantMassBalance, PastMassBalance,
                                   MultipleFlowlineMassBalance)

# Module logger
import logging
log = logging.getLogger(__name__)
from IPython import embed

# RGI Version
rgi_version = '62'

#Specify and RGI region
#rgi_region = '13'

baseline ='custom'

# Initialize OGGM and set up the run parameters
cfg.initialize()
cfg.PARAMS['mu_star_halfperiod']=8.0

# Local paths (where to write the OGGM run output)
WORKING_DIR = '/exports/csce/datastore/geos/groups/geos_iceocean/dgoldber/oggm/oggm_runs/oggm_mb_calibration_ERA5_MSWEP_reference_period_nogeo'
utils.mkdir(WORKING_DIR, reset=True)
cfg.PATHS['working_dir'] = WORKING_DIR

# We are running the calibration ourselves
cfg.PARAMS['run_mb_calibration'] = True

# We are using which baseline data?
cfg.PARAMS['baseline_climate'] = 'CUSTOM'

# change the custom climate data
cfg.PATHS['climate_file'] = '/exports/csce/datastore/geos/groups/geos_iceocean/kinnear/SWARM_files/oggm_SWARM_input_cru_ref_hgts/oggm_cru_hgt_input.nc'
#cfg.PATHS['climate_file'] = '/home/dgoldber/iceOceanShare/dgoldber/oggm/SWARM_inputs/swarm_cru_hgt_input_padmonths.nc'

# Use multiprocessing?
cfg.PARAMS['use_multiprocessing'] = False

# Limit the number of cores, Shin and other geos servers have 64 so limit to half that
cfg.PARAMS['mp_processes'] = 24

# Set to True for operational runs - here we want all glaciers to run
cfg.PARAMS['continue_on_error'] = True

cfg.PARAMS['border'] = 10

cfg.PARAMS['temp_melt'] = -1.25

basin = gpd.read_file('/exports/csce/datastore/geos/groups/geos_iceocean/dgoldber/oggm/SWARM-OGGM-Model/shape_files/swarm_grid.shp')

# Load the glacier ids
#Can dothis either by our own csv file as below or by editing the file downloaded from OGGM
#df = pd.read_csv('/Users/louis/China_test/rids.csv')
#rids = df['RGI{}0_ID'.format(rgi_version[0])]
#Downloading the OGGM data and filtering.
geo_folder_path='/exports/csce/datastore/geos/groups/geos_iceocean/dgoldber/oggm/geodetic_data/'
geo_file_name = 'rgi_wgms_links_with_geo.csv'


#df, _ = utils.get_geodetic_files(geodetic_folder_path=geo_folder_path, geodetic_filename=geo_file_name)
df, _ = utils.get_wgms_files()

#Check this file if needed by printing to a csv (troubleshooting)
#df.to_csv('/exports/csce/datastore/geos/users/s0933963/oggm_mb_test/df.csv')
in_bas = [basin.geometry.contains(shpg.Point(x, y))[0] for
          (x, y) in zip(df.CenLon, df.CenLat)]
df = df.loc[in_bas]
#df = df.loc[df.RGI60_ID == 'RGI60-15.07886']
rids = df['RGI{}0_ID'.format(rgi_version[0])]

#rids.to_csv('/exports/csce/datastore/geos/users/s0933963/oggm_mb_test/rids.csv')

# We have to check which of them actually have enough mb data.
# Let OGGM do it:
from oggm.shop import rgitopo
gdirs = rgitopo.init_glacier_directories_from_rgitopo(rids)
# print(gdirs)

# We need to know which period we have data for
log.info('Process the climate data...')
execute_entity_task(tasks.process_climate_data, gdirs, y1=2018)
#execute_entity_task(tasks.process_climate_data, gdirs)

# Let OGGM decide which of these have enough data
#gdirs = utils.get_ref_mb_glaciers_geodetic(gdirs,temp_geodetic_folder_path=geo_folder_path,temp_geodetic_filename=geo_file_name)
gdirs = utils.get_ref_mb_glaciers(gdirs)

# Save the list of glaciers for later
log.info('For RGIV{} and {} we have {} reference glaciers.'.format(rgi_version,
                                                                   baseline,
                                                                   len(gdirs)))
rgidf = pd.Series(data=[g.rgi_id for g in gdirs])
rgidf.to_csv(os.path.join(WORKING_DIR, 'mb_ref_glaciers.csv'))

# Prepro tasks
task_list = [
    tasks.glacier_masks,
    tasks.compute_centerlines,
    tasks.initialize_flowlines,
    tasks.catchment_area,
    tasks.catchment_intersections,
    tasks.catchment_width_geom,
    tasks.catchment_width_correction,
]
for task in task_list:
    execute_entity_task(task, gdirs)

# Climate tasks
tasks.compute_ref_t_stars(gdirs)
execute_entity_task(tasks.local_t_star, gdirs)
execute_entity_task(tasks.mu_star_calibration, gdirs)

# We store the associated params
mb_calib = gdirs[0].read_json('climate_info')['mb_calib_params']
with open(os.path.join(WORKING_DIR, 'mb_calib_params.json'), 'w') as fp:
    json.dump(mb_calib, fp)

# And also some statistics
utils.compile_glacier_statistics(gdirs)

# Tests: for all glaciers, the mass-balance around tstar and the
# bias with observation should be approx 0
for gd in gdirs:
    mb_mod = MultipleFlowlineMassBalance(gd,
                                         mb_model_class=ConstantMassBalance,
                                         use_inversion_flowlines=True,
                                         bias=0)  # bias=0 because of calib!
    try:
     mb = mb_mod.get_specific_mb()
    except ValueError as a:
     print (gd.rgi_id,file=sys.stderr)
     print (a,file=sys.stderr)
    try:
     print (gd.rgi_id,file=sys.stderr)
     np.testing.assert_allclose(mb, 0, atol=5)  # atol for numerical errors
    except AssertionError as a:
     print (a,file=sys.stderr)

    mb_mod = MultipleFlowlineMassBalance(gd, mb_model_class=PastMassBalance,
                                         use_inversion_flowlines=True)

    try:
     refmb = gd.get_ref_mb_data().copy()
     refmb.index = refmb.index[(refmb.index>=1979)&(refmb.index<=2019)]
     refmb['OGGM'] = mb_mod.get_specific_mb(year=refmb.index)
    except ValueError as a:
     print (gd.rgi_id,file=sys.stderr)
     print (a,file=sys.stderr)
    try:
     print (gd.rgi_id,file=sys.stderr)
     print ('bias: ' + str(refmb.OGGM.mean() - refmb.ANNUAL_BALANCE.mean()), file=sys.stderr)
     np.testing.assert_allclose(refmb.OGGM.mean(), refmb.ANNUAL_BALANCE.mean(),
                               atol=10)  # atol for numerical errors
    except AssertionError as a:
     print (a,file=sys.stderr)

# Log
log.info('Calibration is done!')
