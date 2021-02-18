# Python imports
import logging

# Libs
import geopandas as gpd
import pandas as pd
import shapely.geometry as shpg
import numpy as np
import matplotlib.pyplot as plt
from IPython import embed
import scipy.stats as stats
import os
from shutil import copyfile

# Locals
import oggm.cfg as cfg
from oggm import utils, workflow, tasks

import sys; sys.path.append("/exports/csce/datastore/geos/groups/geos_iceocean/dgoldber/oggm/rgi/glac_thickness");
from readTiff import readTiff


# For timing the run
import time
start = time.time()

# Module logger
log = logging.getLogger(__name__)

# Initialize OGGM and set up the default run parameters
cfg.initialize(logging_level='WORKFLOW')
rgi_version = '61'
#Change default swarm mu Value
cfg.PARAMS['swarm_mu']= 120.0
#Change the default temp_melt
cfg.PARAMS['temp_melt'] = -1.25
# This has thrown errors sometimes when I've reduced it so left at 80
cfg.PARAMS['border'] = 160
# Set this to make the run look for the ref_t* csv in the workflow directory instead of downloading it
cfg.PARAMS['run_mb_calibration'] = True
# Use multiprocessing?
cfg.PARAMS['use_multiprocessing'] = True
# Limit the number of cores, Shin and other geos servers have 64 so limit to half that
cfg.PARAMS['mp_processes'] = 24
# We are using which baseline data?
cfg.PARAMS['baseline_climate'] = 'CUSTOM'
# change the custom climate data
#cfg.PATHS['climate_file'] = '/exports/csce/datastore/geos/groups/geos_iceocean/dgoldber/oggm/SWARM_input/oggm_SWARM_input_cru_ref_hgts/oggm_cru_hgt_input.nc'
cfg.PATHS['climate_file'] = '/exports/csce/datastore/geos/groups/geos_iceocean/dgoldber/oggm/SWARM_inputs/swarm_cru_hgt_input_padmonths.nc'
# Set to True for operational runs - here we want all glaciers to run
cfg.PARAMS['continue_on_error'] = True
# Change the minimum timestep
cfg.PARAMS['cfl_min_dt'] = 10.0
# Local working directory (where OGGM will write its output) Need to create this beforehand and put the mass balance data in (ref_t_stars.csv).
WORKING_DIR = '/exports/csce/datastore/geos/groups/geos_iceocean/dgoldber/oggm/oggm_runs/thickness_test'
if not os.path.exists(WORKING_DIR):
    os.makedirs(WORKING_DIR)
REF_DIR = '/exports/csce/datastore/geos/groups/geos_iceocean/dgoldber/oggm/oggm_runs/oggm_mb_calibration_ERA5_MSWEP_reference_period'
copyfile(REF_DIR + '/ref_tstars.csv',WORKING_DIR + '/ref_tstars.csv')
cfg.PATHS['working_dir'] = WORKING_DIR


# RGI file setup, easiest but very inelegant way to do this atm is to make a list of all 13,14,15 RGI glaciers then filtering
#Start with 13,14 and 15
path = utils.get_rgi_region_file('13', version=rgi_version)
rgidf_13 = gpd.read_file(path)
path = utils.get_rgi_region_file('14', version=rgi_version)
rgidf_14 = gpd.read_file(path)
path = utils.get_rgi_region_file('15', version=rgi_version)
rgidf_15 = gpd.read_file(path)
#Now combine the files
# rgidf_temp = rgidf_13.concat(rgidf_14)
# rgidf = rgidf_temp.concat(rgidf_15)
rgidf = gpd.GeoDataFrame(pd.concat([rgidf_13, rgidf_14, rgidf_15]))


#rgidf = rgidf[rgidf['RGIId']=='RGI60-13.00004']
# Get the shapefile
basin = gpd.read_file('/exports/csce/datastore/geos/groups/geos_iceocean/dgoldber/oggm/SWARM-OGGM-Model/shape_files/swarm_grid.shp')
print('got glacier shp')

# Take all glaciers in the desired areas
in_bas = [basin.geometry.contains(shpg.Point(x, y))[0] for
          (x, y) in zip(rgidf.CenLon, rgidf.CenLat)]
rgidf = rgidf.loc[in_bas]



# Sort for more efficient parallel computing
rgidf = rgidf.sort_values('Area', ascending=False)
# Get rid of smaller glaciers, area is in km^2

#rgidf = rgidf[rgidf.Area >= 0.15]

embed()
#rgidf = rgidf.iloc[np.random.choice ( len(rgidf), size=40, replace=False )]

#This is old code that doesn't work and is less efficient but left here just in case
#rgidf = rgidf.drop(rgidf[rgidf.Area < 1.0].index)

##Added in section to convert to pandas dataframe and print out file if needed/check data
#rgidf.to_file('/Users/louis/workflow_test/rgidf.shp')

"""

##Convert to pandas dataframe then print
#rgipdf = pd.DataFrame(rgidf)
#rgipdf.to_csv('/Users/louis/workflow_test/rgipdf.csv')

#Start the run
log.workflow('Starting OGGM run')
log.workflow('Number of glaciers: {}'.format(len(rgidf)))

#Go get the pre-processed glacier directories, this will eventually be reduced to prepro = 1 but for now need the climate files in the glacier directories.
gdirs = workflow.init_glacier_directories(rgidf, from_prepro_level=1)
#Run the prepro tasks on the dem to get flowlines using the generic processing etc. This could be broken down if we don't want certain tasks.
workflow.gis_prepro_tasks(gdirs)

#Process the climate data
log.info('Process the climate data...')
workflow.execute_entity_task(tasks.process_climate_data, gdirs)
#Run the climate calibrations based on the new mass balance data, this mass balance data should be in the Working directory folder.
workflow.execute_entity_task(tasks.local_t_star, gdirs)
workflow.execute_entity_task(tasks.mu_star_calibration, gdirs)
#Run the inversion tools ahead of the model runs
#First set the paramters we can change, currently OGGM default ones:
# Deformation: from Cuffey and Patterson 2010
glen_a = 2.4e-24
# Sliding: from Oerlemans 1997
fs = 5.7e-20
#Prep the data for inversion
workflow.execute_entity_task(tasks.prepare_for_inversion, gdirs)
#Run the inversion

factors = [3]
slide_med_err = np.empty(len(factors))
noslide_med_err = np.empty(len(factors))
slide_iqr_err = np.empty(len(factors))
noslide_iqr_err = np.empty(len(factors))
for i_f in range(len(factors)):

 suf = '_{:03d}_with_fs'.format(int(factors[i_f] * 10))
 workflow.execute_entity_task(tasks.mass_conservation_inversion, gdirs,
                                 glen_a=glen_a*factors[i_f], fs=fs)
 workflow.execute_entity_task(tasks.filter_inversion_output, gdirs)
 utils.compile_glacier_statistics(gdirs, filesuffix=suf,inversion_only=True)

 df = pd.read_csv(WORKING_DIR + '/glacier_statistics' + suf + '.csv')

 err_vol = np.empty(len(df))

 for i in range(len(df)):
  rgiid = df['rgi_id'][i]
  vol = df['inv_volume_km3'][i]
  reg = rgiid[7]
  folder = '/exports/csce/datastore/geos/groups/geos_iceocean/dgoldber/oggm/rgi/glac_thickness/RGI60-1' + reg + '/'
  filename = folder + rgiid + '_thickness.tif'
  data,xOrigin,yOrigin,pixelWidth,pixelHeight,nX,nY = readTiff(filename)
  volHF = np.sum(data) * np.abs(pixelWidth*pixelHeight) / 1.e9
  err_vol[i] = (volHF - vol)/volHF

 I = ~np.isnan(err_vol)
 mean_err = np.median(err_vol[I])
 std_err = stats.iqr(err_vol[I]) 
 plt.hist(err_vol[I],bins=30);
 plt.savefig(WORKING_DIR + '/volume_err' + suf + '.png') 
 slide_med_err[i_f] = mean_err
 slide_iqr_err[i_f] = std_err
 plt.close('all')

 suf = '_{:03d}_without_fs'.format(int(factors[i_f] * 10))
 workflow.execute_entity_task(tasks.mass_conservation_inversion, gdirs,
                                 glen_a=glen_a*factors[i_f], fs=0)
 workflow.execute_entity_task(tasks.filter_inversion_output, gdirs)
 utils.compile_glacier_statistics(gdirs, filesuffix=suf,inversion_only=True)

 df = pd.read_csv(WORKING_DIR + '/glacier_statistics' + suf + '.csv')

 err_vol = np.empty(len(df))

 for i in range(len(df)):
  rgiid = df['rgi_id'][i]
  vol = df['inv_volume_km3'][i]
  reg = rgiid[7]
  folder = '/exports/csce/datastore/geos/groups/geos_iceocean/dgoldber/oggm/rgi/glac_thickness/RGI60-1' + reg + '/'
  filename = folder + rgiid + '_thickness.tif'
  data,xOrigin,yOrigin,pixelWidth,pixelHeight,nX,nY = readTiff(filename)
  volHF = np.sum(data) * np.abs(pixelWidth*pixelHeight) / 1.e9
  err_vol[i] = (volHF - vol)/volHF
 
 I = ~np.isnan(err_vol)
 mean_err = np.median(err_vol[I])
 std_err = stats.iqr(err_vol[I]) 
 plt.hist(err_vol[I],bins=30);
 plt.savefig(WORKING_DIR + '/volume_err' + suf + '.png') 
 noslide_med_err[i_f] = mean_err
 noslide_iqr_err[i_f] = std_err
 plt.close('all')

#for i in range(len(factors)):
# print('Sliding, factor = ' + str(factors[i]) + ': median = ' + str(slide_med_err[i]) + '; IQR = ' + str(slide_iqr_err[i]))
# print('No Sliding, factor = ' + str(factors[i]) + ': median = ' + str(noslide_med_err[i]) + '; IQR = ' + str(noslide_iqr_err[i]))



#Compile all the previous tasks to give the output ready as input for a model run
workflow.execute_entity_task(tasks.init_present_time_glacier, gdirs)


#Run the model run, this can be changed to 'run_from_climate_data' for our runs
workflow.execute_entity_task(tasks.run_from_climate_data, gdirs,
                             store_monthly_step=True,
                             output_filesuffix='_commitment',ys=1979,ye=2020)

# workflow.execute_entity_task(tasks.compile_run_output, gdirs)
# workflow.execute_entity_task(tasks.compile_climate_input, gdirs)




# Log
m, s = divmod(time.time() - start, 60)
h, m = divmod(m, 60)
log.workflow('OGGM is done! Time needed: %d:%02d:%02d' % (h, m, s))
"""
