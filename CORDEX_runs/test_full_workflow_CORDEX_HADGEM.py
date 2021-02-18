# Python imports
import logging

# Libs
import os
import geopandas as gpd
import pandas as pd
import shapely.geometry as shpg
from shutil import copyfile

# Locals
import oggm.cfg as cfg
from oggm import utils, workflow, tasks
from oggm.core import gcm_climate

from IPython import embed

# For timing the run
import time
start = time.time()

# Module logger
log = logging.getLogger(__name__)

# Initialize OGGM and set up the default run parameters
cfg.initialize(logging_level='WORKFLOW')
rgi_version = '61'

#Change default swarm mu Value
#cfg.PARAMS['swarm_mu']= 90.0000
#Change the default temp_melt
cfg.PARAMS['temp_melt'] = -1.25
# Set to max to avoid eroors from glaciers exceeding boundaries
cfg.PARAMS['border'] = 160
# Set this to make the run look for the ref_t* csv in the workflow directory instead of downloading it
#cfg.PARAMS['run_mb_calibration'] = True
# Use multiprocessing?
cfg.PARAMS['use_multiprocessing'] = True
# Limit the number of cores, Shin and other geos servers have 64 so limit to half that
cfg.PARAMS['mp_processes'] = 40
# We are using which baseline data?
cfg.PARAMS['baseline_climate'] = 'CUSTOM'
# change the custom climate data
#cfg.PATHS['climate_file'] = '/exports/csce/datastore/geos/groups/geos_iceocean/kinnear/SWARM_files/oggm_SWARM_input_cru_ref_hgts/oggm_cru_hgt_input.nc'
# Set to True for operational runs - here we want all glaciers to run
cfg.PARAMS['continue_on_error'] = True
# Change the minimum timestep
cfg.PARAMS['cfl_min_dt'] = 10

cfg.BASENAMES.update({'gcm_data':'gcm_data_hadgem.nc'})
# Local working directory (where OGGM will write its output) Need to create this beforehand and put the mass balance data in (ref_t_stars.csv).
WORKING_DIR = '/exports/csce/datastore/geos/groups/geos_iceocean/dgoldber/oggm/oggm_runs/thickness_test'
if not os.path.exists(WORKING_DIR):
    os.makedirs(WORKING_DIR)
REF_DIR = '/exports/csce/datastore/geos/groups/geos_iceocean/dgoldber/oggm/oggm_runs/oggm_mb_calibration_ERA5_MSWEP_reference_period'
copyfile(REF_DIR + '/ref_tstars.csv',WORKING_DIR + '/ref_tstars.csv')
cfg.PATHS['working_dir'] = WORKING_DIR



# RGI file setup, easiest but very inelegant way to do this atm is to make a list of all 13,14,15 RGI glaciers then filtering
#Start with 13,14 and 15
#path = utils.get_rgi_region_file('13', version=rgi_version)
#rgidf_13 = gpd.read_file(path)
#path = utils.get_rgi_region_file('14', version=rgi_version)
#rgidf_14 = gpd.read_file(path)
#path = utils.get_rgi_region_file('15', version=rgi_version)
#rgidf_15 = gpd.read_file(path)
#Now combine the files
# rgidf_temp = rgidf_13.concat(rgidf_14)
# rgidf = rgidf_temp.concat(rgidf_15)
#rgidf = gpd.GeoDataFrame(pd.concat([rgidf_13, rgidf_14, rgidf_15]))
#print(rgidf)
# Get the shapefile
#basin = gpd.read_file('/exports/csce/datastore/geos/groups/geos_iceocean/kinnear/SWARM-OGGM-Model/shape_files/swarm_grid.shp')
#print('got glacier shp')

# Take all glaciers in the desired areas
#in_bas = [basin.geometry.contains(shpg.Point(x, y))[0] for
#          (x, y) in zip(rgidf.CenLon, rgidf.CenLat)]
#rgidf = rgidf.loc[in_bas]



# Sort for more efficient parallel computing
#rgidf = rgidf.sort_values('Area', ascending=False)
# Get rid of smaller glaciers, area is in km^2

#rgidf = rgidf[rgidf.Area >= 1.0]

#This is old code that doesn't work and is less efficient but left here just in case
#rgidf = rgidf.drop(rgidf[rgidf.Area < 1.0].index)

##Added in section to convert to pandas dataframe and print out file if needed/check data
#rgidf.to_file('/Users/louis/workflow_test/rgidf.shp')


##Convert to pandas dataframe then print
#rgipdf = pd.DataFrame(rgidf)
#rgipdf.to_csv('/Users/louis/workflow_test/rgipdf.csv')

#Start the run
log.workflow('Starting OGGM run')
#log.workflow('Number of glaciers: {}'.format(len(rgidf)))

#Go get the pre-processed glacier directories, this will eventually be reduced to prepro = 1 but for now need the climate files in the glacier directories.
gdirs = workflow.init_glacier_directories()
#Run the prepro tasks on the dem to get flowlines using the generic processing etc. This could be broken down if we don't want certain tasks.
#workflow.gis_prepro_tasks(gdirs)
#Process the climate data
#log.info('Process the climate data...')
#workflow.execute_entity_task(tasks.process_custom_climate_data, gdirs,y1=2018)
#Run the climate calibrations based on the new mass balance data, this mass balance data should be in the Working directory folder.
#workflow.execute_entity_task(tasks.local_t_star, gdirs)
#workflow.execute_entity_task(tasks.mu_star_calibration, gdirs)
#Run the inversion tools ahead of the model runs
#First set the paramters we can change, currently OGGM default ones:
# Deformation: from Cuffey and Patterson 2010
#glen_a = 2.4e-24
## Sliding: from Oerlemans 1997
#fs = 5.7e-20
#Prep the data for inversion
#workflow.execute_entity_task(tasks.prepare_for_inversion, gdirs)
#Run the inversion
#workflow.execute_entity_task(tasks.mass_conservation_inversion, gdirs,
#                                 glen_a=glen_a, fs=fs)
#Apparently can be an issue with the last few grid points being noisy or having a negative slope so can filter these without affecting the total volume
#workflow.execute_entity_task(tasks.filter_inversion_output, gdirs)

#Compile all the previous tasks to give the output ready as input for a model run
workflow.execute_entity_task(tasks.init_present_time_glacier, gdirs)

workflow.execute_entity_task(tasks.run_from_climate_data, gdirs,
                             store_monthly_step=True,ye=2018,
                             output_filesuffix='_spinup')

#The pathways to the CORDEX input Data
#temp_path = '/exports/csce/datastore/geos/groups/geos_iceocean/kinnear/CORDEX/HadGEM2/QM_tas_CORDEX_EA_HadGEM2_swarm_domain_monthly_360-day_calendar_19700101-20991230.nc'
#precip_path = '/exports/csce/datastore/geos/groups/geos_iceocean/kinnear/CORDEX/HadGEM2/QM_pr_CORDEX_EA_HadGEM2_swarm_domain_monthly_360-day_calendar_19700101-20991230.nc'

temp_path = '/exports/csce/datastore/geos/groups/geos_iceocean/kinnear/CORDEX/MPI-M/QM_tas_CORDEX_EA_MPI-M-MPI-ESM-LR_swarm_domain_calendar_monthly_19700101-21001231_k.nc' 
precip_path = '/exports/csce/datastore/geos/groups/geos_iceocean/kinnear/CORDEX/MPI-M/QM_pr_CORDEX_EA_MPI-M-MPI-ESM-LR_swarm_domain_calendar_monthly_19700101-21001231.nc'


workflow.execute_entity_task(gcm_climate.process_swarm_data, gdirs, fpath_temp=temp_path,fpath_precip=precip_path, gcm_anomaly_range=('1971','2100'))
#Run the model run, this can be changed to 'run_from_climate_data' for our runs
workflow.execute_entity_task(tasks.run_from_climate_data, gdirs,store_monthly_step=True,ys=2018, ye=2100, climate_filename='gcm_data',
                             output_filesuffix='_hadgem',init_model_filesuffix='_spinup')

#workflow.execute_entity_task(tasks.compile_run_output, gdirs)
#workflow.execute_entity_task(tasks.compile_climate_input, gdirs)

# Log
m, s = divmod(time.time() - start, 60)
h, m = divmod(m, 60)
log.workflow('OGGM is done! Time needed: %d:%02d:%02d' % (h, m, s))
