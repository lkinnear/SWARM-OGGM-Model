# Python imports
import logging

# Libs
import geopandas as gpd
import pandas as pd
import shapely.geometry as shpg
import sys
from IPython import embed

# Locals
import oggm.cfg as cfg
from oggm import utils, workflow, tasks
from oggm.core import gcm_climate

# For timing the run
import time
start = time.time()

inlist = sys.argv[1:]
if (len(inlist)>0):
  while (len(inlist)>0):
    print(inlist[0])
    flag = inlist[0]
    if (flag=='-precfile'):
     if len(inlist)==1: 
      print ('provide prec file')
      exit()
     else:
      precip_path = inlist[1]
      if len(inlist)==2:
       break
      inlist = inlist[2:]
    elif (flag=='-tempfile'):
     if len(inlist)==1: 
      print ('provide temp file')
      exit()
     else:
      temp_path = inlist[1]
      if len(inlist)==2:
       break
      inlist = inlist[2:]
    elif (flag=='-suf'):
     if len(inlist)==1: 
      print ('provide output suffix')
      exit()
     else:
      gcm_output_suffix = inlist[1]
      if len(inlist)==2:
       break
      inlist = inlist[2:]
    elif (flag=='-mu'):
     if len(inlist)==1: 
      print ('provide mu')
      exit()
     else:
      swarm_mu = inlist[1]
      if len(inlist)==2:
       break
      inlist = inlist[2:]
    else: 
      print ('unrecognised flag')
      exit()
      
print(temp_path)
print(precip_path)
print(swarm_mu)
print(gcm_output_suffix)

#temp_path = '/exports/csce/datastore/geos/groups/geos_iceocean/kinnear/CORDEX/HadGEM2/QM_tas_CORDEX_EA_HadGEM2_swarm_domain_monthly_360-day_calendar_19700101-20991230.nc'
#precip_path = '/exports/csce/datastore/geos/groups/geos_iceocean/kinnear/CORDEX/HadGEM2/QM_pr_CORDEX_EA_HadGEM2_swarm_domain_monthly_360-day_calendar_19700101-20991230.nc'
#gcm_output_suffix = '_hadgem.nc'

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
cfg.PARAMS['border'] = 250
# Set this to make the run look for the ref_t* csv in the workflow directory instead of downloading it
#cfg.PARAMS['run_mb_calibration'] = True
# Use multiprocessing?
cfg.PARAMS['use_multiprocessing'] = True
# Limit the number of cores, Shin and other geos servers have 64 so limit to half that
cfg.PARAMS['mp_processes'] = 24
# We are using which baseline data?
cfg.PARAMS['baseline_climate'] = 'CUSTOM'
# change the custom climate data
#cfg.PATHS['climate_file'] = '/exports/csce/datastore/geos/groups/geos_iceocean/kinnear/SWARM_files/oggm_SWARM_input_cru_ref_hgts/oggm_cru_hgt_input.nc'
# Set to True for operational runs - here we want all glaciers to run
cfg.PARAMS['continue_on_error'] = True
# Change the minimum timestep
cfg.PARAMS['cfl_min_dt'] = 10

# Local working directory (where OGGM will write its output) Need to create this beforehand and put the mass balance data in (ref_t_stars.csv).
cfg.BASENAMES.update({'gcm_data':'gcm_data' + gcm_output_suffix + '.nc'})
# Local working directory (where OGGM will write its output) Need to create this beforehand and put the mass balance data in (ref_t_stars.csv).
WORKING_DIR = '/exports/csce/datastore/geos/groups/geos_iceocean/dgoldber/oggm/oggm_runs/oggm_mswep_era_reference_run_' + str(swarm_mu) + '_geodetic'
#WORKING_DIR = '/exports/csce/datastore/geos/groups/geos_iceocean/dgoldber/oggm/oggm_runs/thickness_test'
cfg.PATHS['working_dir'] = WORKING_DIR


#Start the run
log.workflow('Starting OGGM run')

path = utils.get_rgi_region_file('13', version=rgi_version)
rgidf_13 = gpd.read_file(path)
path = utils.get_rgi_region_file('14', version=rgi_version)
rgidf_14 = gpd.read_file(path)
path = utils.get_rgi_region_file('15', version=rgi_version)
rgidf_15 = gpd.read_file(path)
rgidf = gpd.GeoDataFrame(pd.concat([rgidf_13, rgidf_14, rgidf_15]))
basin = gpd.read_file('/exports/csce/datastore/geos/groups/geos_iceocean/kinnear/SWARM-OGGM-Model/shape_files/swarm_grid.shp')
in_bas = [basin.geometry.contains(shpg.Point(x, y))[0] for
          (x, y) in zip(rgidf.CenLon, rgidf.CenLat)]
rgidf = rgidf.loc[in_bas]
rgidf = rgidf.sort_values('Area', ascending=False)
rgidf = rgidf[rgidf.Area >= 0.15]
gdirs = workflow.init_glacier_directories()
#Compile all the previous tasks to give the output ready as input for a model run
workflow.execute_entity_task(tasks.init_present_time_glacier, gdirs)

workflow.execute_entity_task(tasks.run_from_climate_data, gdirs,
                             store_monthly_step=True,ye=2018,
                             output_filesuffix='_spinup',climate_filename='climate_historical')


workflow.execute_entity_task(gcm_climate.process_swarm_data, gdirs, fpath_temp=temp_path,fpath_precip=precip_path, gcm_anomaly_range=('1971','2098'))
#Run the model run, this can be changed to 'run_from_climate_data' for our runs
workflow.execute_entity_task(tasks.run_from_climate_data, gdirs,store_monthly_step=True,ys=2018, ye=2098, climate_filename='gcm_data',
                             mb_elev_feedback='monthly',output_filesuffix='_forced'+gcm_output_suffix,init_model_filesuffix='_spinup')

#workflow.execute_entity_task(tasks.compile_run_output, gdirs)
#workflow.execute_entity_task(tasks.compile_climate_input, gdirs)

# Log
m, s = divmod(time.time() - start, 60)
h, m = divmod(m, 60)
log.workflow('OGGM is done! Time needed: %d:%02d:%02d' % (h, m, s))
