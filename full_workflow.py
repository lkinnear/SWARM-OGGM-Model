# Python imports
import logging

# Libs
import geopandas as gpd
import pandas as pd
import shapely.geometry as shpg

# Locals
import oggm.cfg as cfg
from oggm import utils, workflow, tasks


# For timing the run
import time
start = time.time()

# Module logger
log = logging.getLogger(__name__)

# Initialize OGGM and set up the default run parameters
cfg.initialize(logging_level='WORKFLOW')
rgi_version = '61'
rgi_region = '15' 

# This has thrown errors sometimes when I've reduced it so left at 80
cfg.PARAMS['border'] = 80
# Set this to make the run look for the ref_t* csv in the workflow directory instead of downloading it
cfg.PARAMS['run_mb_calibration'] = True

# Local working directory (where OGGM will write its output) Need to create this beforehand in order to put the mass balance data in.
WORKING_DIR = '/Users/louis/workflow_test/'
cfg.PATHS['working_dir'] = WORKING_DIR

# RGI file
path = utils.get_rgi_region_file(rgi_region, version=rgi_version)
rgidf = gpd.read_file(path)

# Get the shapefile
basin = gpd.read_file('/Users/louis/china_test.shp')
print('got glacier shp')

# Take all glaciers in the desired areas
in_bas = [basin.geometry.contains(shpg.Point(x, y))[0] for
          (x, y) in zip(rgidf.CenLon, rgidf.CenLat)]
rgidf = rgidf.loc[in_bas]
print('found basins')
# Sort for more efficient parallel computing
rgidf = rgidf.sort_values('Area', ascending=False)
# Get rid of smaller glaciers, area is in km^2
rgidf = rgidf.drop(rgidf[rgidf.Area < 1.0].index)
#Added in section to convert to pandas dataframe and print out file sol only need to run once
rgidf.to_file('/Users/louis/workflow_test/rgidf.shp')
#COnvert to pandas dataframe then print
rgipdf = pd.DataFrame(rgidf)
rgipdf.to_csv('/Users/louis/workflow_test/rgipdf.csv')

log.workflow('Starting OGGM run')
log.workflow('Number of glaciers: {}'.format(len(rgidf)))

# Go - get the pre-processed glacier directories
gdirs = workflow.init_glacier_directories(rgidf, from_prepro_level=2)

workflow.gis_prepro_tasks(gdirs)
workflow.execute_entity_task(tasks.local_t_star, gdirs)
workflow.execute_entity_task(tasks.mu_star_calibration, gdirs)
workflow.inversion_tasks(gdirs)
workflow.execute_entity_task(tasks.init_present_time_glacier, gdirs)
workflow.execute_entity_task(tasks.run_random_climate, gdirs,
                             nyears=300, y0=2000, seed=1,store_monthly_step=True,
                             output_filesuffix='_commitment')

# Log
m, s = divmod(time.time() - start, 60)
h, m = divmod(m, 60)
log.workflow('OGGM is done! Time needed: %d:%02d:%02d' % (h, m, s))
