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
##Added in section to convert to pandas dataframe and print out file if needed/check data
#rgidf.to_file('/Users/louis/workflow_test/rgidf.shp')
##Convert to pandas dataframe then print
#rgipdf = pd.DataFrame(rgidf)
#rgipdf.to_csv('/Users/louis/workflow_test/rgipdf.csv')

#Start the run
log.workflow('Starting OGGM run')
log.workflow('Number of glaciers: {}'.format(len(rgidf)))

#Go get the pre-processed glacier directories, this will eventually be reduced to prepro = 1 but for now need the cliamte files in the glacier directories.
gdirs = workflow.init_glacier_directories(rgidf, from_prepro_level=2)
#Run the prepro tasks on the dem to get flowlines using the generic processing etc.
workflow.gis_prepro_tasks(gdirs)
#Run the climate calibrations based on the new mass balance data
workflow.execute_entity_task(tasks.local_t_star, gdirs)
workflow.execute_entity_task(tasks.mu_star_calibration, gdirs)



#Run the inversion tools ahead of the model runs
#First set the paramters we can change:
# Deformation: from Cuffey and Patterson 2010
glen_a = 2.4e-24
# Sliding: from Oerlemans 1997
fs = 5.7e-20
#Prep the data for inversion
workflow.execute_entity_task(tasks.prepare_for_inversion, gdirs)
#Run the inversion
workflow.execute_entity_task(tasks.mass_conservation_inversion, gdirs,
                                 glen_a=glen_a, fs=fs)
#Apparently can be an issue wioth the last few grid points being noisy or having a negative slope so can filter these without affecting the total volume
workflow.execute_entity_task(tasks.filter_inversion_output, gdirs)



#Compile all the previous tasks to give the output ready for a model run
workflow.execute_entity_task(tasks.init_present_time_glacier, gdirs)


#Run the model run
workflow.execute_entity_task(tasks.run_random_climate, gdirs,
                             nyears=300, y0=2000, seed=1,store_monthly_step=True,
                             output_filesuffix='_commitment')

# Log
m, s = divmod(time.time() - start, 60)
h, m = divmod(m, 60)
log.workflow('OGGM is done! Time needed: %d:%02d:%02d' % (h, m, s))
