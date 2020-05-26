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

#Initialize OGGM and set up the default run parameters
cfg.initialize(logging_level='WORKFLOW')
rgi_version = '61'
rgi_region = '15'  

# No need for a big map here
cfg.PARAMS['border'] = 80

# Local working directory (where OGGM will write its output)
WORKING_DIR = '/Users/louis/Prepro_test/'
utils.mkdir(WORKING_DIR, reset=True)
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
# Get rid of smaller glaciers
rgidf = rgidf.drop(rgidf[rgidf.Area < 0.5].index)
#Added in section to convert to pandas dataframe and print out file sol only need to run once
rgidf.to_file('/Users/louis/Prepro_test/rgidf.shp')
#COnvert to pandas dataframe then print
rgipdf = pd.DataFrame(rgidf)
rgipdf.to_csv('/Users/louis/Prepro_test/rgipdf.csv')

#Loop thorugh and do for each prepro directory
WORKING_DIR = '/Users/louis/Prepro_test/prepro_1'
utils.mkdir(WORKING_DIR, reset=True)
cfg.PATHS['working_dir'] = WORKING_DIR
gdirs = workflow.init_glacier_directories(rgidf, from_prepro_level=1)

WORKING_DIR = '/Users/louis/Prepro_test/prepro_2'
utils.mkdir(WORKING_DIR, reset=True)
cfg.PATHS['working_dir'] = WORKING_DIR
gdirs = workflow.init_glacier_directories(rgidf, from_prepro_level=2)

WORKING_DIR = '/Users/louis/Prepro_test/prepro_3'
utils.mkdir(WORKING_DIR, reset=True)
cfg.PATHS['working_dir'] = WORKING_DIR
gdirs = workflow.init_glacier_directories(rgidf, from_prepro_level=3)

WORKING_DIR = '/Users/louis/Prepro_test/prepro_4'
utils.mkdir(WORKING_DIR, reset=True)
cfg.PATHS['working_dir'] = WORKING_DIR
gdirs = workflow.init_glacier_directories(rgidf, from_prepro_level=4)

