import rasterio as rs
from matplotlib import pyplot as plt

working_dir = '/exports/csce/datastore/geos/groups/geos_iceocean/kinnear/SWARM-OGGM-Model/figure_scripts/'

data_path='/exports/csce/datastore/geos/groups/geos_iceocean/kinnear/oggm_run_data_for_swarm/script_testing/test.nc'

data = rs.open('netcdf:'+data_path+':runoff_area_normalised_data')

raster = data.read(2)

plt.imshow(raster, cmap='jet')
plt.savefig(working_dir+'test_raster.png')
plt.close()