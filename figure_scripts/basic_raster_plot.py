import rasterio as rs
import xarray as xr
import matplotlib as mpl
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

working_dir = '/exports/csce/datastore/geos/groups/geos_iceocean/kinnear/SWARM-OGGM-Model/figure_scripts/'

data_path='/exports/csce/datastore/geos/groups/geos_iceocean/kinnear/oggm_run_data_for_swarm/script_testing/test.nc'

dataset = xr.open_dataset(data_path,decode_times=False)
lat = dataset.lat.values
lon = dataset.lon.values

data = rs.open('netcdf:'+data_path+':runoff_area_normalised_data')


raster = data.read(2)
fig = plt.figure()
ax=fig.add_subplot(1,1,1)

cb = ax.imshow(raster, cmap='RdBu',extent=(min(lon),max(lon),min(lat),max(lat)))

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cb, cax=cax, label='Runoff')


ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
plt.savefig(working_dir+'test_raster.png')
plt.close()