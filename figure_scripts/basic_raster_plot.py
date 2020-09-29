import rasterio as rs
import rioxarray as rio
import xarray as xr
import matplotlib as mpl
import geopandas as gpd
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

working_dir = '/exports/csce/datastore/geos/groups/geos_iceocean/kinnear/SWARM-OGGM-Model/figure_scripts/'
data_path='/exports/csce/datastore/geos/groups/geos_iceocean/kinnear/oggm_run_data_for_swarm/oggm_mswep_era_reference_run_90/oggm_mswep_era_reference_run_90_data.nc'

shape_path = '/exports/csce/datastore/geos/groups/geos_iceocean/kinnear/SWARM-OGGM-Model/shape_files/RiverBasinsMerged.shp'
china = gpd.read_file(shape_path)

dataset = xr.open_dataset(data_path,decode_times=False)
lat = dataset.lat.values
lon = dataset.lon.values
total = dataset.runoff_area_normalised_data.mean(dim='time')
total = total.rio.to_raster('total.tif')
#data = rs.open('netcdf:'+data_path+':runoff_area_normalised_data')

data = rs.open('total.tif')
raster = data.read(1)
fig = plt.figure()
ax=fig.add_subplot(1,1,1)

china.plot(ax=ax,color='none', edgecolor='black',linewidth=0.25)

cb = ax.imshow(raster, cmap='RdBu',extent=(min(lon),max(lon),min(lat),max(lat)), vmax=0.5)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cb, cax=cax, label='Runoff')


ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
plt.savefig(working_dir+'test_raster.png')
plt.close()