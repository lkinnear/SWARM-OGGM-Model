import os as os
import numpy as np
import pandas as pd
import xarray as xr

#Set up the working directory
run_name='cru_rf_hgt_custom_climate_runs'
#os.mkdir('/exports/csce/datastore/geos/groups/geos_iceocean/kinnear/oggm_run_data_for_swarm/'+run_name)
working_dir = '/exports/csce/datastore/geos/groups/geos_iceocean/kinnear/oggm_run_data_for_swarm/'+run_name
#Now locate the raw dataset
run_name = 'oggm_custom_climate_cru_rf_hgt'
raw_data_directory = '/exports/csce/datastore/geos/groups/geos_iceocean/kinnear/oggm_runs/'
raw_data_folder = raw_data_directory+run_name+'/per_glacier/'
#Now give a it an output to make sure it's running properly and you can check
print('Processing the data from '+raw_data_folder+' to output in '+working_dir)
#Set up files to opened

rgi_id = []
num=1
#Get list of RGI-IDs t0 use
for root,dirs,files in os.walk(raw_data_folder, topdown=False):
    for name in files:
        if name == 'model_diagnostics_commitment.nc':
            path = os.path.join(root,name)
            rgi = os.path.basename(os.path.dirname(path))
            rgi = str(rgi)
            rgi_id.append(rgi)

#Hacky way to start the process by making a file to input the data scraped from the loop into.
for root, dirs, files in os.walk(raw_data_folder, topdown=False):
    for name in files:
       if num == 1:
           if name == 'model_diagnostics_commitment.nc':
              fname = os.path.join(root,name)
              with xr.open_dataset(fname) as model_run:
                  time = model_run.time.values
                  yrs = model_run.hydro_year.values
                  months = model_run.hydro_month.values
                  cyrs = model_run.calendar_year.values
                  cmonths = model_run.calendar_month.values

                  #Create the file for writing
                  ds = xr.Dataset()

                  print('opened number: {}'.format(num))
                  #Global attributes
                  ds.attrs['description'] = 'OGGM model output'
                  ds.attrs['calendar'] = '365-day no leap'


                  # Coordinates

                  ds.coords['time'] = ('time', time)
                  ds.coords['time'].astype(float)
                  ds.coords['rgi_id'] = ('rgi_id', rgi_id)
                  ds.coords['hydro_year'] = ('time', yrs)
                  ds.coords['hydro_month'] = ('time', months)
                  ds.coords['calendar_year'] = ('time', cyrs)
                  ds.coords['calendar_month'] = ('time', cmonths)
                  ds['time'].attrs['description'] = 'Floating hydrological year'
                  ds['rgi_id'].attrs['description'] = 'RGI glacier identifier'
                  ds['hydro_year'].attrs['description'] = 'Hydrological year'
                  ds['hydro_month'].attrs['description'] = 'Hydrological month'
                  ds['calendar_year'].attrs['description'] = 'Calendar year'
                  ds['calendar_month'].attrs['description'] = 'Calendar month'

                  #Create the dimensions
                  shape = (len(time), len(rgi_id))
                  print('There are {} glaciers that contain run data'.format(len(rgi_id)))
                  # These variables are always available
                  vol = np.zeros(shape)
                  area = np.zeros(shape)
                  precip = np.zeroes(shape)
                  latitude =np.zeros(len(rgi_id))
                  longitude =np.zeros(len(rgi_id))



                  for i in range (len(rgi_id)):
                                 #Do some trickery to get the paths and speed it up over loops
                                  fold = rgi_id[i][:8]
                                  subfold = rgi_id[i][:11]
                                  folpath = raw_data_folder+'/'+fold+'/'+subfold+'/'+rgi_id[i]+'/model_diagnostics_commitment.nc'
                                  #print('found glacier number: {}'.format(i))
                                  with xr.open_dataset(folpath) as ds_diag:
                                      vol[:, i] = ds_diag.volume_m3.values
                                      area[:, i] = ds_diag.area_m2.values
                                  #open the file to get lat and lon
                                  folpath = raw_data_folder+'/'+fold+'/'+subfold+'/'+rgi_id[i]+'/climate_historical.nc'
                                  with xr.open_dataset(folpath) as ds_diag:
                                      lat = ds_diag.ref_pix_lat
                                      lon = ds_diag.ref_pix_lon

                                      latitude[i] = lat
                                      longitude[i] = lon
                                      latitude.astype(float)
                                      longitude.astype(float)
                                      precip[:, i] = ds_diag.prcp.values

                                 #Now get the rainfall




                  #idx = pd.MultiIndex.from_arrays(arrays=[latitude,longitude], names=["lat","lon"])
                  #volume = pd.Series(data=vol, index=idx)
                  #ds['volume'] = xr.DataArray.from_series(volume)
                  #ds.coords['lat'] = ('lat', latitude)
                  #ds.coords['lon'] = ('lon', longitude)
                  #lat = xr.DataArray(latitude, [('latitude', latitude)])
                  #lon = xr.DataArray(longitude, [('longitude', longitude)])
                  ds['volume'] = ((('time', 'rgi_id'), vol))
                  #ds['volume'] = ds['volume'].expand_dims(lon=lon)
                  #ds['volume'] = ds['volume'].expand_dims(lat=lat)
                  ds['volume'].attrs['description'] = 'Total glacier volume'
                  ds['volume'].attrs['units'] = 'm 3'


                  ds['area'] = (('time','rgi_id'), area)
                  ds['area'].attrs['description'] = 'Total glacier area'
                  ds['area'].attrs['units'] = 'm 2'

                  ds['precipitation'] = (('time','rgi_id'), precip)
                  ds['precipitation'].attrs['description'] = 'Total monthly precipitation'
                  ds['precipitation'].attrs['units'] = 'kg m 2'




                  num = num+1
                  break


       else:
        break
#Create a dataframe with the RGI ID and the lat lon Coordinates.
#Now process the DataFrame
loc_data = np.stack((rgi_id,latitude,longitude),axis=1)
loc_df = pd.DataFrame(loc_data,columns=['rgi_id','lat','lon'])
loc_df["lat"] = pd.to_numeric(loc_df["lat"], downcast="float")
loc_df["lon"] = pd.to_numeric(loc_df["lon"], downcast="float")
#Cut this up for binning
lat_cut = pd.cut(loc_df.lat, np.linspace(16, 57.25, 166),labels=(np.linspace(16, 57, 165)))
lon_cut = pd.cut(loc_df.lon, np.linspace(72, 143.25, 286),labels=(np.linspace(72, 143, 285)))
loc_df['lat_bin'] = lat_cut
loc_df['lon_bin'] = lon_cut
#Now output the File
loc_df.to_csv(working_dir+'/loc_bins_labeled.csv')
#df = ds.to_dask_dataframe()
#ds.to_netcdf(path=working_dir+'/test.nc',mode='w',format='NETCDF4')







