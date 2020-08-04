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

                  lat = np.zeros(shape)
                  lon = np.zeros(shape)

                  for i in range (len(rgi_id)):
                                 #Do some trickery to get the paths and speed it up over loops
                                  fold = rgi_id[i][:8]
                                  subfold = rgi_id[i][:11]
                                  folpath = raw_data_folder+'/'+fold+'/'+subfold+'/'+rgi_id[i]+'/model_diagnostics_commitment.nc'
                                  print('found glacier number: {}'.format(i))
                                  with xr.open_dataset(folpath) as ds_diag:
                                      vol[:, i] = ds_diag.volume_m3.values
                                      area[:, i] = ds_diag.area_m2.values



                  ds['volume'] = (('time', 'rgi_id'), vol)
                  ds['volume'].attrs['description'] = 'Total glacier volume'
                  ds['volume'].attrs['units'] = 'm 3'
                  ds['area'] = (('time', 'rgi_id'), area)
                  ds['area'].attrs['description'] = 'Total glacier area'
                  ds['area'].attrs['units'] = 'm 2'

                  num = num+1
                  break


       else:
        break
#Now output the File
ds.to_netcdf(path=working_dir+'/test.nc')






