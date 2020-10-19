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
filename =
#filename = 'model_diagnostics_commitment.nc'

rgi_id = []
num=1
#Get list of RGI-IDs t0 use
for root,dirs,files in os.walk(raw_data_folder, topdown=False):
    for name in files:
        if name == filename:
            path = os.path.join(root,name)
            rgi = os.path.basename(os.path.dirname(path))
            rgi = str(rgi)
            rgi_id.append(rgi)
#Hacky way to start the process by making a file to input the data scraped from the loop into.
for root, dirs, files in os.walk(raw_data_folder, topdown=False):
    for name in files:
       if num == 1:
           if name == filename:
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
                  length = np.zeros(shape)
                  ela = np.zeros(shape)
                  # These are not
                  calving_m3 = None
                  calving_rate_myr = None
                  volume_bsl_m3 = None
                  volume_bwl_m3 = None
                  for i in range (len(rgi_id)):
                                 #Do some trickery to get the paths and speed it up over loops
                                  fold = rgi_id[i][:8]
                                  subfold = rgi_id[i][:11]
                                  folpath = raw_data_folder+'/'+fold+'/'+subfold+'/'+rgi_id[i]+'/model_diagnostics_commitment.nc'
                                  print('found glacier number: {}'.format(i))
                                  with xr.open_dataset(folpath) as ds_diag:
                                      vol[:, i] = ds_diag.volume_m3.values
                                      area[:, i] = ds_diag.area_m2.values
                                      length[:, i] = ds_diag.length_m.values
                                      ela[:, i] = ds_diag.ela_m.values
                                      if 'calving_m3' in ds_diag:
                                          if calving_m3 is None:
                                              calving_m3 = np.zeros(shape) * np.NaN
                                          calving_m3[:, i] = ds_diag.calving_m3.values
                                      if 'calving_rate_myr' in ds_diag:
                                          if calving_rate_myr is None:
                                              calving_rate_myr = np.zeros(shape) * np.NaN
                                          calving_rate_myr[:, i] = ds_diag.calving_rate_myr.values
                                      if 'volume_bsl_m3' in ds_diag:
                                          if volume_bsl_m3 is None:
                                              volume_bsl_m3 = np.zeros(shape) * np.NaN
                                          volume_bsl_m3[:, i] = ds_diag.volume_bsl_m3.values
                                      if 'volume_bwl_m3' in ds_diag:
                                          if volume_bwl_m3 is None:
                                              volume_bwl_m3 = np.zeros(shape) * np.NaN
                                          volume_bwl_m3[:, i] = ds_diag.volume_bwl_m3.values

                  ds['volume'] = (('time', 'rgi_id'), vol)
                  ds['volume'].attrs['description'] = 'Total glacier volume'
                  ds['volume'].attrs['units'] = 'm 3'
                  ds['area'] = (('time', 'rgi_id'), area)
                  ds['area'].attrs['description'] = 'Total glacier area'
                  ds['area'].attrs['units'] = 'm 2'
                  ds['length'] = (('time', 'rgi_id'), length)
                  ds['length'].attrs['description'] = 'Glacier length'
                  ds['length'].attrs['units'] = 'm'
                  ds['ela'] = (('time', 'rgi_id'), ela)
                  ds['ela'].attrs['description'] = 'Glacier Equilibrium Line Altitude (ELA)'
                  ds['ela'].attrs['units'] = 'm a.s.l'
                  if calving_m3 is not None:
                          ds['calving'] = (('time', 'rgi_id'), calving_m3)
                          ds['calving'].attrs['description'] = ('Total calving volume since '
                                                                'simulation start')
                          ds['calving'].attrs['units'] = 'm3'
                  if calving_rate_myr is not None:
                          ds['calving_rate'] = (('time', 'rgi_id'), calving_rate_myr)
                          ds['calving_rate'].attrs['description'] = 'Instantaneous calving rate'
                          ds['calving_rate'].attrs['units'] = 'm yr-1'
                  if volume_bsl_m3 is not None:
                          ds['volume_bsl'] = (('time', 'rgi_id'), volume_bsl_m3)
                          ds['volume_bsl'].attrs['description'] = ('Total glacier volume below '
                                                                   'sea level')
                          ds['volume_bsl'].attrs['units'] = 'm3'
                  if volume_bwl_m3 is not None:
                          ds['volume_bwl'] = (('time', 'rgi_id'), volume_bwl_m3)
                          ds['volume_bwl'].attrs['description'] = ('Total glacier volume below '
                                                                   'water level')
                          ds['volume_bwl'].attrs['units'] = 'm3'


                  num = num+1
                  break


       else:
        break
#Now output the File
ds.to_netcdf(path=working_dir+'/test.nc')






